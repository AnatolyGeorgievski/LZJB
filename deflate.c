/*  Copyright 2021 Geolab. All Rights Reserved.
    Author: Anatoly Georgievskii <anatoly.georgievski@gmail.com>
*/
/*! \brief алгоритм распаковки Deflate Дифлятор
    Должен попасть в PNG, GZIP и поставить реккорд по скорости распаковки

    \see [RFC 1951] DEFLATE Compressed Data Format Specification version 1.3

    \sa ALISTAIR MOFFAT, Huffman Coding, ACM Comput. Surv., Vol. 1, No. 1, Article 1. Publication date: June 2019.

    101 xxxxx 7 - 1
    0xx x xxxx

    \sa https://github.com/google/zopfli/blob/master/src/zopfli/deflate.c
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "huffman.h"

#define DBG 0

// $ gcc   -dM -E -x c /dev/null -std=c11
#if __BYTE_ORDER__==__ORDER_LITTLE_ENDIAN__
#define BE16(x) __builtin_bswap16(x)
#define LE16(x) (x)
#else
#define LE16(x) __builtin_bswap16(x)
#define BE16(x) (x)
#endif // __BYTE_ORDER__
/*! При фиксированном методе кодирования используются длины и дистанции,
    к ним дополнительные биты для увеличения длинн */
static const struct {
    uint8_t length;
    uint8_t extra ;
} fixed_length[32] ={{0,0}, {1,0}, {2,0}, {3,0}, {4,0}, {5,0}, {6,0}, {7,0}, {8,1}, {10,1}, {12,1}, {14,1},
    {16,2}, {20,2}, {24,2}, {28,2}, {32,3}, {40,3}, {48,3}, {56,3}, {64,4}, {80,4}, {96,4}, {112,4},
    {128,5}, {160,5}, {192,5}, {224,5}, {255,0}, };
static const struct {
    uint16_t offset;
    uint8_t extra ;
} fixed_distance[30] = {{0,0}, {1,0}, {2,0}, {3,0}, {4,1}, {6,1}, {8,2}, {12,2},
    {16,3}, {24,3}, {32,4}, {48,4}, {64,5}, {96,5}, {128,6}, {192,6},
    {256,7}, {384,7}, {512,8}, {768,8}, {1024,9}, {1536,9}, {2048,10}, {3072,10},
    {4096,11}, {6144,11}, {8192,12}, {12288,12}, {16384,13}, {24576,13},};
static const uint8_t cl_order[20] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};

struct _tree_data {
    uint16_t Code;
//    uint8_t Len;
};

#define MAX_BITS 15
/*!
    \sa https://github.com/zlib-ng/zlib-ng/blob/develop/deflate.c
 */
static void gen_codes(uint8_t *tree_len, uint16_t *codes, int t_len, uint8_t *bl_count)
{// построение таблицы алфавита из max_code
    uint16_t code = 0;
    uint16_t next_code[MAX_BITS+1];
    next_code[0]=0;
    bl_count [0]=0;
    int i;
    for (i = 1; i < MAX_BITS; i++) {
        next_code[i] = code;
        code = (code + bl_count[i]) << 1;
    }
    next_code[i] = code;

    for (i=0; i< t_len; i++) {
        int len = tree_len[i];
        if (len != 0) {
                // tree[n].Code = bi_reverse(next_code[len]++, len);
            codes[i] = next_code[len]++;
        }
    }
}

static int btree_min_blen(const uint8_t *bl_count)
{
    int i=1;
    while(bl_count[i]==0) i++;
/*    for(i = 1; i <= MAX_BITS; i++) {
        if (bl_count[i]!=0) break;
    }*/
    return i;
}
/*! возвращает максимальную длину кодов */
static int btree_max_blen(const uint8_t *bl_count, int min_bl, int max_bl)
{
//    int max_blen=min_bl;
    int i=max_bl;
    while (bl_count[i]==0) i--;
/*    for(i = MAX_BITS; i >= 0; --i) {
        if (bl_count[i]!=0) break;
    }*/
    return i;
}
/*! возвращает число не нулевых символов */
static uint32_t btree_size(uint8_t *bl_count, int min_bl, int max_bl)
{
    uint32_t bt_size=0;
    int i;
    for(i = min_bl; i <= max_bl; i++) {
        bt_size += bl_count[i];
    }
    return bt_size;
}
/* Huffman code lookup table entry--this entry is four bytes for machines
   that have 16-bit pointers (e.g. PC's in the small or medium model).
   Valid extra bits are 0..13.  e == 15 is EOB (end of block), e == 16
   means that v is a literal, 16 < e < 32 means that v is a pointer to
   the next table, which codes e - 16 bits, and lastly e == 99 indicates
   an unused code.  If a code with e == 99 is looked up, this implies an
   error in the data. */

typedef struct _btree_codes btree_codes_t;
struct _btree_codes {
    int min_dl, max_dl, tabL;
    uint16_t *first_symbol; // размер таблицы min_dl
    uint16_t *first_code;   // размер таблицы min_dl
    uint8_t  *search_start; // где искать коды с данным префиксом, размер таблицы 1<<min_dl
    uint8_t  *search_end  ; // где искать коды с данным префиксом, размер таблицы 1<<min_dl
};
/*! генерирует таблицы неободимые для декодирования */
//static
void btree_codes(uint8_t *bl_count, int min_bl,int max_bl,
                            struct _btree_codes *bt, int tab)
{
    const int L=MAX_BITS;
    uint32_t bt_code=0, bt_offset=0,code;
    int t = tab - min_bl;
    int i;
    for(i = 0; i <= (max_bl); i++) {
        bt->first_symbol[i] = bt_offset;// номер символа заданной длины (i)
        bt->first_code  [i] = bt_code<<(L- min_bl -i); // первый код данной длины (i)
        uint32_t count = bl_count[i];
        if (count) {
            for (code = bt_code<<(t)>>(i); code<=((((bt_code+count)<<(t))-1)>>(i)); code++) {
                if (bt->search_start[code]==0) {
                    bt->search_start[code] = i;
                }
                bt->search_end[code] = i;
            }
        }
        bt_offset += count;
        bt_code = (bt_code+count)<<1;// расчитываем коды налету
    }
    bt->first_code  [i] = ~0;// не понятно
/*    for(i = 0; i < (1<<tab); i++) {
        if (bt->search_start[i]<=tab && bt->search_start[i]==bt->search_end[i])
            bt->decode_fast[i] = bt->first_symbol[i] + ((code - bt->first_code[i])>>(L-min_bl-bt->search_start[i]));
        else
            bt->decode_fast[i] = bt->search_end[i]<<12;
    }*/

    if (0) {
        printf("CODES: %d, <%d,%d\n", min_bl, max_bl, t);
        for(i = 0; i < (1<<tab); i++) {
            if (bt->search_start[i] != bt->search_end[i])
            printf("%2d SS:%d SE:%d\n", i, bt->search_start[i]+min_bl, bt->search_end[i]+min_bl);
        }
        printf("\n");
    }
}
/*! выполняет сортировку массива по коду */
static uint16_t* btree_code_order(uint16_t *alphabet, int alpha_size, uint8_t *bl_count, uint8_t *code_lengths, int cl_count, int min_bl, int max_bl)
{
    //int min_bl=1;
    uint32_t bt_offset=0;
    uint16_t bt_offs[max_bl-min_bl+1];
    int i;
    for(i = min_bl; i <= max_bl; i++) {
        bt_offs[i-min_bl] = bt_offset;
        bt_offset += bl_count[i];
    }
    for (i=0; i<cl_count; i++) {
        int len = code_lengths[i];
        if (len != 0) {
            alphabet[bt_offs[len-min_bl]++]=i;
        }
    }
    return alphabet;
}

#include <intrin.h>
inline
unsigned int __shrd(unsigned int into, unsigned int from, unsigned int c)
{
   unsigned int res;

   if (__builtin_constant_p(into) &&
       __builtin_constant_p(from) &&
       __builtin_constant_p(c))
   {
      res = (into >> c) | (from << (32 - c));
   }
   else
   {
      __asm ("shrd %b3, %2, %0"
          : "=rm" (res)
          : "0" (into), "r" (from), "ic" (c)
          : "cc");
   }

   return res;
}
static inline
uint64_t __shld(uint64_t res, uint64_t from, unsigned int c)
{
   if (__builtin_constant_p(res) &&
       __builtin_constant_p(from) &&
       __builtin_constant_p(c))
   {
      res = (res << c) | (from >> (64 - c));
   }
   else
   {
      __asm ("shld %b3, %2, %0"
          : "=rm" (res)
          : "0" (res), "r" (from), "ic" (c)
          : "cc");
   }

   return res;
}
static inline uint8_t* __attribute__((__target__("sse2")))
rle_memmove_sse2(uint8_t* dst, uint8_t* s, size_t mlen)
{
    int i;
    if (dst-s>=4){// перекрытие
        if (dst-s>=16) {
            for (i=0; i< (mlen>>4); i++) {
                __m128i v = _mm_loadu_si128((void*)s); s+=16;
                _mm_storeu_si128((void*)dst, v); dst+=16;
            }
            mlen &=15;
        }
        for(i=0; i<(mlen>>2); i++) {
            *(uint32_t*)dst = *(uint32_t*)s;
            dst+=4, s+=4;
        }
        mlen &=3;
    }
    for(i=0; i<(mlen); i++)
        *dst++ = *s++;
    return dst;
}
// Вложенный цикл практически не реализуется потому что последовательности меньше 32
static inline uint8_t* __attribute__((__target__("avx")))
rle_memmove_avx(uint8_t* dst, uint8_t* s, size_t mlen)
{
    int i;
    if (dst-s>=8){// перекрытие
        if (dst-s>=32) {
            for (i=0; i< (mlen>>5); i++) {
                __m256i v = _mm256_loadu_si256((void*)s); s+=32;
                _mm256_storeu_si256((void*)dst, v); dst+=32;
            }
            mlen &=31;
        }
        for(i=0; i<(mlen>>3); i++) {
            *(uint64_t*)dst = *(uint64_t*)s;
            dst+=8, s+=8;
        }
        mlen &=7;
    }
    for(i=0; i<(mlen); i++)
        *dst++ = *s++;
    return dst;
}
// длины последовательностей обычно маленькие, этот вариант вероятно оптимальный
static uint8_t*
__attribute__((optimize("Os","no-tree-loop-distribute-patterns")))
rle_memmove_x64(uint8_t* dst, uint8_t* s, size_t mlen)
{
    if (dst-s>=8){// перекрытие
        int blocks = mlen>>3;
        if (blocks) do {
            *(uint64_t*)dst = *(uint64_t*)s;
            dst+=8, s+=8;
        } while(--blocks);
        mlen&=7;
    }
    if (mlen) do{
        *dst++ = *s++;
    } while (--mlen);
    return dst;
}
#if 0
#if defined(__AVX512F__)
#define rle_memmove rle_memmove_x64
#elif 0//defined(__AVX__)
#define rle_memmove rle_memmove_avx
#elif defined(__SSE2__)
#define rle_memmove rle_memmove_sse2
#else
#define rle_memmove rle_memmove_x64
#endif // 1
#endif // 0

#define rle_memmove rle_memmove_x64
static inline uint16_t __attribute__((__target__("avx512bitalg")))
reverse_16bits(uint64_t stream)
{
    __m128i  bits  = _mm_set_epi64x(0x0001020304050607ULL, 0x0001020304050607ULL+0x0808080808080808ULL);
    __mmask16 mask = _mm_bitshuffle_epi64_mask(_mm_set1_epi64x(stream), bits);
    return mask;
}

/* TODO держать вектор reverse_stream 128 бит из, использовать аффинные преобразования */
// \see https://graphics.stanford.edu/~seander/bithacks.html#BitReverseTable
const uint8_t BitReverseTable16l[16]={0,8,4,0xC,2,0xA,6,0xE,1,9,5,0xD,3,0xB,7,0xF};
const uint8_t BitReverseTable16h[16]={0,0x80,0x40,0xC0,0x20,0xA0,0x60,0xE0,0x10,0x90,0x50,0xD0,0x30,0xB0,0x70,0xF0};
static inline uint8_t revbits_(uint8_t b){
    return BitReverseTable16l[b>>4] | (BitReverseTable16h[b&0xF]);
}
#if defined(__ARM__)
// \see Arm C Language Extensions Documentation, Release ACLE Q3 2020
static inline uint8_t revbits_arm(uint8_t b)
{
    return __builtin_rbit(b);
}
// uint32_t __rbit(uint32_t x);
#endif // defined
static inline uint8_t revbits(uint8_t b){
    // B=1023: 0x80200803 nd=41
    // b = (b * 0x0202020202ULL & 0x010884422010ULL) % 1023;
    // x =
    // b = b - x * 1024 + x;
    return ((b * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
}



static deflate_t context = {0};
uint8_t * deflate (uint8_t *dst, uint8_t *src, size_t s_len, deflate_t* ctx)
{
    if (ctx==NULL) {
        ctx = &context;
        ctx->n_bits =0;
    }
    uint32_t stream = ctx->stream;
    uint32_t n_bits = ctx->n_bits;
    // nested functions
    inline void     stream_preload (int n){
        if (n > n_bits) {
            stream |= (uint32_t)(*src++)<<n_bits;
            n_bits += 8;
            if (n > n_bits) {
                stream |= (uint32_t)(*src++)<<n_bits;
                n_bits += 8;
            }
        }
    }
    inline uint32_t stream_read_bit(int n){
        stream_preload(n);
        n_bits -=n;

        uint32_t code = stream & ((1<<n)-1);// инструкция BZHI
        stream>>=n;
        return code;
    }
   inline uint32_t huffman_read_bit(uint32_t code, int n){
        n_bits -= n;
        do{// это может быть одна инструкция через флаг переноса CF. или SHLD на ARM это может быть RBIT
            code = (code<<1) + (stream &1);
            stream>>=1;
        } while(--n);
        return code;
    }
    /* Декодирование с ипсользованием алфавита
    ACM Comput. Surv., Vol. 1, No. 1, Article 1. Publication date: June 2019

    identify ℓ such that first_code_l[ℓ] ≤ code < first_code_l[ℓ + 1];
    set offset ← (code − first_code_l[ℓ]) >> (L − ℓ);

    set s ← first_symbol[ℓ] + offset;
    set code ← ((code << ℓ) & maskL) + getbits(ℓ);
    output s;

    set ℓ ← search_start[code >> (L − t)];
    while code ≥ first_code_l[ℓ + 1] do
        set ℓ ← ℓ + 1.

    */
    void *bsearch(const void *key, const void *base, size_t num, size_t size,
              int (*cmp)(const void *key, const void *))
    {
            size_t l = 0, u = num;
            while (l < u) {
                    register const size_t mid = (l + u)>>1;
                    register const char* p = (const char*)base + mid * size;
                    register int result = cmp(key, p);
                    if (result < 0)
                            u = mid;
                    else if (result > 0)
                            l = mid + 1;
                    else
                            return (void *)p;
            }
            return NULL;
    }

    /*! \brief Декодирование табличный метод */
    inline uint32_t huffman_decode2(btree_codes_t* bt, int min_dl, int max_dl, int tab)
    {
        const int L=MAX_BITS;
        stream_preload(max_dl);
        uint16_t code  = ((uint16_t)revbits(stream) <<8)|((uint16_t)revbits(stream>>8));// надо быстро развернуть биты _rbits16
        int l = bt->search_start[code>> (16-tab)];
        int u = bt->search_end  [code>> (16-tab)];
        code>>=(16 - L);
        uint16_t *fc = &bt->first_code[0];
        if (l<u) while(code >= fc[l+1]) l++;
        stream>>=l+min_dl; n_bits-=l+min_dl;
        return bt->first_symbol[l] + ((code - fc[l])>>(L-(l+min_dl)));
    }
    /*! \brief Декодирование
        \param min_cwl - минимальное число бит,
     */
    uint32_t huffman_decode(uint8_t* cwl_count, int max_cwl)
    {
        uint32_t bt_code=0, bt_offset=0, code=0;
        stream_preload(max_cwl);
        //__asm volatile("# LLVM-MCA-BEGIN decode");
        int i=1;
        register long count;// = cwl_count[i];
        for(;; i+=1, bt_code <<= 1) {// Алгоритм поиска(cwl_count, max_cwl), заменить на while, убрать проверку
            code = (code<<1) + (stream &1);
            stream>>=1;
            count = cwl_count[i];// число кодов с разрядностью (i)
            if (count==0)  continue;
            // bt_code - первый символ в кодовой таблице с разрядностью (i)
            if (count> code - bt_code) break;
            bt_code   += count;
            bt_offset += count;
        }
        //__asm volatile("# LLVM-MCA-END decode");
        n_bits -= i;
        return bt_offset + (code - bt_code);
    }

    // Each block of compressed data begins with 3 header bits containing the following data:
    uint8_t * s_end = src + s_len;
    int bfinal= stream_read_bit(1);
    int btype = stream_read_bit(2);
    if (btype==0) {// 3.2.4. Non-compressed blocks (BTYPE=00)
        uint16_t len  = LE16(*(uint16_t *)src); src+=2;// LEN is the number of data bytes in the block.
        uint16_t nlen = LE16(*(uint16_t *)src); src+=2;// NLEN is the one's complement of LEN.
        __builtin_memcpy(dst, src, len);
        if(DBG)printf("Non-compressed block len=%d\n", (int)len);
        src+=len, dst+=len;
        stream = 0, n_bits=0;
    } else
    if (btype==1) {// 3.2.6. Compression with fixed Huffman codes (BTYPE=01)
        if(1)printf("Compression with fixed Huffman codes\n");
        uint32_t code, mlen, offset, extra;
        while(src<s_end) {
            // literal/length decode
            stream_preload(9);
            code = huffman_read_bit(0, 5);
            if (code<0b11000) {
                if (code<0b00110){
                    code = huffman_read_bit(code, 2);
                    code = code - 0b0000000 + 256;
                } else {
                    code = huffman_read_bit(code, 3);
                    code = code - 0b00110000 + 0;
                }
            } else {
                if (code<0b11001){
                    code = huffman_read_bit(code, 3);
                    code = code - 0b11000000 + 280;
                } else {
                    code = huffman_read_bit(code, 4);
                    code = code - 0b110010000 + 144;
                }
            }
#if 0 // Читабельная реализация
            switch (code & 0x1F) {// 5-bit prefix
            case 0b00000 ... 0b00101:// 256 - 279
                code = huffman_read_bit(code, 2);
                code = code - 0b0000000 + 256;
                break;
            case 0b00110 ... 0b10111://   0 - 143
                code = huffman_read_bit(code, 3);
                code = code - 0b00110000 + 0;
                break;
            case 0b11000 ... 0b11000:// 280 - 287
                code = huffman_read_bit(code, 3);
                code = code - 0b11000000 + 280;
                break;
            case 0b11001 ... 0b11111:// 144 - 255
                code = huffman_read_bit(code, 4);
                code = code - 0b110010000 + 144;
                break;
            }
#endif // 0
            if (code<256){// literal decode
                *dst++ = code;
                continue;
            }
            if (code==256) {
                break;
            }
            // length decode
            mlen  = fixed_length[code-257].length+3;// минимальная длина кодирования.
            extra = fixed_length[code-257].extra;
            if (extra){
                mlen += stream_read_bit(extra);
            }
            // distance decode
            stream_preload(5);
            code  = huffman_read_bit(0, 5);
            offset= fixed_distance[code].offset+1;// минимальное смещение
            extra = fixed_distance[code].extra;
            if (extra) {
                offset += stream_read_bit(extra);
            }
            {
                uint8_t * s = dst-offset;
                rle_memmove(dst,s,mlen);
                dst+=mlen;
            }
        }
    } else
    if (btype==2) {// 3.2.7. Compression with dynamic Huffman codes (BTYPE=10)
        // read representation of code trees
        uint8_t bl_count[MAX_BITS+1]={0};
        int i;

        uint16_t hlit=257+ stream_read_bit(5);//    5 Bits: HLIT,  # of Literal/Length codes - 257(257 - 286)
        uint8_t hdist= 1 + stream_read_bit(5);//    5 Bits: HDIST, # of Distance codes - 1        (1 - 32)
        uint8_t hclen= 4 + stream_read_bit(4);//    4 Bits: HCLEN, # of Code Length codes - 4     (4 - 19)

        uint8_t tree_len[20]={0};
        if(DBG)printf("hlit %2d, hdist=%2d, hclen=%2d\n", hlit, hdist, hclen);
        for (i=0;i<hclen;i++) {
            uint8_t code = stream_read_bit(3);
            bl_count[code] ++;
            tree_len[cl_order[i]] = code;
        }
        for (;i<19;i++) {
            tree_len[cl_order[i]] = 0;
        }
        // коды по порядку алфавита
        if(DBG) {
            uint16_t tree [20]={0};
            gen_codes(tree_len, tree, 19, bl_count);
            for(i=0; i<19;i++)
                printf("\t%2d: %2d, %x\n", i, tree_len[i], tree[i]);
        }
        // коды по порядку кодирования
        int min_bl = btree_min_blen(bl_count);
        int max_bl = btree_max_blen(bl_count, min_bl,7);
        uint32_t alpha_size = btree_size(bl_count, min_bl, max_bl);
        uint16_t alpha[alpha_size];
        btree_code_order(alpha, alpha_size, bl_count, tree_len, 19, min_bl, max_bl);
        if (DBG) {
            printf("  in order:");
            for(i=0; i<alpha_size;i++)
                printf(" %d", alpha[i]);
            printf("\n");
        }
        uint8_t cwl_count[MAX_BITS+1];
        uint8_t cwl_tree_len[hlit];
        //__builtin_bzero(cwl_tree_len, hlit);
        void huffman_read_alpha(uint8_t* tree, int tlen, uint8_t* tl_count)
        {
            __builtin_bzero(tl_count, MAX_BITS+1);
            uint32_t prev_code=0, code;
            uint8_t* e_tree = tree+tlen;
            while(tree<e_tree){
                code = huffman_decode(bl_count, max_bl);
                code = alpha[code];
                if (code<16) {
                    prev_code = code;
                    tl_count[prev_code] ++;
                    *tree++ = prev_code;
                }
                else {
                    int n;
                    if (code==16) {// Copy the previous code length 3 - 6 times.
                        n = 3 + stream_read_bit(2);
                        tl_count[prev_code] += n;
                    } else {
                        if (code==17){// Repeat a code length of 0 for 3 - 10 times.
                            n = 3 + stream_read_bit(3);
                        } else {// Repeat a code length of 0 for 11 - 138 times
                            n = 11 + stream_read_bit(7);
                        }
                        prev_code = 0;
                    }
                    do {
                        *tree++ = prev_code;
                    } while(--n);
                }
            }
        }
        huffman_read_alpha(cwl_tree_len, hlit, cwl_count);
        if (0) printf("max_bl = %d\n", max_bl);
                // коды алфавита по порядку кодирования
        int min_cwl = btree_min_blen(cwl_count);
        int max_cwl = btree_max_blen(cwl_count,min_cwl,MAX_BITS);
        uint32_t cwl_size = btree_size(cwl_count, min_cwl, max_cwl);
        uint16_t cwl_alpha[cwl_size];// содержит симовлы в порядке кодирования
        btree_code_order(cwl_alpha, cwl_size, cwl_count, cwl_tree_len, hlit, min_cwl, max_cwl);
#if 1
        uint16_t cwl_first_symbol[max_cwl - min_cwl+1];
        uint16_t cwl_first_code  [max_cwl - min_cwl+2];
        //const int cwl_T = 8;
        int tab_cwl = 8;
        // tab_cwl = btree_bound(8, min_cwl, max_cwl);
        if (tab_cwl< min_cwl) tab_cwl = min_cwl;
        else if (tab_cwl> max_cwl) tab_cwl = max_cwl;
        uint8_t  cwl_search_start[1<<tab_cwl];
        uint8_t  cwl_search_end  [1<<tab_cwl];
        uint8_t  cwl_decode_fast [1<<tab_cwl];
        __builtin_bzero(cwl_search_start, 1<<tab_cwl);
        //__builtin_bzero(cwl_search_end,   1<<tab_cwl);
        btree_codes_t cwl_btree = {.search_start = cwl_search_start,.search_end = cwl_search_end, .tabL = 0,
            .first_symbol = cwl_first_symbol, .first_code = cwl_first_code};
        btree_codes(cwl_count+min_cwl, min_cwl, max_cwl-min_cwl, &cwl_btree, tab_cwl);
#endif // 0
        if (0){
            printf("\n  in order:");
            for(i=0; i<cwl_size;i++)
                printf(" %X", cwl_alpha[i]);
            printf("\n");
            uint16_t cwl_codes[hlit+1];
            __builtin_bzero(cwl_codes, hlit);
            gen_codes(cwl_tree_len, cwl_codes, hlit, cwl_count);
            for(i=0; i<=hlit;i++){
                printf("\t%2d: %2d, %x\n", i, cwl_tree_len[i], cwl_codes[i]);
            }
        }

// коды для алфавита смещений
        uint8_t dl_count[MAX_BITS+1]={0};
        uint8_t dl_tree_len[hdist];
        //__builtin_bzero(dl_tree_len, hdist);
        huffman_read_alpha(dl_tree_len, hdist, dl_count);
#if 0
        tree = dl_tree_len;
        e_tree = tree+hdist;
        while(tree<e_tree){
            code = huffman_decode(bl_count, max_bl);
            code = alpha[code];
            if (code<16) {
                prev_code = code;
                dl_count[prev_code] ++;
                *tree++ = prev_code;
            }
            else {
                if (code==16) {// Copy the previous code length 3 - 6 times.
                    n = 3 + stream_read_bit(2);
                    dl_count[prev_code] += n;
                } else {
                    if (code==17){// Repeat a code length of 0 for 3 - 10 times.
                        n = 3 + stream_read_bit(3);
                    } else {// Repeat a code length of 0 for 11 - 138 times
                        n = 11 + stream_read_bit(7);
                    }
                    prev_code = 0;
                }
                int count = n;
                do {
                    *tree++ = prev_code;
                } while(--count);
            }
        }
#endif
                // коды дистанций по порядку кодирования
        int min_dl = btree_min_blen(dl_count);
        int max_dl = btree_max_blen(dl_count, min_dl, MAX_BITS);
        uint32_t dl_size = btree_size(dl_count, min_dl, max_dl);
        uint16_t dl_alpha[dl_size];
        btree_code_order(dl_alpha, dl_size, dl_count, dl_tree_len, hdist, min_dl, max_dl);
#if 1
        uint16_t dl_first_symbol[max_dl - min_dl+1];
        uint16_t dl_first_code  [max_dl - min_dl+2];
        //const int dl_T = 5;
        int tab_dl = 4;
        if (min_dl>tab_dl) tab_dl = min_dl;
        else
            if (max_dl<tab_dl) tab_dl = max_dl;
        uint8_t  dl_search_start[1<<tab_dl];
        uint8_t  dl_search_end  [1<<tab_dl];
        __builtin_bzero(dl_search_start, 1<<tab_dl);
        //__builtin_bzero(dl_search_end,   1<<tab_dl);
        btree_codes_t dl_btree = {.search_start = dl_search_start,.search_end = dl_search_end,
            .first_symbol = dl_first_symbol, .first_code = dl_first_code};
        btree_codes(dl_count+min_dl, min_dl, max_dl-min_dl, &dl_btree, tab_dl);
#endif // 0
        if(0) {// отладка
            printf("\n  in order:");
            for(i=0; i<dl_size;i++)
                printf(" %d:%d(%d)", i, dl_alpha[i], dl_tree_len[dl_alpha[i]]);
            printf("\n");
            uint16_t dl_codes[hdist+1];
            __builtin_bzero(dl_codes, (hdist+1)*2);
            gen_codes(dl_tree_len, dl_codes, hdist, dl_count);
            char code[MAX_BITS+1];
            for(i=0; i<hdist;i++){
                int j;
                for(j=0;j<dl_tree_len[i]; j++)
                    code[j] = dl_codes[i]&(1<<(dl_tree_len[i]-j-1))?'1':'0';

                printf("\t%2d: %2d, %-.*s\n", i, dl_tree_len[i], dl_tree_len[i], code);
            }
        }
        //printf("max_cwl %2d, max_dl=%2d, max_bl=%2d\n", max_cwl, max_dl, max_bl);
        // выполнить само декодирование
        int mlen, extra;// match length, extra bits
        uint32_t offset;// distance offset
//        uint64_t ts=0, ts2=0;
        while(src<s_end){
            //ts-= __builtin_ia32_rdtsc();
            uint32_t code;
            code = huffman_decode2(&cwl_btree, min_cwl, max_cwl, tab_cwl);
            //code = huffman_decode(cwl_count, max_cwl);
            code = cwl_alpha[code];
            //ts+= __builtin_ia32_rdtsc();
            if (code<256){
                *dst++ = code;
            } else {
                if (code==256){
                    if(DBG)printf("EOB\n");
                    break;
                }
                mlen = fixed_length[code-257].length + 3;
                extra= fixed_length[code-257].extra;
                if (extra){
                    mlen += stream_read_bit(extra);
                }
                // distance decode
                code = huffman_decode2(&dl_btree, min_dl, max_dl, tab_dl);
                //code = huffman_decode(dl_count, max_dl);
                //if (0 && code1 != code) _Exit(0xC0DE);
                code = dl_alpha[code];// distance length alphabet
                offset= fixed_distance[code].offset+1;
                extra = fixed_distance[code].extra;
                if (extra){
                    offset += stream_read_bit(extra);
                }
                // статический выбор стратегии: CPU::target, ...
                //if (dst+mlen>d_end)
                rle_memmove(dst, dst-offset, mlen);
                dst+=mlen;
            }
        }
        //if(1) printf("ts1=%lld ts2=%lld %d %d\n", ts, ts2, max_cwl, max_dl);
    }
    ctx->s_end  = src;
    ctx->stream = stream;
    ctx->n_bits = n_bits;
    return dst;
}

#if defined(TEST_DEFLATE)
int main()
{
    uint8_t buf[2048];
    uint8_t * dst = buf;
    uint8_t test[] =  "\x73\x49\x4D\xCB\x49\x2C\x49\x55\x00\x11\x00";
    uint8_t tst2[] =  {
        0b00001100, 0b11001000, 0b01000001, 0b00001010, 0b10000000, 0b00100000, 0b00010000, 0b00000101,
        0b11010000, 0b01111101, 0b11010000, 0b00011101, 0b11111110, 0b00001001, 0b10111010, 0b10000100,
        0b11101011, 0b10100000, 0b00101011, 0b01001100, 0b11111010, 0b10110101, 0b00000001, 0b00011101,
        0b00100001, 0b00100111, 0b10100001, 0b11011011, 0b11010111, 0b01011011, 0b10111110, 0b11010000,

        0b10101101, 0b11011100, 0b11100010, 0b01001111, 0b00010101, 0b11010111, 0b01101110, 0b00000011,
        0b11011101, 0b01110000, 0b00110010, 0b11110110, 0b10100110, 0b01010110, 0b00100000, 0b10000110,
        0b00111101, 0b00011100, 0b00011011, 0b10001110, 0b01001010, 0b00011001, 0b11111100, 0b00011111,
        0b10010010, 0b10100110, 0b00001110, 0b00100110, 0b11111000, 0b00100101, 0b00001110, 0b11100110,

        0b11001100, 0b11101000, 0b00111010, 0b00001001, 0b01101101, 0b10001101, 0b01001001, 0b11000101,
        0b01011001, 0b11011111, 0b01110101, 0b11111001, 0b00000110, 0b00000000,
        };
    int d_len = deflate(dst, test, sizeof(test)-1, NULL) - dst;
    printf("dlen=%d: %s\n", d_len, dst);
    d_len = deflate(dst, tst2, sizeof(tst2), NULL) - dst;
    printf("dlen=%d: %s\n", d_len, dst);
    return 0;
}
#endif
