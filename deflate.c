/*! \brief алгоритм распаковки Deflate Дифлятор

    [RFC 1951] DEFLATE Compressed Data Format Specification version 1.3
    Должен попасть в PNG

    Кодирование тесктов: 6 бит на каждый символ 0xxxxx -
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define BE16(x) __builtin_bswap16(x)
/*! При фиксированном методе кодирования используются длины и дистанции,
    к ним дополнительные биты для увеличения длинн */
static const uint8_t fixed_length_mlen[32] =
    {0,1,2,3,4,5,6,7,8,10,12,14,16,20,24,28,32,40,48,56,64,80,96,112,128,160,192,224,255};
static const uint8_t fixed_length_extra[32]=
    {0,0,0,0,0,0,0,0,1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4,  4,  5,  5,  5,  5,  0};
static const uint16_t fixed_distance_offset[32]=
    {0,1, 2,3, 4,6, 8,12, 16,24, 32,48, 64,96, 128,192, 256,384, 512,768, 1024,1536,
    2048,3072, 4096,6144, 8192,12288, 16384,24576};
static const uint8_t fixed_distance_extra[32]=
    {0,0, 0,0, 1,1, 2,2, 3,3, 4,4, 5,5, 6,6, 7,7, 8,8, 9,9, 10,10, 11,11, 12,12, 13,13};
static const uint8_t cl_order[20] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};

struct _tree_data {
    uint16_t Code;
//    uint8_t Len;
};

#define MAX_BITS 15
/*!
    \sa https://github.com/zlib-ng/zlib-ng/blob/develop/deflate.c
 */
static void gen_codes(uint8_t *tree_len, struct _tree_data *tree, int max_code, uint8_t *bl_count)
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

    for (i=0; i<=max_code; i++) {
        int len = tree_len[i];
        if (len != 0) {
                // tree[n].Code = bi_reverse(next_code[len]++, len);
            tree[i].Code = next_code[len]++;
        }
    }
}

int btree_max_blen(uint8_t *bl_count)
{
    int max_blen=0;
    int i;
    for(i = 1; i <= MAX_BITS; i++) {
        if (bl_count[i]!=0) max_blen = i;
    }
    return max_blen;
}
static uint32_t btree_size(uint8_t *bl_count)
{
    uint32_t bt_offset=0;
    int i;
    for(i = 1; i <= MAX_BITS; i++) {
        bt_offset += bl_count[i];
    }
    return bt_offset;
}
/*! выполняет сортировку массива по коду */
uint16_t* btree_code_order(uint16_t *alphabet, int alpha_size, uint8_t *bl_count, uint8_t *code_lengths, int cl_count)
{
    uint32_t bt_offset=0;
    uint16_t bt_offs[MAX_BITS+1];
    int i;
    for(i = 1; i <= MAX_BITS; i++) {
        bt_offs[i] = bt_offset;
        bt_offset += bl_count[i];
    }
    for (i=0; i<cl_count; i++) {
        int len = code_lengths[i];
        if (len != 0) {
            alphabet[bt_offs[len]++]=i;
        }
    }
    return alphabet;
}

#include <intrin.h>
uint8_t huffman_read_test(uint32_t code, long* src, int offset, int n){
//    code = (code+code) + (src[offset>>3] >> (offset&7))&1;
    int i;
    for(i=0;i<n;i++) {
        code = (code<<1) + _bittest(src, offset++);
    }
    return code;
}
uint8_t huffman_read_test2(uint32_t code, long* src, int offset){
    code = (code<<1) | _bittest64 (src, offset);//, code, code, &code);
    return code;
}
uint8_t * deflate (uint8_t *dst, uint8_t *src, int s_len)
{
    uint32_t stream = 0;
    uint8_t n_bits=0;
    // nested functions
    uint32_t stream_read_bit(int n){
        while (n > n_bits) {
            stream |= (uint32_t)(*src++)<<n_bits;
            n_bits += 8;
        }
        n_bits -=n;

        uint32_t code = stream & ((1<<n)-1);// инструкция BZHI
        stream>>=n;
        return code;
    }
    uint32_t huffman_read_bit(uint32_t code, int n){
        while (n > n_bits) {
            stream |= (uint32_t)(*src++)<<n_bits;
            n_bits += 8;
        }
        n_bits -= n;
        do{// это может быть одна инструкция через флаг переноса CF. или SHLD на ARM это может быть RBIT
            code = (code<<1) | (stream &1);
            stream>>=1;
        } while(--n);
        return code;
    }
    /* Декодирование с ипсользованием алфавита */
    uint32_t huffman_decode(uint8_t* cwl_count, int max_cwl)
    {
        uint32_t bt_code=0, bt_offset=0, code=0;
        int i;
        for(i=1; i<=max_cwl; i+=1, bt_code <<= 1) {// Алгоритм поиска(cwl_count, max_cwl), заменить на while, убрать проверку
            code = huffman_read_bit(code, 1);
            if (cwl_count[i]==0) continue;
            if (code>=bt_code && code - bt_code < cwl_count[i]) break;
            bt_code   += cwl_count[i];
            bt_offset += cwl_count[i];
        }
        return bt_offset + (code - bt_code);
    }

    // Each block of compressed data begins with 3 header bits containing the following data:
    uint8_t * s_end = src + s_len;
    int bfinal= stream_read_bit(1);
    int btype = stream_read_bit(2);
    if (btype==0) {// 3.2.4. Non-compressed blocks (BTYPE=00)
        uint16_t len  = (*(uint16_t *)src); src+=2;// LEN is the number of data bytes in the block.
        uint16_t nlen = (*(uint16_t *)src); src+=2;// NLEN is the one's complement of LEN.
        __builtin_memcpy(dst, src, len);
        printf("\tNon-compressed block len=%d\n", (int)len);
        src+=len, dst+=len;
    } else
    if (btype==1) {// 3.2.6. Compression with fixed Huffman codes (BTYPE=01)
        printf("Compression with fixed Huffman codes\n");
        uint32_t code, mlen, offset, extra;
        while(src<s_end) {
            // literal/length decode
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
            if (code==256) break;
            // length decode
            mlen  = fixed_length_mlen  [code-257]+3;// минимальная длина кодирования.
            extra = fixed_length_extra [code-257];
            if (extra){
                mlen += huffman_read_bit(0, extra);
            }
            // distance decode
            code  = huffman_read_bit(0, 5);
            offset= fixed_distance_offset[code]+1;// минимальное смещение
            extra = fixed_distance_extra [code];
            if (extra) {
                offset += huffman_read_bit(0, extra);
            }
            {
                uint8_t * s = dst-offset;
                int i;
                for(i=0; i<mlen; i++)
                    *dst++ = *s++;
            }
        }
    } else
    if (btype==2) {// 3.2.7. Compression with dynamic Huffman codes (BTYPE=10)
        // read representation of code trees
        uint8_t * d_start = dst;
        uint8_t bl_count[MAX_BITS+1]={0};
        int i, n;

        uint16_t hlit=257+ stream_read_bit(5);//    5 Bits: HLIT,  # of Literal/Length codes - 257(257 - 286)
        uint8_t hdist= 1 + stream_read_bit(5);//    5 Bits: HDIST, # of Distance codes - 1        (1 - 32)
        uint8_t hclen= 4 + stream_read_bit(4);//    4 Bits: HCLEN, # of Code Length codes - 4     (4 - 19)

        uint8_t tree_len[hclen];
        __builtin_bzero(tree_len, hclen);
        void code_length_repeat(uint8_t code, int n){
            bl_count[code] += n;
            do {
                tree_len[cl_order[i++]] = code;
            }while (--n);
        }

        printf("hlit %2d, hdist=%2d, hclen=%2d\n", hlit, hdist, hclen);
        for (i=0;i<hclen;) {// i++ внутри функции code_length_repeat
            uint8_t code = stream_read_bit(3);
            if (code<16) {// Represent code lengths of 0 - 15
                code_length_repeat(code, 1);
            } else {
                if (code==16) {// Copy the previous code length 3 - 6 times.
                    code = tree_len[cl_order[i-1]];
                    n = 3 + stream_read_bit(2);
                } else
                {
                    if (code==17)// Repeat a code length of 0 for 3 - 10 times.
                        n = 3 + stream_read_bit(3);
                    else // Repeat a code length of 0 for 11 - 138 times
                        n = 11 + stream_read_bit(7);
                    code = 0;
                }
                code_length_repeat(code, n);
            }
        }
        // коды по порядку алфавита
        if(1) {
            struct _tree_data tree [20]={0};
            gen_codes(tree_len, tree, 18, bl_count);
            for(i=0; i<19;i++)
                printf("\t%2d: %2d, %x\n", i, tree_len[i], tree[i].Code);
        }
        // коды по порядку кодирования
        int max_blen = btree_max_blen(bl_count);
        uint32_t alpha_size = btree_size(bl_count);
        uint16_t alpha[alpha_size];
        btree_code_order(alpha, alpha_size, bl_count, tree_len, 19);
        if (1) {
            printf("  in order:");
            for(i=0; i<alpha_size;i++)
                printf(" %d", alpha[i]);
            printf("\n");
        }
        uint32_t prev_code=0, code;
        uint8_t cwl_count[MAX_BITS+1]={0};
        uint8_t cwl_tree_len[hlit];
        __builtin_bzero(cwl_tree_len, hlit);
        int alpha_idx=0;
        while(alpha_idx<hlit){
            code = huffman_decode(bl_count, max_blen);
            code = alpha[code];
            if (code<16) {
                prev_code = code;
                cwl_count[prev_code] ++;
                cwl_tree_len[alpha_idx++] = prev_code;
            }
            else {
                if (code==16) {// Copy the previous code length 3 - 6 times.
                    n = 3 + stream_read_bit(2);
                    cwl_count[prev_code] += n;
                    int count = n;
                    do {
                        cwl_tree_len[alpha_idx++] = prev_code;
                    } while(--count);
                } else {
                    if (code==17){// Repeat a code length of 0 for 3 - 10 times.
                        n = 3 + stream_read_bit(3);
                    } else {// Repeat a code length of 0 for 11 - 138 times
                        n = 11 + stream_read_bit(7);
                    }
                    alpha_idx+=n;
                }
            }
        }
                // коды алфавита по порядку кодирования
        int max_cwl = btree_max_blen(cwl_count);
        uint32_t cwl_size = btree_size(cwl_count);
        uint16_t cwl_alpha[cwl_size];
        btree_code_order(cwl_alpha, cwl_size, cwl_count, cwl_tree_len, hlit);
        if (1){
            printf("\n  in order:");
            for(i=0; i<cwl_size;i++)
                printf(" %X", cwl_alpha[i]);
            printf("\n");
        }

// коды для алфавита смещений
        prev_code=0;
        uint8_t dl_count[MAX_BITS+1]={0};
        uint8_t *dl_tree_len=malloc(hdist+16);//[hdist];
        __builtin_bzero(dl_tree_len, hdist+16);
        alpha_idx=0;
        while(alpha_idx<hdist){
            code = huffman_decode(bl_count, max_blen);
            code = alpha[code];
            if (code<16) {
                prev_code = code;
                dl_count[prev_code] ++;
                dl_tree_len[alpha_idx++] = prev_code;
            }
            else {
                if (code==16) {// Copy the previous code length 3 - 6 times.
                    n = 3 + stream_read_bit(2);
                    dl_count[prev_code] += n;
                    int count = n;
                    do {
                        dl_tree_len[alpha_idx++] = prev_code;
                    } while(--count);
                } else {
                    if (code==17){// Repeat a code length of 0 for 3 - 10 times.
                        n = 3 + stream_read_bit(3);
                    } else {// Repeat a code length of 0 for 11 - 138 times
                        n = 11 + stream_read_bit(7);
                    }
                    alpha_idx+=n;
                }
            }
        }
                // коды дистанций по порядку кодирования
        int max_dl = btree_max_blen(dl_count);
        uint32_t dl_size = btree_size(dl_count);
        uint16_t dl_alpha[dl_size];
        btree_code_order(dl_alpha, dl_size, dl_count, dl_tree_len, hdist);
        if(1) {// отладка
            printf("\n  in order:");
            for(i=0; i<dl_size;i++)
                printf(" %X", dl_alpha[i]);
            printf("\n");
        }
        // выполнить само декодирование
        while(src<s_end){
            uint32_t code = huffman_decode(cwl_count, max_cwl);
            code = cwl_alpha[code];
            if (code<256){
                *dst++ = code;
            } else {
                if (code==256)
                    break;
                int mlen = fixed_length_mlen [code-257]+3;
                int extra= fixed_length_extra[code-257];
                if (extra){
                    mlen += stream_read_bit(extra);
                }
                // distance decode
                code = huffman_decode(dl_count, max_dl);
                code = dl_alpha[code];
                uint32_t offset= fixed_distance_offset[code]+1;
                extra = fixed_distance_extra [code];
                if (extra){
                    offset += stream_read_bit(extra);
                }
                {// выделить фукнцию копирования, второй раз
                    uint8_t * s = dst-offset;
                    int i;
                    for(i=0; i<mlen; i++)
                        *dst++ = *s++;
                }
            }

        }
    }
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
    int d_len = huffan_static_decode(dst, test, sizeof(test)-1) - dst;
    printf("dlen=%d: %s\n", d_len, dst);
    d_len = huffan_static_decode(dst, tst2, sizeof(tst2)) - dst;
    printf("dlen=%d: %s\n", d_len, dst);
    return 0;
}
#endif
