/*
https://www.kernel.org/doc/Documentation/lzo.txt

      0 1 L D D D S S  (64..127)
           Copy 3-4 bytes from block within 2kB distance
           state = S (copy S literals after this block)
           length = 3 + L
         Always followed by exactly one byte : H H H H H H H H
           distance = (H << 3) + D + 1

      1 L L D D D S S  (128..255)
           Copy 5-8 bytes from block within 2kB distance
           state = S (copy S literals after this block)
           length = 5 + L
         Always followed by exactly one byte : H H H H H H H H
           distance = (H << 3) + D + 1

      L L L D D D S S | H H H H H H H H
           Copy 3-10 bytes from block within 2kB distance
           state = S (copy S literals after this block)
           length = 3 + L
         Always followed by exactly one byte : H H H H H H H H
           distance = (H << 3) + D + 1
 */
#include <stdint.h>
#include <stdio.h>

#define	MATCH_BITS	6
#define	MATCH_MIN	3
#define	MATCH_MAX	((1 << MATCH_BITS) + (MATCH_MIN - 1))
#define	OFFSET_MASK	((1 << (16 - MATCH_BITS)) - 1)// =1023
#define	LEMPEL_SIZE	2048

static int  lzjb_memcmp(uint8_t *a, uint8_t *b, int max_len) {
	int i;
	for (i=0; i<max_len; i++) {
		if(a[i]!=b[i]) break;
	}
	return i;
#if 0
	__mmask64 mask = ~0ULL>>(-max_len & 63);
	mask = _mm512_mask_cmpne_epi8_mask (mask, a, b);
	return __builtin_ctzll(mask);
#endif
}
static void lzjb_memcpy64(uint8_t *dst, uint8_t *src, uint32_t mlen){
	__builtin_memcpy(dst, src, mlen);
#if 0
	__m512i v = _mm512_loadu_epi8(src);
	_mm512_storeu_epi8(dst, v);
#endif
}
static void lzjb_memcpy(uint8_t *dst, uint8_t *src, uint32_t mlen){
	__builtin_memcpy(dst, src, mlen);
#if 0
	__mmask64 mask = ~0ULL>>(-mlen & 63);
	__m512i v = _mm512_maskz_loadu_epi8(mask, cpy);
	_mm512_mask_storeu_epi8(dst, mask, v);
#endif
}

static inline uint32_t lz1_hash_1(uint8_t* src, int min_len){
	uint32_t hash;// = *(uint16_t *)src;// %509;
	hash = (src[0] *33) + (src[1] << 0);
	return hash + hash%31;// + hash%509;
}
static inline uint32_t lz1_hash_2(uint8_t* src, int min_len){
	uint32_t hash;// = *(uint16_t *)src;// %509;
	hash = (src[0] *33) + (src[1] << 0)+5381;
	return hash;
}
static inline uint32_t lz1_hash(uint8_t* src, int min_len){
	uint32_t hash;// = *(uint16_t *)src;// %509;
	hash = (src[0] *66) + src[1];
	return hash;
}
static inline uint32_t lz1_hash_0(uint8_t* src, int min_len){
	uint32_t hash;
	hash = (src[0] << 8) + (src[1] << 0);
	hash += hash >> 9;
	hash += hash >> 5;
	return hash;
}
static uint32_t lzjb_hash(uint8_t* src, int min_len)
{
	uint32_t hash;
	hash = (src[0] << 16) + (src[1] << 8) + src[2];
	hash += hash >> 9;
	hash += hash >> 5;
	return hash;
}
#define STN_UNDEF 0
typedef struct _Hash _Hash_t;
struct _Hash {
    uint16_t *bucket;
    uint16_t nbucket;
    uint16_t nchain;
};
// параметры форматов упаковки смещения
#define LZ1_DEPTH 11
#define LZ2_DEPTH 6
#define HT_CHAIN (1<<11)
#define HT_MASK (HT_CHAIN-1)
static void _hashtable_init(_Hash_t *htable,uint32_t nbucket)
{
	htable->nbucket = nbucket;
	htable->nchain = 0;
	//htable->bucket = malloc((nchain+nbucket)*sizeof(uint16_t));
	int i;
	for (i=0; i<nbucket; i++){
		htable->bucket[i]=STN_UNDEF;
	}
	for (i=0; i<=HT_MASK; i++){
		htable->bucket[nbucket+i]=STN_UNDEF;
	}
}
static uint16_t _hashtable_insert(_Hash_t *htable, uint32_t key)
{
	uint16_t *chain = htable->bucket + htable->nbucket;
	uint16_t y = htable->nchain++;
	uint16_t* head = &htable->bucket[key % htable->nbucket];
    uint16_t prev = *head;
	chain[y & HT_MASK] = prev;//(y - prev)>2048? STN_UNDEF: prev;
	return *head = y;
}
static uint16_t _hashtable_next(_Hash_t *htable, uint16_t ref)
{
    uint16_t *chain = htable->bucket + htable->nbucket;
    return chain[ref & HT_MASK];
}
/*!
Разработка алгоритма потокового сжатия данных
1. Форматы кодирования. Битовые последовательности
Идея заимствованная из LZJB - карусель copymap, добавлена аналогичная идея для выбора формата
2. Ограничения: дистанция 2048 b, длина последовательности до 32 b
 */
uint8_t* lz1_compress(uint8_t *dst, uint8_t *src, size_t s_len)
{
	uint8_t* s_end  =src+s_len;
	uint8_t* s_start  =src;
//	uint16_t offsets[4096];
	uint16_t bucket [4096+512];
	_Hash_t htable;
	htable.bucket = bucket;
	_hashtable_init(&htable, 512);
//	uint16_t *chain = htable.bucket+ htable.nbucket;
	uint8_t copymask=0;
	uint8_t *copy_map;

    int total_depth=0;
    int total_subst=0;
    int max_depth =0;
    int skip=0;
	while (src < s_end) {
        if((copymask<<=1)==0) {
            copymask = 1;
            copy_map = dst++;//=sizeof(copymask);
            *copy_map = 0;
        }

		uint32_t hash = lz1_hash(src, MATCH_MIN);
		uint32_t y = _hashtable_insert(&htable, hash);
//		offsets[y] = src - s_start;
//		if (y!=offsets[y]) printf("$");

        int max_len = s_end-src;
        // ограничение на длину последовательности
        if(max_len>256+33) max_len=256+33;
		int mlen = 0;

		uint16_t next = y;
		uint16_t moffset;
        int depth =0;
        //int count =64;// ограничение на глубину поиска
		while ((next = _hashtable_next(&htable, next))!=STN_UNDEF /* && count-- */) {// поиск последовательности максимальной длины
            if (1&& (y - next) >= (1<<LZ1_DEPTH)){// ограничение глубины поиска по offset.
                //chain[next] = STN_UNDEF;
                break;
            }
            int len = lzjb_memcmp(s_start + next, src, max_len);
            if (mlen<len) {
                mlen = len;
                moffset = (y - next);
            }
            depth++;
		}
#if 1
		total_depth+=depth;
		total_subst++;
		if (depth>max_depth) max_depth= depth;
#endif
        if(mlen>=2 && mlen<4+2 && moffset<=(1<<LZ2_DEPTH)){// 1 байт на mlen<6
            *copy_map |= copymask;
            if((copymask<<=1)==0) {// кодировние формата упаковки 1байт/2 байта
                copymask = 1;
                copy_map = dst++;
                *copy_map = 0;
            }
//            if (mlen<)
            *copy_map |= copymask;
            // формат кодирования mlen(2) | offset(6)
            *dst++ = ((mlen-2)<<LZ2_DEPTH) | (moffset-1);
		} else
		if(mlen>=MATCH_MIN) {
            //printf("F2:%d-%d-%d\n", skip, mlen, moffset);
            *copy_map |= copymask;
            if((copymask<<=1)==0) {// кодировние формата упаковки 1байт/2-3 байта
                copymask = 1;
                copy_map = dst++;
                *copy_map = 0;
            }
            if (mlen>=34) {// формат кодирования x1F | x11 | mlen(8) //34..256+33
                *(uint16_t *)dst = ((0x1F) << LZ1_DEPTH) | (moffset-1);
                dst+=2;
                *dst++ = mlen-34;

            } else {
            // формат кодирования mlen(5)| x11
                *(uint16_t *)dst = ((mlen - MATCH_MIN) << LZ1_DEPTH) | (moffset-1);
                dst+=2;
            }
		} else {
		    skip++;
		    *dst++ = *src++;
            continue;
		}
		skip=0;
        int i;
        for (i=1; i<mlen; i++) {
            uint32_t hash = lz1_hash(src+i, MATCH_MIN);// может ли заступить за границу
            uint32_t y = _hashtable_insert(&htable, hash);
        }
        src+=mlen;

	}
	printf ("avg depth=%1.1f max %d\n", (float)total_depth/total_subst, max_depth);// число подстановок
	return dst;
}
uint8_t* lz1_decompress(uint8_t *dst, uint8_t *src, size_t s_len)
{
    uint32_t mlen, offset;
    uint32_t len;
    uint8_t*  s_end = src+ s_len;
    uint8_t copymask=0;
    uint8_t copy_map;
    while(src<s_end) {
        if((copymask<<=1)==0) {
            copymask = 1;
            copy_map = *src++;
        }
        if (copy_map & copymask) {// используется кодирование
            if((copymask<<=1)==0) {
                copymask = 1;
                copy_map = *src++;
            }
            if (copy_map & copymask) {// формат кодирования 1 байт + 2бита
                uint32_t  data = *src++;
                offset = (data&((1<<LZ2_DEPTH)-1)) + 1;
                mlen   = (data>>LZ2_DEPTH)+2;// 2-5 байт
            } else {// формат кодирования 2 байта 5-11
                uint32_t data = *(uint16_t* )src; src+=2;
                offset = (data&((1<<LZ1_DEPTH)-1)) + 1;
                mlen   = (data>>LZ1_DEPTH)+3;
                if (mlen==34) {
                    //printf("&");
                    mlen += *src++;
                }
            }
            //printf(":%d-%d:'%-.*s'\n", mlen, offset, mlen, dst-offset);
            if(0) printf("%s:%d-%d-%3d:\n", (copy_map & copymask)?"F1":"F2", len, mlen, offset);
            len=0;
            uint8_t* s = (offset==0)? src: dst-offset;
            int i;
            for(i=0;i<mlen;i++)
                dst[i] = s[i];//*(dst-offset+i);
            dst+= mlen;
        } else {// literal
            len++;
            *dst++=*src++;
        }
    }
    return dst;
}


extern size_t
lzjb_compress(void *s_start, void *d_start, size_t s_len, size_t d_len, int n);
#include <locale.h>
int main(int argc, char *argv[]){
    setlocale(LC_ALL, "");
    setlocale(LC_NUMERIC, "C");

	char* filename = NULL;
	if (argc>1) filename = argv[1];
	if (filename==NULL) return -1;
	FILE* fp = fopen(filename, "rb");
	if (fp==NULL) return -1;
	uint8_t buf[4096];
	uint8_t buf2[4096+512];
	uint8_t out[4096+512];
	size_t len, tlen=0;
	uint32_t cavg=0, zavg=0;
	while((len = fread(buf, 1, 4096, fp))!=0) {
        size_t clen = lzjb_compress(buf, out, len, 4096, 0);
        size_t zlen = lz1_compress(out, buf, len) - out;
        size_t xlen = lz1_decompress(buf2, out, zlen) - buf2;
        cavg += clen;
        zavg += zlen;
        tlen += len;
        if (xlen==len && __builtin_memcmp(buf2, buf, len)==0) printf("..ok ");
        printf("compress ratio =%1.3f / %1.3f\n", (float)clen/len, (float)zlen/len);

	}
	fclose(fp);
	printf("Avg. compress ratio =%1.3f / %1.2f%%\n", (float)cavg/tlen, (float)zavg/tlen*100.0);
	//lz_decompress(buf, out, clen);
	return 0;
}
