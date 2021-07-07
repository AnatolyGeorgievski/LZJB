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

static int lzjb_memcmp(uint8_t *a, uint8_t *b, int max_len) {
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

static inline uint32_t lz1_hash(uint8_t* src, int min_len){
	return *(uint8_t *)src;
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

static void _hashtable_init(_Hash_t *htable,uint32_t nbucket)
{
	htable->nbucket = nbucket;
	htable->nchain = 0;
	//htable->bucket = malloc((nchain+nbucket)*sizeof(uint16_t));
	int i;
	for (i=0; i<nbucket; i++){
		htable->bucket[i]=STN_UNDEF;
	}
}
static uint16_t _hashtable_insert(_Hash_t *htable, uint32_t key)
{
	uint16_t *chain = htable->bucket + htable->nbucket;
	uint16_t y = htable->nchain++;
	uint16_t* head = &htable->bucket[key % htable->nbucket];
	chain[y] = *head;
	*head = y;
	return y;
}
// параметры форматов упаковки смещения
#define LZ1_DEPTH 11
#define LZ2_DEPTH 6
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
	uint16_t offsets[4096];
	uint16_t bucket [4096+256];
	_Hash_t htable;
	htable.bucket = bucket;
	_hashtable_init(&htable, 256);
	uint16_t *chain = htable.bucket+ htable.nbucket;
	uint8_t copymask=0;
	uint8_t *copy_map;
	while (src < s_end) {
        if((copymask<<=1)==0) {
            copymask = 1;
            copy_map = dst++;//=sizeof(copymask);
            *copy_map = 0;
        }
		uint32_t hash = lz1_hash(src, MATCH_MIN);
		uint32_t y = _hashtable_insert(&htable, hash);
		offsets[y] = src - s_start;

        int max_len = s_end-src;
        // ограничение на длину последовательности
        if(max_len>32) max_len=32;
		int mlen = 0;

		uint16_t next = y;
		uint16_t moffset;
		// ограничение на глубину поиска
		while ((next = chain[next])!=STN_UNDEF) {// поиск последовательности максимальной длины
		    if ((offsets[y]-offsets[next])>(1<<LZ1_DEPTH)) {// вращение
                chain[next] = STN_UNDEF;
                break;
		    }
            int len = lzjb_memcmp(s_start + offsets[next], src, max_len);
            if (mlen<len) {
                mlen = len;
                moffset = offsets[y]-offsets[next];
            }
		}
        if(mlen>=2 && mlen<4+2 && moffset<=(1<<LZ2_DEPTH)){// 1 байт на mlen<6
            *copy_map |= copymask;
            if((copymask<<=1)==0) {// кодировние формата упаковки 1байт/2 байта
                copymask = 1;
                copy_map = dst++;
                *copy_map = 0;
            }
            *copy_map |= copymask;
            // формат кодирования mlen(2) | offset(6)
            *dst++ = ((mlen-2)<<LZ2_DEPTH) | (moffset-1);
            int i;
            if (1) // пополняем таблицу
            for (i=1; i<mlen; i++) {
                uint32_t hash = lz1_hash(src+i, MATCH_MIN);// может ли заступить за границу
                uint32_t y = _hashtable_insert(&htable, hash);
                offsets[y] = src - s_start+i;
            }
            src+=mlen;
		} else
		if(mlen>=MATCH_MIN) {
            *copy_map |= copymask;
            if((copymask<<=1)==0) {// кодировние формата упаковки 1байт/2-3 байта
                copymask = 1;
                copy_map = dst++;
                *copy_map = 0;
            }
            // формат кодирования mlen(5)| x11
            if ((mlen - MATCH_MIN)>=32 || (moffset-1)>=2048) _Exit(78);
            *(uint16_t *)dst = ((mlen - MATCH_MIN) << LZ1_DEPTH) | (moffset-1); // = старшие 3 бита, offset
            dst+=2;

            int i;
            if (1)for (i=1; i<mlen; i++) {
                uint32_t hash = lz1_hash(src+i, MATCH_MIN);// может ли заступить за границу
                uint32_t y = _hashtable_insert(&htable, hash);
                offsets[y] = src - s_start+i;
            }
            src+=mlen;
		} else {
            *dst++ = *src++;
		}

	}
	return dst;
}
uint8_t* lz1_decompress(uint8_t *dst, uint8_t *src, size_t s_len)
{
    uint32_t mlen, offset;
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
            if (copy_map & copymask) {// формат кодирования 1 байт
                uint32_t data = *src++;
                offset = (data&((1<<LZ2_DEPTH)-1)) + 1;
                mlen   = (data>>LZ2_DEPTH)+2;// 2-5 байт
            } else {// формат кодирования 2 байта 5-11
                uint32_t data = *(uint16_t* )src; src+=2;
                offset = (data&((1<<LZ1_DEPTH)-1)) + 1;
                mlen   = (data>>LZ1_DEPTH)+3;
            }
            //printf(":%d-%d:'%-.*s'\n", mlen, offset, mlen, dst-offset);
            int i;
            for(i=0;i<mlen;i++)
                dst[i] = *(dst-offset+i);
            dst+= mlen;
        } else {
            //printf("%c",src[0]);
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
	uint8_t buf[4096+512];
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
        printf("compress ratio =%1.3f / %1.3f\n", (float)clen/len, (float)zlen/len);

	}
	printf("Avg. compress ratio =%1.3f / %1.3f\n", (float)cavg/tlen, (float)zavg/tlen);
	//lz_decompress(buf, out, clen);
	return 0;
}
