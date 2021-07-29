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

//#define	MATCH_BITS	6
#define	MATCH_MIN	3
//#define	MATCH_MAX	((1 << MATCH_BITS) + (MATCH_MIN - 1))
//#define	OFFSET_MASK	((1 << (16 - MATCH_BITS)) - 1)// =1023
//#define	LEMPEL_SIZE	2048

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

// Hash four bytes starting a p. Используется в
//
// This is Fibonacci hashing, also known as Knuth's multiplicative hash. The
// constant is a prime close to 2^32/phi.


static inline uint32_t lz1_hash_(uint8_t* src, int bits) {
	uint32_t val = (uint32_t) src[0] | ((uint32_t) src[1] << 8);
	return (val * 2654435761ULL) >> (32-16);
	//11400714819323198485ULL
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
// Этот вариант хеша лучше подходит для наших условий
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
static inline uint32_t lzjb_hash(uint8_t* src, int min_len)
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

#define LZ1_MAX_LEN (256+(1<<(16-LZ1_DEPTH))+1)
#define LZ1_LEN_EXT ((1<<(16-LZ1_DEPTH))+2)
#define LZ2_MAX_OFFSET (1<<LZ2_DEPTH)
#define HT_CHAIN (1<<12)
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
	uint16_t bucket [HT_CHAIN+512];
	_Hash_t htable;
	htable.bucket = bucket;
	_hashtable_init(&htable, 512);

	uint8_t copymask=0;
	uint8_t *copy_map;

    int total_depth=0;
    int total_subst=0;
    int max_depth =0;
	while (src < s_end) {
        if((copymask<<=1)==0) {
            copymask = 1;
            copy_map = dst++;//=sizeof(copymask);
            *copy_map = 0;
        }

		uint32_t hash = lz1_hash(src, MATCH_MIN);
		uint32_t y = _hashtable_insert(&htable, hash);

        int max_len = s_end-src;
        // ограничение на длину последовательности
        if(max_len>LZ1_MAX_LEN) max_len=LZ1_MAX_LEN;
		int mlen = 0;

		uint16_t next = y;
		uint16_t moffset;
        int depth =0;
        //int count =64;// ограничение на глубину поиска
		while ((next = _hashtable_next(&htable, next))!=STN_UNDEF /* && count-- */) {// поиск последовательности максимальной длины
            if ((y - next) >= (1<<LZ1_DEPTH)){// ограничение глубины поиска по offset.
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
        //int lz2_depth = 32 - __builtin_clz(src-s_start);
        //if (lz2_depth>LZ2_DEPTH)lz2_depth =LZ2_DEPTH;
        if(mlen>=2 && mlen<(1<<(8-LZ2_DEPTH))+2 && moffset<=LZ2_MAX_OFFSET){// 1 байт на mlen<6
            *copy_map |= copymask;
            if((copymask<<=1)==0) {// кодировние формата упаковки 1байт/2 байта
                copymask = 1;
                copy_map = dst++;
                *copy_map = 0;
            }
            *copy_map |= copymask;
            // формат кодирования mlen(2) | offset(6)
            *dst++ = ((mlen-2)<<LZ2_DEPTH) | (moffset-1);
		} else
		if(mlen>=MATCH_MIN) {
            *copy_map |= copymask;
            if((copymask<<=1)==0) {// кодировние формата упаковки 1байт/2-3 байта
                copymask = 1;
                copy_map = dst++;
                *copy_map = 0;
            }
            if (mlen>=(1<<(16-LZ1_DEPTH))+2) {// формат кодирования 3 байта x1F | x11 | mlen(8) //34..256+33
                *(uint16_t *)dst = ((~0) << LZ1_DEPTH) | (moffset-1);
                dst+=2;
                *dst++ = mlen-LZ1_LEN_EXT;

            } else {// формат кодирования mlen(5)| x11
                *(uint16_t *)dst = ((mlen - MATCH_MIN) << LZ1_DEPTH) | (moffset-1);
                dst+=2;
            }
		} else {
		    *dst++ = *src++;
            continue;
		}
        int i;
        for (i=1; i<mlen; i++) {
            uint32_t hash = lz1_hash(src+i, MATCH_MIN);// может ли заступить за границу
            (void) _hashtable_insert(&htable, hash);
        }
        src+=mlen;

	}
	printf ("avg depth=%1.1f max %d\n", (float)total_depth/total_subst, max_depth);// число подстановок
	return dst;
}
struct _Smap {
    uint32_t moffset:12;// смещение
    uint32_t mlen:11;// длина кода
    uint32_t nlit:9;// длина литерала
};
struct _Smap * lz1_compress_(struct _Smap *map, uint8_t *src, size_t s_len)
{
	uint8_t* s_end  =src+s_len;
	uint8_t* s_start  =src;
	uint16_t bucket [HT_CHAIN+512];
	_Hash_t htable;
	htable.bucket = bucket;
	_hashtable_init(&htable, 512);

    int nlit=0;
	while (src < s_end) {

		uint32_t hash = lz1_hash(src, MATCH_MIN);
		uint32_t y = _hashtable_insert(&htable, hash);

        int max_len = s_end-src;
        // ограничение на длину последовательности
        if(max_len>LZ1_MAX_LEN) max_len=LZ1_MAX_LEN;
		int mlen = 0;

		uint16_t next = y;
		uint16_t moffset;
        //int count =64;// ограничение на глубину поиска
		while ((next = _hashtable_next(&htable, next))!=STN_UNDEF /* && count-- */) {// поиск последовательности максимальной длины
            if ((y - next) >= (1<<LZ1_DEPTH)){// ограничение глубины поиска по offset.
                //chain[next] = STN_UNDEF;
                break;
            }
            int len = lzjb_memcmp(s_start + next, src, max_len);
            if (mlen<len) {
                mlen = len;
                moffset = (y - next);
            }
		}
        if(mlen>=2 && mlen<(1<<(8-LZ2_DEPTH))+2 && moffset<=LZ2_MAX_OFFSET){// 1 байт на mlen<6
            map->nlit = nlit;
		    map->mlen = mlen;
		    map->moffset = moffset-1;
		    map++;
		    nlit= 0;
		} else
		if(mlen>=MATCH_MIN)
        {
            map->nlit = nlit;
		    map->mlen = mlen;
		    map->moffset = moffset-1;
            map++;
		    nlit= 0;
		} else {
		    nlit++;
		    src++;
            continue;
		}
        int i;
        for (i=1; i<mlen; i++) {
            uint32_t hash = lz1_hash(src+i, MATCH_MIN);// может ли заступить за границу
            (void)_hashtable_insert(&htable, hash);
        }
        src+=mlen;

	}
	if(nlit) {
        map->nlit=nlit;
        map->mlen=0;
        map->moffset=0;
        map++;
	}
	return map;
}
uint8_t* lz1_encode(uint8_t *dst, uint8_t *src, struct _Smap *map, int map_len)
{
    uint8_t copymask=0;
    uint8_t *copy_map;
    // nested
    uint8_t _copymask_rotate(uint8_t copymask) {
        if((copymask<<=1)==0) {
            copymask = 1;
            copy_map = dst++;
            *copy_map = 0;
        }
        return copymask;
    }
    struct _Smap *map_end = map+map_len;
    while (map < map_end) {
        if (map->nlit) {
            int i;
            for (i=0;i< map->nlit; i++) {
                copymask = _copymask_rotate(copymask);
                *dst++ = *src++;
            }
            // условие выхода,
            if (map->mlen==0) break;
        }
        copymask = _copymask_rotate(copymask);
        *copy_map |= copymask;
        copymask = _copymask_rotate(copymask);
        uint32_t mlen = map->mlen;
        uint32_t moffset = map->moffset;// единицу вычли уже чтобы помещалось в 11 бит.
        if(mlen>=2 && mlen<(1<<(8-LZ2_DEPTH))+2 && moffset<(1<<LZ2_DEPTH)){// 1 байт на mlen<6
            *copy_map |= copymask;
            // формат кодирования mlen(2) | offset(6)
            *dst++ = ((mlen-2)<<LZ2_DEPTH) | (moffset);
		} else
//		if(mlen>=MATCH_MIN) //-- эти условия выполнены на этапе разбора
        {
            if (mlen>=LZ1_LEN_EXT) {// формат кодирования 3 байта x1F | x11 | mlen(8) //34..256+33
                *(uint16_t *)dst = ((~0) << LZ1_DEPTH) | (moffset);
                dst+=2;
                *dst++ = mlen - LZ1_LEN_EXT;
            } else {// формат кодирования mlen(5)| x11
                *(uint16_t *)dst = ((mlen - MATCH_MIN) << LZ1_DEPTH) | (moffset);
                dst+=2;
            }
		}
		src+=mlen;
		map++;
    }
    return dst;
}
uint8_t* lz1_compress_1(uint8_t *dst, uint8_t *src, size_t s_len)
{
    struct _Smap map[4096];
    int m_count = lz1_compress_(map, src, s_len) - map;
    dst = lz1_encode(dst, src, map, m_count);

    return dst;
}

uint8_t* lz1_decompress(uint8_t *dst, uint8_t *src, size_t s_len)
{
    uint32_t mlen, offset;
    uint32_t nlit=0;
    uint8_t*  s_dst = dst;
    uint8_t*  s_end = src+ s_len;
    uint8_t copymask=0;
    uint8_t copy_map;
    int _copymask_bit_test(){
        if((copymask<<=1)==0) {
            copymask = 1;
            copy_map = *src++;
        }
        return copy_map & copymask;
    }
    while(src<s_end) {
        if (_copymask_bit_test()) {// используется кодирование
            if (_copymask_bit_test()) {// формат кодирования 1 байт + 2бита
                uint32_t  data = *src++;
                offset = (data&((1<<LZ2_DEPTH)-1)) + 1;
                mlen   = (data>>LZ2_DEPTH)+2;// 2-5 байт
            } else {// формат кодирования 2 байта 5-11
                uint32_t data = *(uint16_t* )src; src+=2;
                offset = (data&((1<<LZ1_DEPTH)-1)) + 1;
                mlen   = (data>>LZ1_DEPTH)+3;
                if (mlen==LZ1_LEN_EXT) {
                    mlen += *src++;
                }
            }
            //printf(":%d-%d:'%-.*s'\n", mlen, offset, mlen, dst-offset);
            if(0) printf("%s:%d-%d-%3d:%d\n", (copy_map & copymask)?"F1":"F2", nlit, mlen, offset,
                         32 - __builtin_clz(dst-s_dst));
            nlit=0;
            uint8_t* s = dst-offset;
            int i;
            for(i=0;i<mlen;i++)
                dst[i] = s[i];//*(dst-offset+i);
            dst+= mlen;
        } else {// literal
            nlit++;
            *dst++=*src++;
        }
    }
    if (0&& nlit) printf("F0:%d:\n", nlit);
//    printf ("LZ4 toks =%d: size=%d ratio=%1.2f\n", (int)lz4_nlit, (int)(lz4_tokens+lz4_nlit+lz4_rle), (float)(lz4_tokens+lz4_nlit+lz4_rle)/(dst-s_dst));
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
        size_t zlen = lz1_compress_1(out, buf, len) - out;
        size_t xlen = lz1_decompress(buf2, out, zlen) - buf2;
        cavg += clen;
        zavg += zlen;
        tlen += len;
        if (xlen==len && __builtin_memcmp(buf2, buf, len)==0) printf("..ok ");
        printf("compress ratio =%1.3f / %1.3f\n", (float)clen/len, (float)zlen/len);

	}
	fclose(fp);
	printf("Avg. compress ratio =%1.3f / %1.2f%% (%d B)\n", (float)cavg/tlen, (float)zavg/tlen*100.0, (uint32_t)tlen);
	//lz_decompress(buf, out, clen);
	return 0;
}
