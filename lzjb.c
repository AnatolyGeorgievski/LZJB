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
//#include <stdlib.h>// qsort
#include <stdio.h>

extern float huffman_estimate(const uint8_t *code_lengths, const uint16_t *weights, int cl_len);

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

#include <intrin.h>
struct _stream {
    uint8_t copymask;
    uint8_t* stream;
    uint8_t* dst;
};
struct _stream stream_set_bits(struct _stream s, uint8_t bits) {
    //if (__builtin_add_overflow((s.stream<<1), bits&1, &s.stream))
    *s.stream ^= bits!=0? s.copymask: 0;
    //if (_addcarry_u64(bit, s.stream, s.stream, &s.stream))
    if ((s.copymask<<=1)==0)
    {
        s.stream = s.dst; s.dst+=sizeof(*s.stream);
        s.copymask = 1;
        *s.stream=0;
    }
    return s;
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
/*! Алгоритм сжатия в два прохода. На первом составляет таблицы длин. */
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
    uint8_t dummy;
    uint8_t *copy_map= &dummy;
    register uint8_t copymask=0;
    register uint32_t stream=0;
    // nested
    uint8_t _copymask_rotate(uint8_t copymask) {
        if((copymask<<=1)==0) {
            copymask = 1;
            *copy_map = stream;
            copy_map = dst++;
            stream >>=8;
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
        stream |= copymask;
        copymask = _copymask_rotate(copymask);
        uint32_t mlen = map->mlen;
        uint32_t moffset = map->moffset;// единицу вычли уже чтобы помещалось в 11 бит.
        if(mlen>=2 && mlen<(1<<(8-LZ2_DEPTH))+2 && moffset<(1<<LZ2_DEPTH)){// 1 байт на mlen<6
            stream |= copymask;
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
    *copy_map = stream;
    return dst;
}

static uint32_t btree_size(uint8_t *bl_count, int max_bl)
{
    uint32_t bt_size=0;
    int i;
    for(i = 1; i <= max_bl; i++) {
        bt_size += bl_count[i];
    }
    return bt_size;
}

int quick_count=0;
int quick_depth=0;
// Сортировка Шелла (англ. Shell sort) — алгоритм сортировки, являющийся усовершенствованным вариантом сортировки вставками.
static
void shell_sort_0(uint8_t *a, uint16_t * cl_count, int nmemb)//size_t nmemb, size_t size, int (*compar)(const void *, const void *))
{
    int i,j,s;
    for(s=1/*nmemb/2 */; s>0; s/=2){
        for(i=0; i<nmemb; i++){
            uint16_t temp = cl_count[a[i]];
            for(j=i+s; j<nmemb; j+=s){
                //if(cl_count[*pl] > temp || (cl_count[*pl] == temp && pl[0]<pl[1])
                if(cl_count[a[j]] > cl_count[a[i]] ||(cl_count[a[j]] == cl_count[a[i]] && a[j]<a[i])){//array[i] > array[j]){
                    uint8_t t= a[j];
                    a[j] = a[i];
                    a[i] = t;
                }
                quick_depth++;
            }
        }
    }
}
//void __mask_memmove (void* d, void* s, int len) __attribute__((__target__("avx512vl","avx512bw")));
static
void __mask_memcopy (uint8_t* d, uint8_t* s, size_t len) {
#if !defined(__AVX512F__)
    int i;
    for(i=0;i<len;i++)
        d[i] = s[i];//*(dst-offset+i);
#else // 0
#endif
}
static
void __mask_memmove (uint8_t* d, uint8_t* s, size_t len) {
#if !defined(__AVX512F__)
    s+=len, d+=len;
    do{
        *(--d) = *(--s);//*(dst-offset+i);
    } while(--len);
//    __builtin_memmove(d,s,len);
#else // 0
//    if (d==s) return;
    if (0 && d<s) {
        int i, blocks = len>>6;
        for(i=0; i<blocks; i++){
            __m512i v = _mm512_loadu_si512(s);
            _mm512_storeu_si512(d, v);
            s+=64, d+=64;
        }
        if (len&63){
            __mmask64 mask = ~0ULL>>((-len)&63);
            __m512i v = _mm512_maskz_loadu_epi8(mask, s);
            _mm512_mask_storeu_epi8(d, mask, v);
        }
    } else {
/*        int i, blocks = len>>6;
        s+=len, d+=len;
        for(i=0; i<blocks; i++){
            s-=64, d-=64;
            __m512i v = _mm512_loadu_si512(s);
            _mm512_storeu_si512(d, v);
        }*/
        //if (len&63)
        {
            //s-=len&63, d-=len&63;
            __mmask64 mask = ~0ULL>>((-len)&63);
            __m512i v = _mm512_maskz_loadu_epi8(mask, s);
            _mm512_mask_storeu_epi8(d, mask, v);
        }
    }
#endif // 0

}
/*! Сортировка Шелла. Без реккурсии, похоже что быстрая.  */
static void shell_sort_5(uint8_t *a, uint16_t * cl_count, int size)
{
  int inc, i, j, seq[2]={1,4,};// хорошие шаги 2^n*3^m
  int s=1;
  //s = increment(seq, size);
    do {
        inc = seq[s--];
        for (i = inc; i < size; i++) {
            uint8_t m = a[i];
            uint8_t temp = cl_count[m];
            for (j = i-inc; j >= 0; j -= inc){
                // условие надо подбирать или менять направление
                if (cl_count[a[j]] < temp)//cl_count[a[j]] > temp || (cl_count[a[j]] == temp && a[j]<m))
                    break;
                a[j+inc] = a[j];
                quick_count++;
            }
            a[j+inc] = m;
        }
    }while (s >= 0);
}
// https://en.wikipedia.org/wiki/Insertion_sort
// Хорошо работает на упорядоченном списке. Хорошо - значит мало операций копирования.
// Этот алгоритм можно применять на спиках
// Insertion sort для size<=16;
static
void shell_sort(uint8_t *a, uint16_t * cl_count, int size)
{
    uint8_t *pm, *pl;
    const int inc=1;
    for (pm = a + inc; pm < a + size; pm ++){
        uint16_t m = pm[0];
        uint16_t temp = cl_count[pm[0]];
        for (pl = pm-inc; pl >= a; pl-=inc){
            if (cl_count[*pl] > temp || (cl_count[*pl] == temp && pl[0]<m))
                break;
        }
        if(pm>pl+inc){
            quick_count++;
            __mask_memmove(pl+inc+1, pl+inc, pm-(pl+inc));
            pl[inc] = m;
        }
    }
}

static
void quick_sort(uint8_t* array, uint16_t * cl_count, int l, int r)
{
    quick_depth++;
    int m = l;
    int i;
    for (i = l; i <= r; i++) {
        if (cl_count[array[i]] >= cl_count[array[r]]){//mas[i] <= mas[r] )
            uint8_t t= array[m];
            array[m] = array[i];
            array[i] = t;
            quick_count++;
            m++;
        }
    }
    if (l < m-2) quick_sort(array,cl_count, l,m-2);
    if (m < r  ) quick_sort(array,cl_count, m,r);
}
/* Этот вариант для коротких кодов годится. Для маленькх MAX_VALUE

    void counting_sort(vector<int> &mas) {
      vector<int> amount(MAX_VALUE,0);
      for (int i=0;i<mas.size();i++)
        amount[mas[i]]++;
      int pos = -1;
      for (int i=0;i<MAX_VALUE;i++)
        for (int j=0;j<amount[i];j++)
          mas[++pos] = i;
    }
*/


/*! \brief строит гистограмму */
void lz1_hist(struct _Smap *map, int map_len)
{
    uint16_t ml_count[LZ1_DEPTH+1]={0};// распределение длин
    uint16_t dl_count[LZ1_DEPTH+1]={0};// распределение дистанций
    uint16_t cl_count[LZ1_DEPTH+1]={0};// распределение кодов
    int ml_max=0, dl_max=0, cl_max=0;
    int ml_len=0, dl_len=0, cl_len=0;
    int i;
    for(i=0;i<map_len; i++){
        if (map[i].nlit>0) {
            int cl = 32-__builtin_clz(map[i].nlit);
            if (cl_max< cl) cl_max = cl;
            cl_count[cl]++;
            ml_count[0]+=map[i].nlit;
        }
        if (map[i].mlen>0) {
            int ml = 32-__builtin_clz(map[i].mlen-1);
            if (ml_max< ml) ml_max = ml;
            ml_count[ml]++;
            //if (ml)
            {
                int dl = map[i].moffset==0?0:32-__builtin_clz(map[i].moffset);
                if (dl_max< dl) dl_max = dl;
                dl_count[dl]++;
            }
        }
    }
    printf("max: %3d |%3d |%3d\n",ml_max, dl_max, cl_max);
    printf("bits|mlen|dist|nlit\n");
    for(i=0;(i<=ml_max || i<=dl_max || i<=cl_max) && i<LZ1_DEPTH+1; i++){
        printf("%3d |%3d |%3d |%3d\n", i, ml_count[i], dl_count[i], cl_count[i]);
    }
    //alpha_size = nz_count(ml_count, ml_max);// число не нулевых элементов
    for (i=0; i<=ml_max;i++){// вычисляет длину (число символов алфавита)
        if (ml_count[i]>0) ml_len++;
    }
    uint8_t ml_alpha[ml_len+1];
    uint8_t*a = ml_alpha;
    for (i=0; i<=ml_max;i++){// упорядочено по алфавиту
         if(ml_count[i]>0) *a++=i;
    }
    // упорядочить по частоте использования
    shell_sort(ml_alpha, ml_count, ml_len);
    //quick_sort(ml_alpha, ml_count, 0, ml_len-1);
    if (1){
        uint16_t weights[ml_len];
        printf("mlen in order:");
        for(i=0; i<ml_len;i++){
            printf(" %d", ml_alpha[i]);
            weights[i] = ml_count[ml_alpha[i]];
        }
        printf(" (%d,%d)\n", quick_count,quick_depth);
        quick_count=quick_depth=0;
        uint8_t ml_code_lengths[ml_len];//
        extern void huffman_tree(uint8_t *code_lengths, const uint16_t *weights, int cl_len);
        huffman_tree(ml_code_lengths, weights, ml_len);
        float ratio =
        huffman_estimate(ml_code_lengths, weights, ml_len);
        printf("mlen lengths :");
        for(i=0; i<ml_len;i++){
            printf(" %d(%d)", ml_code_lengths[i], weights[i]);
//            weights[i] = ml_count[ml_alpha[i]];
        }
        printf(" ratio=%1.2f%%\n", ratio/3);
    }

    for (i=0; i<=dl_max;i++){
        if (dl_count[i]>0) dl_len++;
    }
    uint8_t dl_alpha[dl_len];
    a = dl_alpha;
    for (i=0; i<=dl_max;i++){// упорядочено по алфавиту
         if(dl_count[i]>0) *a++=i;
    }
    shell_sort(dl_alpha, dl_count, dl_len);
    //quick_sort(dl_alpha, dl_count, 0, dl_len-1);

    if (1){
        printf("dist in order:");
        for(i=0; i<dl_len;i++)
            printf(" %d", dl_alpha[i]);
        printf(" (%d,%d)\n", quick_count,quick_depth);
        quick_count=quick_depth=0;

        printf("freq in order:");
        for(i=0; i<dl_len;i++)
            printf(" %d", dl_count[dl_alpha[i]]);
        printf("\n");
    }

    for (i=0; i<=cl_max;i++){
        if (cl_count[i]>0) cl_len++;
    }
    uint8_t cl_alpha[cl_len];
    a = cl_alpha;
    for (i=0; i<=cl_max;i++){// упорядочено по алфавиту
        if(cl_count[i]>0) *a++=i;
    }
    shell_sort(cl_alpha, cl_count, cl_len);
    //quick_sort(cl_alpha, cl_count, 0, cl_len-1);
    if (1){// длина литерала
        printf("nlit in order:");
        for(i=0; i<cl_len;i++)
            printf(" %d", cl_alpha[i]);
        printf(" (%d,%d)\n", quick_count,quick_depth);
        quick_count=quick_depth=0;
    }


}
extern uint8_t* huffman_fixed_encode(uint8_t *dst, uint8_t *src, struct _Smap *map, int map_len);

uint8_t* lz1_compress_1(uint8_t *dst, uint8_t *src, size_t s_len)
{
    struct _Smap map[4096];
    int m_count = lz1_compress_(map, src, s_len) - map;
    lz1_hist(map, m_count);
    // выбор стратегии на базе гистограммы
    uint8_t buf[4096];
    int h_len = huffman_fixed_encode(buf, src, map, m_count)-buf;
    int d_len = lz1_encode(dst, src, map, m_count)-dst;
    printf("Huffman fixed size=%1.2f%% / LZJB2=%1.2f%%\n",
           (float)h_len*100.f/s_len, (float)d_len*100.f/s_len);

    return dst+d_len;
}
// длины последовательностей обычно маленькие, этот вариант вероятно оптимальный
static inline void _memmove_x64(uint8_t* dst, uint8_t* s, size_t mlen)
{
    int i;
    if(0)
    if (dst-s>=8){// перекрытие
        for(i=0; i<(mlen>>3); i++) {
            *(uint64_t*)dst = *(uint64_t*)s;
            dst+=8, s+=8;
        }
        mlen&=7;
    }
    for(i=0; i<(mlen); i++)
        *dst++ = *s++;
}
// длины последовательностей обычно маленькие, этот вариант вероятно оптимальный
static inline void _memcpy_x64(uint8_t* dst, uint8_t* s, size_t mlen)
{
    int i;
    if(1) {
    for(i=0; i<(mlen>>3); i++) {
        *(uint64_t*)dst = *(uint64_t*)s;
        dst+=8, s+=8;
    }
    mlen&=7;
    }
    for(i=0; i<(mlen); i++)
        *dst++ = *s++;
}
uint8_t* lz1_decompress(uint8_t *dst, uint8_t *src, size_t s_len)
{
    uint32_t mlen, offset;
    uint32_t nlit=0;
    uint8_t*  s_dst = dst;
    uint8_t*  s_end = src+ s_len;
    uint8_t copymask=0;
    int bit_offset=8;
    uint32_t copy_map;
    int stream_read_bits(int n){

        return copy_map & ((1U<<(n))-1);
    }
    int stream_bit_test(){
        if (bit_offset == 8) {
            bit_offset -= 8;
            copy_map = *src++;
        }
        return copy_map & (1U<<(bit_offset++));
        if((copymask<<=1)==0) {
            copymask = 1;
            copy_map = *src++;
        }
        return copy_map & copymask;
    }
    while(src<s_end) {
        if (stream_bit_test()) {// используется кодирование
            if (stream_bit_test()) {// формат кодирования 1 байт + 2бита
                uint32_t  data = *src++;
                offset = (data&((1<<LZ2_DEPTH)-1)) + 1;
                mlen   = (data>>LZ2_DEPTH)+2;// 2-5 байт
            } else {// формат кодирования 2 байта 5-11
                uint32_t data = *(uint16_t* )src; src+=2;
                offset = (data&((1<<LZ1_DEPTH)-1)) + 1;
                mlen   = (data>>LZ1_DEPTH)+3;
                    //printf("&(%d:%d)", offset, mlen);
                if (mlen==LZ1_LEN_EXT) {
                    mlen += *src++;
                }
            }
            if (1 && mlen>0)
            {
                _memmove_x64(dst, dst-offset, mlen);
                dst += mlen;
            } else {
                _memcpy_x64(dst, src, offset);
                src += offset, dst+=offset;
            }
        } else {// literal
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
