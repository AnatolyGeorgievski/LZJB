/*! \brief Реализация распаковки графики PNG
    [RFC 2083] PNG: Portable Network Graphics, March 1997
    https://datatracker.ietf.org/doc/html/rfc2083
*/
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "huffman.h"


#include <intrin.h>

#define ADLER_PRIME 65521// (1<<16)-(1<<4)+(1)
#define ADLER_MASK 0xFFFFUL;
#if 1
/*! обновление контрольной суммы ADLER32. Начальное значение должно быть 0x1.
    Алгоритм используется в графическом формате PNG.
    */
static uint32_t ADLER32_update_(uint32_t adler, uint8_t *p, size_t len){
	const uint32_t poly = ADLER_PRIME;
	uint32_t s1 = (adler      ) & ADLER_MASK;
	uint32_t s2 = (adler >> 16) & ADLER_MASK;
	if (len) do{
		s1 += *p++;
		if (s1 >= poly) s1 -= poly;
		s2 += s1;
		if (s2 >= poly) s2 -= poly;
	} while (--len);
	return (s2<<16) + s1;
}
#endif // 0
#define mod65521(x) ((x) - 0xFFF1*(unsigned long)(((x)*0x80078071ULL)>>47))
static uint32_t
__attribute__((__target__("sse2","ssse3")))
ADLER32_update_sse2(uint32_t adler, uint8_t *p, size_t len)
{
    const int Nmax = 1024;// Евгенская магия 5552/32
	uint32_t  s1 = (adler      ) & ADLER_MASK;
	uint32_t  s2 = (adler >> 16) & ADLER_MASK;
    int blocks = len>>5;
    if (blocks>0){
        //s2+=s1*(blocks<<5);
        __m128i vs1={0}, vs2={0}, vs3={0}, Z={0};
        __m128i M = _mm_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
        __m128i M1 = _mm_set_epi8(16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31);
        do {
            int n = blocks>Nmax? Nmax: blocks;
            blocks -= n;
            s2+=s1*(n<<5);
            s2-=(s2>>16)*0xFFF1uL;
__asm volatile("# LLVM-MCA-BEGIN adler_sse2");
            do{
                vs3 = _mm_add_epi32(vs3, vs1);
                __m128i v = _mm_lddqu_si128((void*)p); p+=16;
                __m128i v0= _mm_lddqu_si128((void*)p); p+=16;
                __m128i v1 = _mm_sad_epu8(v, Z);
                __m128i v3 = _mm_sad_epu8(v0, Z);
//                __m128i v1 = _mm_maddubs_epi16(v, _mm_set1_epi8(1));
//                __m128i v3 = _mm_maddubs_epi16(v0, _mm_set1_epi8(1));
                __m128i v2 = _mm_maddubs_epi16(v, M1);
                __m128i v4 = _mm_maddubs_epi16(v0, M);
//                v1 = _mm_madd_epi16 (v1, _mm_set1_epi16(1));
//                v3 = _mm_madd_epi16 (v3, _mm_set1_epi16(1));
                v2 = _mm_madd_epi16 (v2, _mm_set1_epi16(1));
                v4 = _mm_madd_epi16 (v4, _mm_set1_epi16(1));
                vs1 = _mm_add_epi32(vs1,v1);
                vs1 = _mm_add_epi32(vs1,v3);
                vs2 = _mm_add_epi32(vs2,v2);
                vs2 = _mm_add_epi32(vs2,v4);
            } while(--n);
__asm volatile("# LLVM-MCA-END adler_sse2");
            // быстрое и неполное редуцирование
            const __m128i P = _mm_set1_epi32(0xFFF1uL);
            vs3 = _mm_sub_epi32 (vs3, _mm_mullo_epi32(_mm_srli_epi32(vs3, 16), P));
            vs2 = _mm_sub_epi32 (vs2, _mm_mullo_epi32(_mm_srli_epi32(vs2, 16), P));
            vs1 = _mm_sub_epi32 (vs1, _mm_mullo_epi32(_mm_srli_epi32(vs1, 16), P));
        } while (blocks>0);

        vs2 = _mm_add_epi32(vs2,vs1);
        vs2 = _mm_add_epi32(vs2,_mm_slli_epi32(vs3,5));
// вставить сюда свертку по последнему фрагменту

        vs1 = _mm_hadd_epi32 (vs1,vs1);
        vs1 = _mm_hadd_epi32 (vs1,vs1);
        vs2 = _mm_hadd_epi32 (vs2,vs2);
        vs2 = _mm_hadd_epi32 (vs2,vs2);
        s1+= _mm_extract_epi32(vs1, 0);
        s2+= _mm_extract_epi32(vs2, 0);
        s1 = mod65521(s1);
        s2 = mod65521(s2);
    }
    int n = len&31;
    if(n){
        s2+=s1*n;
        do {// в таком варианте быстрее, за счет паралельного исполнения
            s1 += p[0];//*p++;// s1 + p0 + p1 + p2 + p3 ... MAX(s1) = 0xFFF0 + 0xFF*n
            s2 += p[0]*n;//s1;// s2 + s1*n + p0*(n) + p1*(n-1) ... MAX(s2) = (n+1)*(0xFFF0 + 0xFF*n/2)
            p++;
        }while (--n);
    }
    s1 = mod65521(s1);
    s2 = mod65521(s2);
    return (s2 << 16) | s1;
}
static uint32_t
__attribute__((__target__("avx2")))
ADLER32_update_avx2(uint32_t adler, uint8_t *p, size_t len)
{
    const int Nmax = 512;// Евгенская магия 5552/32
	uint32_t  s1 = (adler      ) & ADLER_MASK;
	uint32_t  s2 = (adler >> 16) & ADLER_MASK;
    int blocks = len>>6;
    if (blocks>0){
        //s2+=s1*(blocks<<5);
        __m256i vs1={0}, vs2={0}, vs3={0};
//        __m256i E = _mm256_set1_epi8(1);
        __m256i E = _mm256_set1_epi16(1);
        __m256i M = _mm256_set_epi8(
                        0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                        16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31);
        __m256i M1 = _mm256_set_epi8(
            32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,
            48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63);

        do {
            int n = blocks>Nmax? Nmax: blocks;
            blocks -= n;
            s2+=s1*(n<<6);
            s2-=(s2>>16)*0xFFF1uL;
__asm volatile("# LLVM-MCA-BEGIN adler_avx2");
            do{
                vs3 = _mm256_add_epi32(vs3, vs1);
                __m256i v = _mm256_lddqu_si256((void*)p); p+=32;
                __m256i v0 = _mm256_lddqu_si256((void*)p); p+=32;
//                __m256i v1 = _mm256_sad_epu8 (v, Z);
//                __m256i v3 = _mm256_sad_epu8 (v0, Z);
                __m256i v1 = _mm256_maddubs_epi16(v, _mm256_set1_epi8(1));
                __m256i v3 = _mm256_maddubs_epi16(v0, _mm256_set1_epi8(1));
                __m256i v2 = _mm256_maddubs_epi16(v, M1);
                __m256i v4 = _mm256_maddubs_epi16(v0, M);
//                v1 = _mm256_madd_epi16 (v1, _mm256_set1_epi16(1));
                v2 = _mm256_madd_epi16 (v2, E);
                v4 = _mm256_madd_epi16 (v4, E);
                vs1 = _mm256_add_epi32(vs1,v1);
                vs1 = _mm256_add_epi32(vs1,v3);
                vs2 = _mm256_add_epi32(vs2,v2);
                vs2 = _mm256_add_epi32(vs2,v4);
            } while(--n);
__asm volatile("# LLVM-MCA-END adler_avx2");
            // быстрое и неполное редуцирование
            const __m256i P = _mm256_set1_epi32(0xFFF1uL);
            vs3 = _mm256_sub_epi32 (vs3, _mm256_mullo_epi32(_mm256_srli_epi32(vs3, 16), P));
            vs2 = _mm256_sub_epi32 (vs2, _mm256_mullo_epi32(_mm256_srli_epi32(vs2, 16), P));
            vs1 = _mm256_sub_epi32 (vs1, _mm256_mullo_epi32(_mm256_srli_epi32(vs1, 16), P));
        } while (blocks>0);

        vs2 = _mm256_add_epi32(vs2,vs1);
        vs2 = _mm256_add_epi32(vs2,_mm256_slli_epi32(vs3,6));
// вставить сюда свертку по последнему фрагменту
        __m128i v1 = _mm_add_epi32(_mm256_extracti128_si256 (vs1, 0),_mm256_extracti128_si256 (vs1, 1));
        __m128i v2 = _mm_add_epi32(_mm256_extracti128_si256 (vs2, 0),_mm256_extracti128_si256 (vs2, 1));

        v1 = _mm_hadd_epi32 (v1,v1);
        v1 = _mm_hadd_epi32 (v1,v1);
        v2 = _mm_hadd_epi32 (v2,v2);
        v2 = _mm_hadd_epi32 (v2,v2);
        s1+= _mm_extract_epi32(v1, 0);
        s2+= _mm_extract_epi32(v2, 0);
        s1 = mod65521(s1);
        s2 = mod65521(s2);
    }
    int n = len&63;
    if(n){
        s2+=s1*n;
        do {// в таком варианте быстрее, за счет паралельного исполнения
            s1 += p[0];//*p++;// s1 + p0 + p1 + p2 + p3 ... MAX(s1) = 0xFFF0 + 0xFF*n
            s2 += p[0]*n;//s1;// s2 + s1*n + p0*(n) + p1*(n-1) ... MAX(s2) = (n+1)*(0xFFF0 + 0xFF*n/2)
            p++;
        }while (--n);
    }
    s1 = mod65521(s1);
    s2 = mod65521(s2);
    return (s2 << 16) | s1;
}
static uint32_t
__attribute__((__target__("avx512bw")))
ADLER32_update_avx512(uint32_t adler, uint8_t *p, size_t len)
{
    const int Nmax = 512;// Евгенская магия 5552/32
	uint32_t  s1 = (adler      ) & ADLER_MASK;
	uint32_t  s2 = (adler >> 16) & ADLER_MASK;
    int blocks = len>>6;
    if (blocks>0){
        //s2+=s1*(blocks<<5);
        __m512i vs1={0}, vs2={0}, vs3={0}, Z={0};
//        __m256i E = _mm256_set1_epi8(1);
        __m512i E = _mm512_set1_epi16(1);
        __m512i M = _mm512_set_epi8(
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
            16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
            32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,
            48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63);
        //__m512i M1 = _mm512_add_epi8(_mm512_set1_epi8(64),M);

        do {
            int n = blocks>Nmax? Nmax: blocks;
            blocks -= n;
            s2+=s1*(n<<6);
            s2-=(s2>>16)*0xFFF1uL;
__asm volatile("# LLVM-MCA-BEGIN adler_avx512");
            do{
                vs3 = _mm512_add_epi32(vs3, vs1);
                __m512i v = _mm512_loadu_si512((void*)p); p+=64;
                //__m512i v0 = _mm512_loadu_si512((void*)p); p+=64;
                __m512i v1 = _mm512_sad_epu8 (v, Z);
                //__m512i v3 = _mm512_sad_epu8 (v0, Z);
                __m512i v2 = _mm512_maddubs_epi16(v, M);
                //__m512i v4 = _mm512_maddubs_epi16(v0, M);
//                v1 = _mm256_madd_epi16 (v1, _mm256_set1_epi16(1));
                v2 = _mm512_madd_epi16 (v2, E);
//                v4 = _mm512_madd_epi16 (v4, E);
//                vs1 = _mm512_add_epi32(v1,v3);
//                vs2 = _mm512_add_epi32(v2,v4);
                vs1 = _mm512_add_epi32(vs1,v1);
                vs2 = _mm512_add_epi32(vs2,v2);
            } while(--n);
__asm volatile("# LLVM-MCA-END adler_avx512");
            // быстрое и неполное редуцирование
            const __m512i P = _mm512_set1_epi32(0xFFF1uL);
            vs3 = _mm512_sub_epi32 (vs3, _mm512_mullo_epi32(_mm512_srli_epi32(vs3, 16), P));
            vs2 = _mm512_sub_epi32 (vs2, _mm512_mullo_epi32(_mm512_srli_epi32(vs2, 16), P));
            vs1 = _mm512_sub_epi32 (vs1, _mm512_mullo_epi32(_mm512_srli_epi32(vs1, 16), P));
        } while (blocks>0);

        vs2 = _mm512_add_epi32(vs2,vs1);
        vs2 = _mm512_add_epi32(vs2,_mm512_slli_epi32(vs3,6));
// вставить сюда свертку по последнему фрагменту
        __m128i v1 = _mm_add_epi32(
                _mm_add_epi32(_mm512_extracti32x4_epi32(vs1, 0),_mm512_extracti32x4_epi32(vs1, 1)),
                _mm_add_epi32(_mm512_extracti32x4_epi32(vs1, 2),_mm512_extracti32x4_epi32(vs1, 3)));
        __m128i v2 = _mm_add_epi32(
                _mm_add_epi32(_mm512_extracti32x4_epi32(vs2, 0),_mm512_extracti32x4_epi32(vs2, 1)),
                _mm_add_epi32(_mm512_extracti32x4_epi32(vs2, 2),_mm512_extracti32x4_epi32(vs2, 3)));


        v1 = _mm_hadd_epi32 (v1,v1);
        v1 = _mm_hadd_epi32 (v1,v1);
        v2 = _mm_hadd_epi32 (v2,v2);
        v2 = _mm_hadd_epi32 (v2,v2);
        s1+= _mm_extract_epi32(v1, 0);
        s2+= _mm_extract_epi32(v2, 0);
        s1 = mod65521(s1);
        s2 = mod65521(s2);
    }
    int n = len&63;
    if(n){
        s2+=s1*n;
        do {// в таком варианте быстрее, за счет паралельного исполнения
            s1 += p[0];//*p++;// s1 + p0 + p1 + p2 + p3 ... MAX(s1) = 0xFFF0 + 0xFF*n
            s2 += p[0]*n;//s1;// s2 + s1*n + p0*(n) + p1*(n-1) ... MAX(s2) = (n+1)*(0xFFF0 + 0xFF*n/2)
            p++;
        }while (--n);
    }
    s1 = mod65521(s1);
    s2 = mod65521(s2);
    return (s2 << 16) | s1;
}
static uint32_t
__attribute__((__target__("avx512vl","avx512vnni")))
ADLER32_update_vnni(uint32_t adler, uint8_t *p, size_t len)
{
    const int Nmax = 1024;// Евгенская магия 5552/32
	uint32_t  s1 = (adler      ) & ADLER_MASK;
	uint32_t  s2 = (adler >> 16) & ADLER_MASK;
    int blocks = len>>5;
    if (blocks>0){
        //s2+=s1*(blocks<<5);
        __m256i vs1={0}, vs2={0}, vs3={0};
        __m256i E = _mm256_set1_epi8(1);
        __m256i M = _mm256_set_epi8(
                        0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                        16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31);
        do {
            int n = blocks>Nmax? Nmax: blocks;
            blocks -= n;
            s2+=s1*(n<<5);
            s2-=(s2>>16)*0xFFF1uL;
__asm volatile("# LLVM-MCA-BEGIN adler_vnni");
            do{
                vs3 = _mm256_add_epi32(vs3, vs1);
                __m256i v = _mm256_lddqu_si256((void*)p); p+=32;
                vs1 = _mm256_dpbusd_epi32(vs1, v, E);
                vs2 = _mm256_dpbusd_epi32(vs2, v, M);
            } while(--n);
__asm volatile("# LLVM-MCA-END adler_vnni");
            // быстрое и неполное редуцирование
            const __m256i P = _mm256_set1_epi32(0xFFF1uL);
            vs3 = _mm256_sub_epi32 (vs3, _mm256_mullo_epi32(_mm256_srli_epi32(vs3, 16), P));
            vs2 = _mm256_sub_epi32 (vs2, _mm256_mullo_epi32(_mm256_srli_epi32(vs2, 16), P));
            vs1 = _mm256_sub_epi32 (vs1, _mm256_mullo_epi32(_mm256_srli_epi32(vs1, 16), P));
        } while (blocks>0);

        vs2 = _mm256_add_epi32(vs2,vs1);
        vs2 = _mm256_add_epi32(vs2,_mm256_slli_epi32(vs3,5));
// вставить сюда свертку по последнему фрагменту
        __m128i v1 = _mm_add_epi32(_mm256_extracti32x4_epi32(vs1, 0),_mm256_extracti32x4_epi32(vs1, 1));
        __m128i v2 = _mm_add_epi32(_mm256_extracti32x4_epi32(vs2, 0),_mm256_extracti32x4_epi32(vs2, 1));

        v1 = _mm_hadd_epi32 (v1,v1);
        v1 = _mm_hadd_epi32 (v1,v1);
        v2 = _mm_hadd_epi32 (v2,v2);
        v2 = _mm_hadd_epi32 (v2,v2);
        s1+= _mm_extract_epi32(v1, 0);
        s2+= _mm_extract_epi32(v2, 0);
        s1 = mod65521(s1);
        s2 = mod65521(s2);
    }
    int n = len&31;
    if(n){
        s2+=s1*n;
        do {// в таком варианте быстрее, за счет паралельного исполнения
            s1 += p[0];//*p++;// s1 + p0 + p1 + p2 + p3 ... MAX(s1) = 0xFFF0 + 0xFF*n
            s2 += p[0]*n;//s1;// s2 + s1*n + p0*(n) + p1*(n-1) ... MAX(s2) = (n+1)*(0xFFF0 + 0xFF*n/2)
            p++;
        }while (--n);
    }
    s1 = mod65521(s1);
    s2 = mod65521(s2);
    return (s2 << 16) | s1;
}
static uint32_t
__attribute__((__target__("avx512vnni")))
ADLER32_update_512vnni(uint32_t adler, uint8_t *p, size_t len)
{
    const int Nmax = 512;// Евгенская магия 5552/32
	uint32_t  s1 = (adler      ) & ADLER_MASK;
	uint32_t  s2 = (adler >> 16) & ADLER_MASK;
    int blocks = len>>6;
    if (blocks>0){
        //s2+=s1*(blocks<<5);
        __m512i vs1={0}, vs2={0}, vs3={0};
        __m512i E = _mm512_set1_epi8(1);
        __m512i M = _mm512_set_epi8(
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
            16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
            32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,
            48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63);
        do {
            int n = blocks>Nmax? Nmax: blocks;
            blocks -= n;
            s2+=s1*(n<<6);
            s2-=(s2>>16)*0xFFF1uL;
__asm volatile("# LLVM-MCA-BEGIN adler_512vnni");
            do{
                vs3 = _mm512_add_epi32(vs1, vs3);
                __m512i v = _mm512_loadu_epi32((void*)p); p+=64;
                vs1 = _mm512_dpbusd_epi32(vs1, v, E);// задержка 5 тактов, нужно сделать 5 таких команд с разными аргументами
                vs2 = _mm512_dpbusd_epi32(vs2, v, M);
            } while(--n);
__asm volatile("# LLVM-MCA-END adler_512vnni");
            // быстрое и неполное редуцирование
            const __m512i P = _mm512_set1_epi32(0xFFF1uL);
            vs3 = _mm512_sub_epi32 (vs3, _mm512_mullo_epi32(_mm512_srli_epi32(vs3, 16), P));
            vs2 = _mm512_sub_epi32 (vs2, _mm512_mullo_epi32(_mm512_srli_epi32(vs2, 16), P));
            vs1 = _mm512_sub_epi32 (vs1, _mm512_mullo_epi32(_mm512_srli_epi32(vs1, 16), P));
        } while (blocks>0);

        vs2 = _mm512_add_epi32(vs2,vs1);
        vs2 = _mm512_add_epi32(vs2,_mm512_slli_epi32(vs3,6));
// вставить сюда свертку по последнему фрагменту
        __m128i v1 = _mm_add_epi32(
                _mm_add_epi32(_mm512_extracti32x4_epi32(vs1, 0),_mm512_extracti32x4_epi32(vs1, 1)),
                _mm_add_epi32(_mm512_extracti32x4_epi32(vs1, 2),_mm512_extracti32x4_epi32(vs1, 3)));
        __m128i v2 = _mm_add_epi32(
                _mm_add_epi32(_mm512_extracti32x4_epi32(vs2, 0),_mm512_extracti32x4_epi32(vs2, 1)),
                _mm_add_epi32(_mm512_extracti32x4_epi32(vs2, 2),_mm512_extracti32x4_epi32(vs2, 3)));

        v1 = _mm_hadd_epi32 (v1,v1);
        v1 = _mm_hadd_epi32 (v1,v1);
        v2 = _mm_hadd_epi32 (v2,v2);
        v2 = _mm_hadd_epi32 (v2,v2);
        s1+= _mm_extract_epi32(v1, 0);
        s2+= _mm_extract_epi32(v2, 0);
        s1 = mod65521(s1);
        s2 = mod65521(s2);
    }
    int n = len&63;
    if(n){
        s2+=s1*n;
        do {// в таком варианте быстрее, за счет паралельного исполнения
            s1 += p[0];//*p++;// s1 + p0 + p1 + p2 + p3 ... MAX(s1) = 0xFFF0 + 0xFF*n
            s2 += p[0]*n;//s1;// s2 + s1*n + p0*(n) + p1*(n-1) ... MAX(s2) = (n+1)*(0xFFF0 + 0xFF*n/2)
            p++;
        }while (--n);
    }
    s1 = mod65521(s1);
    s2 = mod65521(s2);
    return (s2 << 16) | s1;
}
static
uint32_t ADLER32_update_small(uint32_t adler, uint8_t *p, size_t len)
{
    const size_t Nmax = 5552;// Евгенская магия
	uint32_t  s1 = (adler      ) & ADLER_MASK;
	uint32_t  s2 = (adler >> 16) & ADLER_MASK;
    do {
        int n = len>Nmax? Nmax: len;
        len-=n;
        do {
            s1 += *p++;// s1 + p0 + p1 + p2 + p3 ... MAX(s1) = 0xFFF0 + 0xFF*n
            s2 += s1;// s2 + s1*n + p0*(n) + p1*(n-1) ... MAX(s2) = (n+1)*(0xFFF0 + 0xFF*n/2)
        }while (--n);
        s1-= (s1>>16)*0xFFF1u;// не полное редуцирование
        s2-= (s2>>16)*0xFFF1u;
    } while (len);
    s1 = mod65521(s1);
    s2 = mod65521(s2);
    return (s2 << 16) | s1;
}
static
uint32_t ADLER32_update_default(uint32_t adler, uint8_t *p, size_t len)
{
	uint32_t  s1 = (adler      ) & ADLER_MASK;
	uint32_t  s2 = (adler >> 16) & ADLER_MASK;
    inline uint8_t* adler_sum_4(uint8_t *p, uint32_t n){
        s2+=s1*n;
        uint32_t a1=0, a2=0, a3=0;
        do {
           a3 += a1;
           a1 += p[0]   + p[1]   + p[2]   + p[3];
           a2 += p[0]*3 + p[1]*2 + p[2]*1 + p[3]*0;
           p+=4;
        } while ((n-=4)>0);
        s2+=a3*4 + a2 + a1;
        s1+=a1;
        return p;
    }
    inline uint8_t* adler_sum_2(uint8_t *p, uint32_t n){
        s2+=s1*n;
        uint32_t a1=0, a2=0, a3=0;
        do {
           a3 += a1;
           a1 += p[0]   + p[1];
           a2 += p[0]*1 + p[1]*0;
           p+=2;
        } while ((n-=2)>0);
        s2+=a3*2 + a2 + a1;
        s1+=a1;
        return p;
    }
    const int N=12;
    const int L=2;
    int blocks = len>>N;
    if(blocks) {
        uint32_t a1=0, a2=0, a3=0;
        do {
            s2+=s1<<N;
            s2-= (s2>>16)*0xFFF1u;
            int n = 1<<N;
            do {
               a3 += a1;
               int i;
               for (i=0;i<L; i++) {
                    a1 += p[i];
                    a2 += p[i]*(L-i-1);
               }
               p+=L;
            } while ((n-=L)>0);
            a1-= (a1>>16)*0xFFF1u;
            a2-= (a2>>16)*0xFFF1u;
            a3-= (a3>>16)*0xFFF1u;
        } while(--blocks);
        s2+=a3*L + a2 + a1;
        s1+=a1;
        s1-= (s1>>16)*0xFFF1u;
        s2-= (s2>>16)*0xFFF1u;
    }
    int n = len &((1<<N)-1);
    if(n) {
        s2+=s1*n;
        do {// в таком варианте быстрее, за счет паралельного исполнения
            s1 += p[0];//*p++;// s1 + p0 + p1 + p2 + p3 ... MAX(s1) = 0xFFF0 + 0xFF*n
            s2 += p[0]*n;//s1;// s2 + s1*n + p0*(n) + p1*(n-1) ... MAX(s2) = (n+1)*(0xFFF0 + 0xFF*n/2)
            p++;
        } while (--n);
    }
    s1 = mod65521(s1);
    s2 = mod65521(s2);
    return (s2 << 16) | s1;
}
#if defined(__AVX512VL__) && defined(__AVX512VNNI__)
 #define ADLER32_update ADLER32_update_vnni
#elif defined(__AVX512BW__)
 #define ADLER32_update ADLER32_update_avx512
#elif defined(__AVX2__)
 #define ADLER32_update ADLER32_update_avx2
#elif defined(__SSE2__) && defined(__SSSE3__)
 #define ADLER32_update ADLER32_update_sse2
#else
 #define ADLER32_update ADLER32_update_default
#endif // defined

#define BE32(x) __builtin_bswap32(x)
static const uint8_t magic[] = {137, 80, 78, 71, 13, 10, 26, 10};
struct _PNG_Hdr {
    uint32_t width;     //Width:              4 bytes
    uint32_t height;    //Height:             4 bytes
    uint8_t bit_depth;  //Bit depth:          1 byte
    uint8_t color_type; //Color type:         1 byte
    uint8_t compression;//Compression method: 1 byte
    uint8_t filter;     //Filter method:      1 byte
    uint8_t interlace;  //Interlace method:   1 byte
};
struct _Zlib_Hdr {
    uint8_t compression;// Compression method/flags code: 1 byte
    uint8_t flags;  // Additional flags/check bits:   1 byte
//      Compressed data blocks:        n bytes
//      Check value:                   4 bytes
};

int png_to_image(uint8_t *src, size_t s_len)
{
    if (__builtin_memcmp(src, magic, 8)!=0) return -1;
    struct _PNG_Hdr *hdr = NULL;
    int chunk=0;
    uint8_t *dst=NULL;
    uint8_t *s_end = src+s_len;
    src+=8;
    if(1) printf("PNG:\n");
    while(src<s_end) {
        uint32_t length = BE32(*(uint32_t *)src); src+=4;
        char *ctype  = (char*)src; src+=4;
        uint8_t *cdata  = src; src+=length;
        if (src>s_end) break;
        uint32_t crc    = BE32(*(uint32_t *)(src)); src+=4;
        if (crc != crc32_from_block(cdata-4, length+4)) return -1;
        if (1) printf ("%-.*s crc=%08X  len=%d\n", 4, ctype, crc, length);
        if (strncasecmp(ctype,"ihdr", 4)==0) {
            hdr = (void*)cdata;
            printf ("\twidth=%d height=%d bits=%d colors=%d compression=%d\n", BE32(hdr->width), BE32(hdr->height),
                    hdr->bit_depth, hdr->color_type, hdr->compression);
            dst = malloc(BE32(hdr->width)*BE32(hdr->height)*4+2048);// кто то пишет мимо
        } else
        if (strncasecmp(ctype,"text", 4)==0) {
            int len = strlen((char*)cdata)+1;
            printf("\t%s: %-.*s\n", (char*)cdata, length-len, (char*)cdata+len);
        } else
        if (strncasecmp(ctype,"idat", 4)==0) {
            if (chunk==0 && cdata[0]==0){
                uint32_t len = *(uint16_t*)(cdata+1);
                printf("\tchunk length=%d\n", len);

            } else {
                deflate_t ctx={0};
                size_t d_len = deflate(dst, cdata+2, length-6, &ctx)-dst;
                printf("\ncompression =%1.2f%%\n", (float)(length-6)*100.f/(d_len));
                //size_t ilen = BE32(hdr->width)*BE32(hdr->height)*4;
                //if (ilen<d_len)  d_len = ilen;
                uint32_t crc = crc32_from_block(dst, d_len);//BE32(hdr->width)*BE32(hdr->height)*4);
                uint32_t adler = ADLER32_update(1, dst, d_len);
                if(adler==BE32(*(uint32_t *)(cdata+length-4)))
                    printf("ADLER32 Check sum ..ok  %08X, len=%d crc=%08X\n", adler, (int)d_len, crc);
                else
                    printf("ADLER32 Check sum ..fail %08X, len=%d crc=%08X\n", adler, (int)d_len, crc);
            }
        } else
        if (strncasecmp(ctype,"iend", 4)==0) {
            break;
        }
    }
    return 0;
}

static int _get_contents(char* filename, char** contents, size_t *length, void* error)
{
    struct stat     statbuf;
    int res = stat(filename, &statbuf);
    if (res==0) {
        char* data = malloc(statbuf.st_size);
        FILE * f = fopen(filename, "rb");
        if (f!=NULL) {
            *length = fread(data,1,statbuf.st_size, f);
            *contents = data;
            fclose(f);
        }
    }
    return res==0;
}
#if defined(TEST_PNG)
int main()
{
    char* filename = "test6.png";
    uint8_t *contents=NULL;
    size_t length=0;
    _get_contents(filename, (char**)&contents, &length, NULL);
    png_to_image(contents, length);
    free(contents);
    return 0;
}
#endif // defined

