/*! \file crc.c
    Copyright (C) 2021 Anatoly Georgievskii <anatoly.georgievski@gmail.com>

    \see Intel@ Fast CRC Computation Using PCLMULQDQ Instruction
Algorithm 1. Barrett Reduction Algorithm for a degree-32 polynomial modulus
(polynomials defined over GF(2))
Input: degree-63 polynomial R(x), degree-32 polynomial P(x), μ = [(x^64 / P(x))]
Output: C(x) = R(x) mod P(x)
Step 1: T1(x) = [ (R(x)/x^32) ] • μ
Step 2: T2(x) = [ (T1(x)/x^32) ] • P(x)
Step 3: C(x) = R(x) &oplus; T2(x) mod x^32
After step 3, the 32 high-order coefficients of C will be 0.

 */

#include <stdint.h>
#if defined(__PCLMUL__)

#include <intrin.h>
typedef uint64_t poly64x2_t __attribute__((__vector_size__(16)));
typedef uint64_t poly64x4_t __attribute__((__vector_size__(32)));
typedef  int64_t v2di __attribute__((__vector_size__(16)));
static inline poly64x2_t LOAD128U(uint8_t* p) {
    return (poly64x2_t)_mm_lddqu_si128((void*)p);
}
static inline poly64x2_t LOAD128U_maskz(uint8_t* p, int len) {
    poly64x2_t v={0};
	__builtin_memcpy(&v, p, len);
    return (poly64x2_t)v;
}

#define CL_MUL128(a,b,c) (poly64x2_t) _mm_clmulepi64_si128((__m128i)a,(__m128i)b,c)
// Структура коэффициентов
struct _CRC_ctx {
	poly64x2_t K34[16];
	poly64x2_t KBP;//!< final reduction: Barrett's constant and Polynom
	poly64x2_t K12;//!< fold by 1 (128 bits)
	poly64x2_t KF2;//!< fold by 2
	poly64x2_t KF3;//!< fold by 3
	poly64x2_t KF4;//!< fold by 4
	poly64x2_t KF5;//!< fold by 5
};
static const struct _CRC_ctx CRC32B_ctx= {
.KBP = {0xB4E5B025F7011641, 0x1DB710641},// Barrett reduction constant and polynom 33bit
.KF5 = {0x1C279815, 0xAE0B5394},
.KF4 = {0x8F352D95, 0x1D9513D7},// x^{543}, x^{479}
.KF3 = {0x3DB1ECDC, 0xAF449247},// x^{415}, x^{351}
.KF2 = {0xF1DA05AA, 0x81256527},// x^{287}, x^{223}
.K12 = {0xAE689191, 0xCCAA009E},// x^{159}, x^{ 95}
.K34 = {
[ 1] = {0x3F036DC2, 0x40B3A940},// x^{-25}, x^{-89}
[ 2] = {0x7555A0F1, 0x769CF239},// x^{-17}, x^{-81}
[ 3] = {0xCACF972A, 0x5F7314FA},// x^{-9}, x^{-73}
[ 4] = {0xDB710641, 0x5D376816},// x^{-1}, x^{-65}
[ 5] = {0x01000000, 0xF4898239},// x^{ 7}, x^{-57}
[ 6] = {0x00010000, 0x5FF1018A},// x^{15}, x^{-49}
[ 7] = {0x00000100, 0x0D329B3F},// x^{23}, x^{-41}
[ 8] = {0x00000001, 0xB66B1FA6},// x^{31}, x^{-33}
[ 9] = {0x77073096, 0x3F036DC2},// x^{39}, x^{-25}
[10] = {0x191B3141, 0x7555A0F1},// x^{47}, x^{-17}
[11] = {0x01C26A37, 0xCACF972A},// x^{55}, x^{-9}
[12] = {0xB8BC6765, 0xDB710641},// x^{63}, x^{-1}
[13] = {0x3D6029B0, 0x01000000},// x^{71}, x^{ 7}
[14] = {0xCB5CD3A5, 0x00010000},// x^{79}, x^{15}
[15] = {0xA6770BB4, 0x00000100},// x^{87}, x^{23}
[ 0] = {0xCCAA009E, 0x00000001},// x^{95}, x^{31}
}};


uint32_t
__attribute__((__target__("vpclmulqdq","avx512vl")))
	CRC64B_update_N_(const struct _CRC_ctx * ctx,  uint32_t crc, uint8_t *data, int len){
	poly64x4_t c = {crc};
	int blocks = (len+15) >> 4;
    if (blocks>3) {// fold by 256 bits
        __m256i KF2= _mm256_set_m128i ((__m128i)ctx->KF2,(__m128i)ctx->KF2);
        do {
            c^= (poly64x4_t)_mm256_lddqu_si256((void*)(data));
            c = (poly64x4_t)_mm256_clmulepi64_epi128((__m256i)c, KF2, 0x00) ^ _mm256_clmulepi64_epi128((__m256i)c, KF2, 0x11);
            blocks-=1, data+=32;
        } while(blocks>3);
        poly64x2_t c1 = (poly64x2_t)_mm256_extracti128_si256 ((__m256i)c, 1);
        poly64x2_t c0 = (poly64x2_t)_mm256_castsi256_si128((__m256i)c);
        c0^= (poly64x2_t)_mm_lddqu_si128((void*)(data));
        c0  = c1 ^ CL_MUL128(c0, ctx->K12, 0x00) ^ CL_MUL128(c0, ctx->K12, 0x11);
    }
	len &= 15;
	return c[1];
}

static uint32_t
__attribute__((__target__("pclmul")))
CRC64B_update_N(const struct _CRC_ctx * ctx,  uint32_t crc, uint8_t *data, int len){
	poly64x2_t c = {crc};
	int blocks = (len+15) >> 4;
    if (1 && blocks>9) {// fold by 5x128 bits
        poly64x2_t c1 = {0}, c2 = {0}, c3 = {0}, c4 = {0};
        const poly64x2_t KF5 = ctx->KF5;
__asm volatile("# LLVM-MCA-BEGIN clmul");
        do {
            c ^= (poly64x2_t)_mm_lddqu_si128((void*)(data   ));
            c1^= (poly64x2_t)_mm_lddqu_si128((void*)(data+16));
            c2^= (poly64x2_t)_mm_lddqu_si128((void*)(data+32));
            c3^= (poly64x2_t)_mm_lddqu_si128((void*)(data+48));
            c4^= (poly64x2_t)_mm_lddqu_si128((void*)(data+64));
            c  = CL_MUL128(c , KF5, 0x00) ^ CL_MUL128(c , KF5, 0x11);
            c1 = CL_MUL128(c1, KF5, 0x00) ^ CL_MUL128(c1, KF5, 0x11);
            c2 = CL_MUL128(c2, KF5, 0x00) ^ CL_MUL128(c2, KF5, 0x11);
            c3 = CL_MUL128(c3, KF5, 0x00) ^ CL_MUL128(c3, KF5, 0x11);
            c4 = CL_MUL128(c4, KF5, 0x00) ^ CL_MUL128(c4, KF5, 0x11);
            blocks-=5, data+=80;
        } while(blocks>9);
__asm volatile("# LLVM-MCA-END clmul");
        c ^= (poly64x2_t)_mm_lddqu_si128((void*)(data   ));
        c1^= (poly64x2_t)_mm_lddqu_si128((void*)(data+16));
        c2^= (poly64x2_t)_mm_lddqu_si128((void*)(data+32));
        c3^= (poly64x2_t)_mm_lddqu_si128((void*)(data+48));
        c  = c4
           ^ CL_MUL128(c , ctx->KF4, 0x00) ^ CL_MUL128(c , ctx->KF4, 0x11)
           ^ CL_MUL128(c1, ctx->KF3, 0x00) ^ CL_MUL128(c1, ctx->KF3, 0x11)
           ^ CL_MUL128(c2, ctx->KF2, 0x00) ^ CL_MUL128(c2, ctx->KF2, 0x11)
           ^ CL_MUL128(c3, ctx->K12, 0x00) ^ CL_MUL128(c3, ctx->K12, 0x11);
        blocks-=4, data+=64;
    }
    if (0 && blocks>7) {// fold by 4x128 bits
        poly64x2_t c1 = {0}, c2 = {0}, c3 = {0};
__asm volatile("# LLVM-MCA-BEGIN clmul");
        do {
            c ^= (poly64x2_t)_mm_lddqu_si128((void*)(data   ));
            c1^= (poly64x2_t)_mm_lddqu_si128((void*)(data+16));
            c2^= (poly64x2_t)_mm_lddqu_si128((void*)(data+32));
            c3^= (poly64x2_t)_mm_lddqu_si128((void*)(data+48));
            c  = CL_MUL128(c , ctx->KF4, 0x00) ^ CL_MUL128(c , ctx->KF4, 0x11);
            c1 = CL_MUL128(c1, ctx->KF4, 0x00) ^ CL_MUL128(c1, ctx->KF4, 0x11);
            c2 = CL_MUL128(c2, ctx->KF4, 0x00) ^ CL_MUL128(c2, ctx->KF4, 0x11);
            c3 = CL_MUL128(c3, ctx->KF4, 0x00) ^ CL_MUL128(c3, ctx->KF4, 0x11);
            blocks-=4, data+=64;
        } while(blocks>7);
__asm volatile("# LLVM-MCA-END clmul");
        c ^= (poly64x2_t)_mm_lddqu_si128((void*)(data   ));
        c1^= (poly64x2_t)_mm_lddqu_si128((void*)(data+16));
        c2^= (poly64x2_t)_mm_lddqu_si128((void*)(data+32));
        c  = c3
           ^ CL_MUL128(c , ctx->KF3, 0x00) ^ CL_MUL128(c , ctx->KF3, 0x11)
           ^ CL_MUL128(c1, ctx->KF2, 0x00) ^ CL_MUL128(c1, ctx->KF2, 0x11)
           ^ CL_MUL128(c2, ctx->K12, 0x00) ^ CL_MUL128(c2, ctx->K12, 0x11);
        blocks-=3, data+=48;
    }
    if (0 && blocks>3) {// fold by 2x128 bits
        poly64x2_t c1 = {0};
        do {
//            _mm_prefetch(data,  _MM_HINT_NTA);
            c ^= (poly64x2_t)_mm_lddqu_si128((void*)(data   ));
            c1^= (poly64x2_t)_mm_lddqu_si128((void*)(data+16));
            //c  = CL_MUL128(c, ctx->K12, 0x00) ^ CL_MUL128(c, ctx->K12, 0x11);
            c  = CL_MUL128(c, ctx->KF2, 0x00) ^ CL_MUL128(c, ctx->KF2, 0x11);
            c1 = CL_MUL128(c1, ctx->KF2, 0x00) ^ CL_MUL128(c1, ctx->KF2, 0x11);
            blocks-=2, data+=32;
        } while(blocks>3);
        c ^= (poly64x2_t)_mm_lddqu_si128((void*)(data));
        c  = c1 ^ CL_MUL128(c, ctx->K12, 0x00) ^ CL_MUL128(c, ctx->K12, 0x11);
        blocks-=1, data+=16;
    }
    if (blocks>1) {// fold by 128 bits
        do {
            c^= (poly64x2_t)_mm_lddqu_si128((void*)(data));
            c = CL_MUL128(c, ctx->K12, 0x00) ^ CL_MUL128(c, ctx->K12, 0x11);
            blocks-=1, data+=16;
        } while(blocks>1);
    }
	len &= 15;
	if (len){
		c^= (poly64x2_t)LOAD128U_maskz(data, len);
	} else
		c^= (poly64x2_t)LOAD128U(data);
	c = CL_MUL128(c, ctx->K34[len], 0x00) // 15+64
	  ^ CL_MUL128(c, ctx->K34[len], 0x11);// 15
	poly64x2_t t;
	t  = CL_MUL128(c, ctx->KBP, 0x00);
	c ^= CL_MUL128(t, ctx->KBP, 0x10);
	return c[1];
}

uint32_t crc32_from_block(uint8_t *src, size_t len)
{
    uint32_t crc = ~0UL;
    crc = CRC64B_update_N(&CRC32B_ctx, crc, src, len);
    //crc = CRC32B_update_N(crc, src, len);
    return ~crc;
}
#else
static const uint32_t CRC32B_Lookup4[16]={
0x00000000, 0x1DB71064, 0x3B6E20C8, 0x26D930AC,
0x76DC4190, 0x6B6B51F4, 0x4DB26158, 0x5005713C,
0xEDB88320, 0xF00F9344, 0xD6D6A3E8, 0xCB61B38C,
0x9B64C2B0, 0x86D3D2D4, 0xA00AE278, 0xBDBDF21C
};
/*! \brief маленькая функция с маленькой табличкой */
static uint32_t CRC32B_update(uint32_t crc, unsigned char val){
	crc^= val;
	crc = (crc>>4) ^ CRC32B_Lookup4[crc & 0xF];
	crc = (crc>>4) ^ CRC32B_Lookup4[crc & 0xF];
	return crc;
}
uint32_t crc32_from_block(uint8_t *src, size_t len)
{
    int i;
    uint32_t crc = ~0UL;
    for(i=0; i<len; i++)
        crc = CRC32B_update(crc, src[i]);
    return ~crc;
}
#endif
