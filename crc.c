#include <stdint.h>
#include <intrin.h>
typedef uint64_t poly64x2_t __attribute__((__vector_size__(16)));
typedef  int64_t v2di __attribute__((__vector_size__(16)));
static inline poly64x2_t LOAD128U(uint8_t* p) {
    return (poly64x2_t)__builtin_ia32_loaddqu((void*)p);
}
#define CL_MUL128(a,b,c) (poly64x2_t)__builtin_ia32_pclmulqdq128 ((v2di)a,(v2di)b,c)
// Структура коэффициентов
struct _CRC_ctx {
	poly64x2_t K34[16];
	poly64x2_t K12;
	poly64x2_t KBP;
};

static const struct _CRC_ctx CRC32B_ctx= {
.KBP = {0xB4E5B025F7011641, 0x1DB710641},
.K12 = {0xAE689191, 0xCCAA009E},
.K34 = {
[ 1] = {0x3F036DC2, 0x40B3A940},// x^{-25}, x^{-89}
[ 2] = {0x7555A0F1, 0x769CF239},// x^{-17}, x^{-81}
[ 3] = {0xCACF972A, 0x5F7314FA},// x^{-9}, x^{-73}
[ 4] = {0xDB710641, 0x5D376816},// x^{-1}, x^{-65}
[ 5] = {0x01000000, 0xF4898239},// x^{7}, x^{-57}
[ 6] = {0x00010000, 0x5FF1018A},// x^{15}, x^{-49}
[ 7] = {0x00000100, 0x0D329B3F},// x^{23}, x^{-41}
[ 8] = {0x00000001, 0xB66B1FA6},// x^{31}, x^{-33}
[ 9] = {0x77073096, 0x3F036DC2},// x^{39}, x^{-25}
[10] = {0x191B3141, 0x7555A0F1},// x^{47}, x^{-17}
[11] = {0x01C26A37, 0xCACF972A},// x^{55}, x^{-9}
[12] = {0xB8BC6765, 0xDB710641},// x^{63}, x^{-1}
[13] = {0x3D6029B0, 0x01000000},// x^{71}, x^{7}
[14] = {0xCB5CD3A5, 0x00010000},// x^{79}, x^{15}
[15] = {0xA6770BB4, 0x00000100},// x^{87}, x^{23}
[ 0] = {0xCCAA009E, 0x00000001},// x^{95}, x^{31}
}};


static uint64_t 	CRC64B_update_N(const struct _CRC_ctx * ctx,  uint64_t crc, uint8_t *data, int len){
	poly64x2_t c = {crc};
	int blocks = (len+15) >> 4;
//__asm volatile("# LLVM-MCA-BEGIN clmul");
	while (--blocks>0){
		c^= (poly64x2_t)LOAD128U(data); data+=16;
		c = CL_MUL128(c, ctx->K12, 0x00) // 128+15 143
		  ^ CL_MUL128(c, ctx->K12, 0x11);//  64+15  79
	}
//__asm volatile("# LLVM-MCA-END");
	len &= 15;
	if (len){
		poly64x2_t v={0};
		__builtin_memcpy(&v, data, len);
		c^= v;
	} else
		c^= (poly64x2_t)LOAD128U(data);
	c = CL_MUL128(c, ctx->K34[len], 0x00) // 15+64
	  ^ CL_MUL128(c, ctx->K34[len], 0x11);// 15
	poly64x2_t t;
	t  = CL_MUL128(c, ctx->KBP, 0x00);
	c ^= CL_MUL128(t, ctx->KBP, 0x10);
	return c[1];
}

#if defined(__PCLMUL__)
uint32_t crc32_from_block(uint8_t *src, size_t len)
{
    uint32_t crc = 0xFFFFFFFF;
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
uint32_t CRC32B_update(uint32_t crc, unsigned char val){
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
