#include <stdint.h>
static const uint64_t Prime1 = 11400714785074694791ULL;
static const uint64_t Prime2 = 14029467366897019727ULL;
static const uint64_t Prime3 =  1609587929392839161ULL;
static const uint64_t Prime4 =  9650029242287828579ULL;
static const uint64_t Prime5 =  2870177450012600261ULL;

typedef uint64_t uint64x4_t __attribute__((__vector_size__(32)));

typedef struct _HashCtx HashCtx_t;
struct _HashCtx {
    uint64_t state[4];
    uint8_t  block[32];
    int offset;
};

#define ROTL64(x, n) ((x)<<n | (x)>>(64-n))
void xxh64_init(HashCtx_t *ctx, uint64_t seed)
{
    uint64x4_t state = (uint64x4_t){Prime1 + Prime2, Prime2, 0, -Prime1};
    state += seed;
    __builtin_memcpy(ctx->state, &state, 16);
    __builtin_bzero (ctx->block, 16);
    ctx->offset = 0;
//    bufferSize  = 0;
//    totalLength = 0;
}
uint64x4_t xxh64_update(uint64x4_t state, uint8_t * data, int length)
{
    uint64x4_t block;
    int  i;
    int blocks = length>>5;
    for (i=0; i<blocks; i++) {
        __builtin_memcpy(&block, data, 32); data+=32;
        state = ROTL64(state + block*Prime2, 31) * Prime1;
    }
    length &= 0x1F;
    if (length) {
    }
    return state;
}
// \sa https://github.com/erthink/t1ha
