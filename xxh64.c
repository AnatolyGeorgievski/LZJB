#include <stdint.h>
#include <stddef.h>
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
static inline uint64_t ROUND64(uint64_t x){
    return ROTL64(x*Prime2,31)*Prime1;
}
static inline uint64_t MERGE64(uint64_t hash, uint64_t x){
    return (hash ^ ROUND64(x))*Prime1 + Prime4;
}
uint64_t xxh64(uint64_t hash, uint8_t* data, size_t data_len)
{
    if (data_len>=32){
        uint64x4_t state = (uint64x4_t){Prime1 + Prime2, Prime2, 0, -Prime1};
        state+=hash;
        int blocks = data_len>>5;
        int  i;
        for (i=0; i<blocks; i++) {// вектор 128 бит
            uint64x4_t block;
            __builtin_memcpy(&block, data, 32); data+=32;
            state = ROTL64(state + block*Prime2, 31) * Prime1;
        }
        hash  = ROTL64(state[0],  1) +
                ROTL64(state[1],  7) +
                ROTL64(state[2], 12) +
                ROTL64(state[3], 18);
        hash  = MERGE64(hash, state[0]);
        hash  = MERGE64(hash, state[1]);
        hash  = MERGE64(hash, state[2]);
        hash  = MERGE64(hash, state[3]);
    } else {
        hash = hash + Prime5;
    }
    hash += data_len;
    data_len &= 0x1F;
    int i;
    for (i=0; i < data_len>>3; i++, data+=8)
        hash = ROTL64(hash ^ ROUND64(*(uint64_t*)data), 27) * Prime1 + Prime4;
    for (i*=2; i < data_len>>2; i++, data+=4)
        hash = ROTL64(hash ^ *(uint32_t*)data * Prime1, 23) * Prime2 + Prime3;
    for (i*=4; i < data_len; i++, data++)
        hash = ROTL64(hash ^ *data * Prime5, 11) * Prime1;
    hash ^= hash >> 33;
    hash *= Prime2;
    hash ^= hash >> 29;
    hash *= Prime3;
    hash ^= hash >> 32;
    return hash;
}

// \sa https://github.com/erthink/t1ha
