#include <stdint.h>
/*!

 */
static const uint32_t Prime1 = 2654435761U;
static const uint32_t Prime2 = 2246822519U;
static const uint32_t Prime3 = 3266489917U;
static const uint32_t Prime4 =  668265263U;
static const uint32_t Prime5 =  374761393U;

typedef uint32_t uint32x4_t __attribute__((__vector_size__(16)));
typedef struct _HashCtx HashCtx_t;
struct _HashCtx {
    uint32x4_t state;
    uint32x4_t block;
    int offset;
    uint32_t total_length;
};

#define ROTL32(x, n) ((x)<<n | (x)>>(32-n))
void xxh32_init(HashCtx_t *ctx, uint32_t seed)
{
    uint32x4_t state = (uint32x4_t){Prime1 + Prime2, Prime2, 0, -Prime1};
    state+=seed;
    ctx->state = state;
    ctx->block = (uint32x4_t){0};
    ctx->offset=0;
}
void xxh32_update(HashCtx_t *ctx, uint8_t * data, int length)
{
    uint32x4_t state = ctx->state;
    if (ctx->offset){
        int len = (16-ctx->offset)<length?(16-ctx->offset): length;
        __builtin_memcpy((uint8_t*)&ctx->block + ctx->offset, data, len);
        ctx->offset+=len;
        if(ctx->offset!=16) return;
        data  +=len;
        length-=len;

// static U32 XXH32_round(U32 seed, U32 input)
        state = ROTL32(state + ctx->block*Prime2, 13) * Prime1;

        ctx->block = (uint32x4_t){0};
        ctx->offset=0;
    }
    int blocks = length>>4;
    int  i;
    for (i=0; i<blocks; i++) {
        uint32x4_t block;
        __builtin_memcpy(&block, data, 16); data+=16;
        state = ROTL32(state + block*Prime2, 13) * Prime1;
    }
    ctx->state = state;
    length &= 0xF;
    if (length) {
        ctx->block = (uint32x4_t){0};
        __builtin_memcpy(&ctx->block, data, length);
    }
}
uint32_t xxh32_final(HashCtx_t *ctx)
{
    uint32_t hash = (uint32_t)ctx->total_length;
    if (hash >= 16)
        hash += ROTL32(ctx->state[0],  1) +
                ROTL32(ctx->state[1],  7) +
                ROTL32(ctx->state[2], 12) +
                ROTL32(ctx->state[3], 18);
    else
        hash += ctx->state[2] + Prime5;// seed+prime5
    // process remaining bytes in temporary buffer
    if(ctx->offset) {
        uint32_t* s4 = (uint32_t*)&ctx->block;
        int blocks = ctx->offset>>2;
        int i;
        for (i=0; i<blocks; i++)
            hash = ROTL32(hash + s4[i] * Prime3, 17) * Prime4;
        uint8_t* s = (uint8_t*)s4;
        for (i*=4; i<ctx->offset; i++, s++)
            hash = ROTL32(hash + s[i] * Prime5, 11) * Prime1;
    }
/* mix all bits */
    hash ^= hash >> 15;
    hash *= Prime2;
    hash ^= hash >> 13;
    hash *= Prime3;
    hash ^= hash >> 16;
    return hash;
}
/*! \brief Расчет некриптографического хеш, используется в качестве контрольной суммы

 */
uint32_t xxh32(uint32_t hash, uint8_t* data, size_t data_len)
{
    if (data_len>=16){
        uint32x4_t state = (uint32x4_t){Prime1 + Prime2, Prime2, 0, -Prime1};
        state+=hash;
        int blocks = data_len>>4;
        int  i;
        for (i=0; i<blocks; i++) {// вектор 128 бит
            uint32x4_t block;
            __builtin_memcpy(&block, data, 16); data+=16;
            state = ROTL32(state + block*Prime2, 13) * Prime1;
        }
        hash  = ROTL32(state[0],  1) +
                ROTL32(state[1],  7) +
                ROTL32(state[2], 12) +
                ROTL32(state[3], 18);
    } else {
        hash = hash + Prime5;
    }
    hash += data_len;
    data_len &= 0xF;
    int i;
    for (i=0; i < data_len>>2; i++, data+=4)
        hash = ROTL32(hash + *(uint32_t*)data * Prime3, 17) * Prime4;
    for (i*=4; i < data_len; i++, data++)
        hash = ROTL32(hash + *data * Prime5, 11) * Prime1;
    hash ^= hash >> 15;
    hash *= Prime2;
    hash ^= hash >> 13;
    hash *= Prime3;
    hash ^= hash >> 16;
    return hash;
}

/*
Результат компиляции, в цикле получается такой код, цикл за два такта.
        vmovdqa .LC0(%rip), %xmm0
        vpmulld (%rdx), %xmm0, %xmm0
        vpaddd  (%rcx), %xmm0, %xmm0
        vprold  $13, %xmm0, %xmm0
        vpmulld .LC1(%rip), %xmm0, %xmm0

*/
