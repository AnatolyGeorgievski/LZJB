#ifndef HUFFMAN_H
#define HUFFMAN_H
#include <stdint.h>
typedef struct _deflate deflate_t;
struct _deflate {
    uint32_t stream;
    uint32_t n_bits;
    uint8_t * s_end;
};
struct _Smap {
    uint32_t nlit:16;//!< длина литерала
    uint32_t dist:16;//!< смещение
    uint32_t mlen:16;//!< длина кода
};

uint8_t* deflate (uint8_t *dst, uint8_t *src, size_t s_len, deflate_t* ctx);
uint32_t crc32_from_block(uint8_t *src, size_t len);
uint8_t* huffman_fixed_encode(uint8_t *dst, uint8_t *src, struct _Smap *map, int map_len);
float    huffman_estimate(const uint8_t *code_lengths, const uint16_t *weights, int cl_len);
#endif// HUFFMAN_H
