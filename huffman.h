#ifndef HUFFMAN_H
#define HUFFMAN_H
#include <stdint.h>
typedef struct _deflate deflate_t;
struct _deflate {
    uint32_t stream;
    uint32_t n_bits;
    size_t s_end;
};
uint8_t * deflate (uint8_t *dst, uint8_t *src, int s_len, deflate_t* ctx);
uint32_t crc32_from_block(uint8_t *src, size_t len);
#endif// HUFFMAN_H
