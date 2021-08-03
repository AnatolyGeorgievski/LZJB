// https://www.ietf.org/rfc/rfc1952.txt
/*! Распаковывает GZip формат
    [RFC 1952]
    https://www.ietf.org/rfc/rfc1952.txt

    Понадобилось для отладки Deflate на больших файлах
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "huffman.h"

#define FTEXT   (1<<0)
#define FHCRC   (1<<1)
#define FEXTRA  (1<<2)
#define FNAME   (1<<3)
#define FCOMMENT (1<<4)

#if __BYTE_ORDER__==__ORDER_LITTLE_ENDIAN__
#define LE32(x) (x)
#define LE16(x) (x)
#define BE32(x) __builtin_bswap32(x)
#define BE16(x) __builtin_bswap16(x)
#else
#define BE32(x) (x)
#define BE16(x) (x)
#endif
#define ASSERT(x) if (!(x)) { printf("ERR: ASSERT "#x "\n"); return (-1);}
int gunzip(uint8_t*src, size_t s_len)
{
    uint8_t * s_end = src+s_len;
    uint8_t buff[4096];
    // header
    const uint16_t magic=0x8B1F;
    uint8_t compression_method;// 8=Deflate
    uint8_t flg;//флаги
    //FILE *fp =fopen(filename, "rb");
//    fread((char*)buff, 1, 4096, fp);
//    uint8_t *src = buff;
    ASSERT (*(uint16_t*)src == magic);

    src+=2;
    compression_method = *src++;
    flg = *src++;
    uint32_t mtime=*(uint32_t*) src; src+=4;
    // xfg OS
    src+=2;
    if (flg&FEXTRA) {
        uint16_t len = (*(uint16_t*)src); src+=2;
        src+=len;
    }
    if (flg&FNAME) {
        char* name = (char*)src;
        src+= strlen(name)+1;
        printf("filename: %s\n",name);
    }
    if (flg&FCOMMENT) {
        char* comment = (char*)src;
        src+= strlen(comment)+1;
        printf("comment: %s\n",comment);
    }
    if (flg&FHCRC) {
        uint16_t crc16 = (*(uint16_t*)src); src+=2;
        //ASSERT (*(uint16_t*)src == crc16);
    }
    // CRC32 ISIZE
    uint32_t crc32 = LE32(*(uint32_t*)(s_end-8));// src+=4;
    uint32_t isize = LE32(*(uint32_t*)(s_end-4));// src+=4;
    printf("initial size: %d, compression ratio: %1.2f%%\n", isize, (float)(s_end-src-8)*100.f/isize);
    uint8_t *dst  = malloc(isize);

    // compressed blocks
    deflate_t ctx={0};
    size_t d_size = 0;

    do {
        size_t chunk_size = deflate(dst+d_size, src, s_end-src-8, &ctx)-(dst+d_size);
        d_size+=chunk_size;
        src = ctx.s_end;
        //printf("chunk size: %d %1.1f\n", chunk_size, (float)(d_size)*100.f/isize);
    } while (d_size < isize);

    ASSERT((uint32_t)d_size == isize);
    uint32_t crc = crc32_from_block(dst, d_size);
    if (crc == crc32)
        printf("CRC32: %08X ..ok\n", crc32);
    else
        printf("CRC32: %08X<>%08X ..fail\n", crc32, crc);

    ASSERT(crc == crc32);
    return 0;
}
#if defined(TEST_GZIP)
#include <sys/stat.h>
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
int main()
{
    char* filename = "test2.gz";
    uint8_t *contents=NULL;
    size_t length=0;
    // заменить на memory map
    if(_get_contents(filename, (char**)&contents, &length, NULL)){
        gunzip(contents, length);
        free(contents);
    }
    return 0;
}
#endif // defined TEST_GZIP

