/*
https://github.com/lz4/lz4/blob/dev/doc/lz4_Block_format.md
https://github.com/jibsen/blz4/blob/master/lz4_leparse.h

*/
#include <stdint.h>
#define LE16(x) (x)
static inline void lz4_memcpy(uint8_t *dst, uint8_t *src, uint32_t len) {
    __builtin_memcpy(dst,src, len);
}
static uint8_t * lsic_decode(uint8_t* src, uint32_t*t1)
{
    uint32_t n=*t1;
    uint32_t e;
    do{
        e = *src++;
        n+=e;
    } while(e==0xFF);
    *t1 = n;
    return src;
}
/*! \brief
    \param token - 4 бита кодирования токена */
static uint8_t * lsic_encode(uint8_t* dst, uint32_t n, uint32_t *token)
{
    if(n>=15) {
        while(n>=255+15){
            *dst++ = 255;
            n-=255;
        }
        *dst++ = n-15;
        n=15;
    }
    *token = n;
    return dst;
}
uint8_t* lz4_decode(uint8_t* dst, uint8_t* src, size_t slen)
{
    uint8_t *s_end = src+slen;
    while(1) {
        uint32_t token = *src++;
        uint32_t t1 = (token>> 4);// literal length
        if(t1 != 0) {
            if(t1==15) src = lsic_decode(src, &t1);
            lz4_memcpy(dst, src, t1);// copy literal
            dst+=t1, src+=t1;
        }
        /* Check for last incomplete sequence */
        if (src==s_end) break;
        uint32_t t2 = (token&0xF);// rle length
        //if(t2 != 0)
        {
            uint32_t offset = LE16(*(uint16_t *)src);
            src+=2;
            if(t2==15) src = lsic_decode(src, &t2);
            t2+=4;
            if (offset<t2){

            } else
                lz4_memcpy(dst, dst-offset, t2);// смешение может быть меньше длины
            dst+=t2;
        }
    }
    return dst;
}
/*! \brief декодирование LZF

http://dist.schmorp.de/liblzf/
*/
uint8_t* lzf_decode(uint8_t* dst, uint8_t* src, size_t slen)
{
    uint8_t *s_end = src+slen;
    while(src<s_end) {
        uint32_t token = *src++;
        if (token<(1<<5)) {// literal run 1..32
            uint32_t nlit = token+1;
            lz4_memcpy(dst, src, nlit);
        } else {// два формата: 3-5 и 8-13
            uint32_t mlen = token>>5;
            uint32_t offset = ((token&0x1F)<<8)+1;
            if(mlen==7) {
                mlen+= *src++;
            }
            offset += *src++; //+1
            mlen+=2;
            lz4_memcpy(dst, dst-offset, mlen);// возможно перекрытие
            dst+=mlen;
        }
    }
    return dst;
}
struct _Smap {
    uint32_t nlit:16;// длина литерала
    uint32_t moffset:11;// смещение
    uint32_t mlen:5;// длина кода-4
};
/*! \brief кодирование отделено от разбора */
uint8_t* lz4_encode(uint8_t* dst, uint8_t* src, size_t slen, struct _Smap* map)
{
    uint32_t t1, t2;
    int i=0;
    while(1) {
        uint8_t *token = dst++;
        uint32_t nlit = map->nlit; // длина литерала
        dst = lsic_encode(dst, nlit, &t1);
        lz4_memcpy(dst, src+i, nlit);
        dst+=nlit;
        i  +=nlit;
        if (i==slen) {
            *token = t1<<4;
            break;
        }
        uint32_t offset = map->moffset;
        *(uint16_t*)dst = LE16(offset);
        t2  = map->mlen;//-4; уже вычли в процессе разбора
        dst = lsic_encode(dst, t2, &t2);
        *token = (t1<<4) | t2;
    }
    return dst;
}
/* LZ4_LEGACY_MAGIC (LE32)? |  hdr_packedsize (LE32) | lz4_decode */
