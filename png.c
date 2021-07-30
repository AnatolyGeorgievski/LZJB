/*! \brief Реализация распаковки графики PNG
    [RFC 2083] PNG: Portable Network Graphics, March 1997
    https://datatracker.ietf.org/doc/html/rfc2083
*/
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
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

#define ADLER_PRIME 65521// (1<<16)-(1<<4)+(1)
#define ADLER_MASK 0xFFFFUL;
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
#define mod65521(x) ((x) - 0xFFF1*(unsigned long)(((x)*0x80078071ULL)>>47))
/*! Евгенская версия алгоритма, в цикле операций меньше.
*/
uint32_t ADLER32_update(unsigned long adler, uint8_t *p, size_t len)
{
    const size_t Nmax = 5552;// Евгенская магия
	unsigned long  s1 = (adler      ) & ADLER_MASK;
	unsigned long  s2 = (adler >> 16) & ADLER_MASK;
    while (len)
    {
        size_t n = len > Nmax ? Nmax : len;
        len -= n;
        do {
            s1 += *p++;// s1 + p0 + p1 + p2 + p3 ... MAX(s1) = 0xFFF0 + 0xFF*n
            s2 += s1;// s2 + s1*n + p0*(n) + p1*(n-1) ... MAX(s2) = (n+1)*(0xFFF0 + 0xFF*n/2)
        }while (--n);

        s1 = mod65521(s1);
        s2 = mod65521(s2);
    }
    return (s2 << 16) | s1;
}
static uint32_t crc_from_block(uint8_t *src, size_t len)
{
    int i;
    uint32_t crc = ~0UL;
    for(i=0; i<len; i++)
        crc = CRC32B_update(crc, src[i]);
    return ~crc;
}
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
extern uint8_t* deflate(uint8_t* dst, uint8_t* src, size_t s_len);
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
        if (crc != crc_from_block(cdata-4, length+4)) return -1;
        if (1) printf ("%-.*s crc=%08X  len=%d\n", 4, ctype, crc, length);
        if (strncasecmp(ctype,"ihdr", 4)==0) {
            hdr = (void*)cdata;
            printf ("\twidth=%d height=%d bits=%d colors=%d compression=%d\n", BE32(hdr->width), BE32(hdr->height),
                    hdr->bit_depth, hdr->color_type, hdr->compression);
            dst = malloc(BE32(hdr->width)*BE32(hdr->height)*4+256);// кто то пишет мимо
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
                int d_len = deflate(dst, cdata+2, length-6)-dst;
                printf("\ncompression =%1.2f%%\n", (float)(length-6)*100.f/(d_len));
                if(ADLER32_update(1, dst, d_len)==BE32(*(uint32_t *)(cdata+length-4)))
                    printf("ADLER32 Check sum ..ok\n");
            }
        } else
        if (strncasecmp(ctype,"iend", 4)==0) {
            break;
        }
    }
    printf("=done\n");
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
    char* filename = "test1.png";
    uint8_t *contents=NULL;
    size_t length=0;
    _get_contents(filename, (char**)&contents, &length, NULL);
    png_to_image(contents, length);
    free(contents);
    return 0;
}
#endif // defined

