#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/param.h>

#define	NBBY	8 // число бит в байте
#define	MATCH_BITS	6
#define	MATCH_MIN	3
#define	MATCH_MAX	((1 << MATCH_BITS) + (MATCH_MIN - 1))
#define	OFFSET_MASK	((1 << (16 - MATCH_BITS)) - 1)
#define	LEMPEL_SIZE	1024

/*ARGSUSED*/
size_t
lzjb_compress(void *s_start, void *d_start, size_t s_len, size_t d_len, int n)
{
	uint8_t *src = s_start;
	uint8_t *dst = d_start;
	uint8_t *cpy;
	uint8_t *copymap = NULL;
	int copymask = 1 << (NBBY - 1);
	int mlen, offset, hash;
	uint16_t *hp;
	uint16_t lempel[LEMPEL_SIZE] = { 0 };

	while (src < (uint8_t *)s_start + s_len) {
		if ((copymask <<= 1) == (1 << NBBY)) {
			if (dst >= (uint8_t *)d_start + d_len - 1 - 2 * NBBY)
				return (s_len);
			copymask = 1;
			copymap = dst;
			*dst++ = 0;
		}
		if (src > (uint8_t *)s_start + s_len - MATCH_MAX) {
			*dst++ = *src++;
			continue;
		}
		hash = (src[0] << 16) + (src[1] << 8) + src[2];
		hash += hash >> 9;
		hash += hash >> 5;
		hp = &lempel[hash & (LEMPEL_SIZE - 1)];
		offset = (intptr_t)(src - *hp) & OFFSET_MASK;
		*hp = (uint16_t)(uintptr_t)src;
		cpy = src - offset;
		if (cpy >= (uint8_t *)s_start && cpy != src &&
		    src[0] == cpy[0] && src[1] == cpy[1] && src[2] == cpy[2]) {
			*copymap |= copymask;
			//printf("copymap = %02X\n", *copymap);
			for (mlen = MATCH_MIN; mlen < MATCH_MAX; mlen++)
				if (src[mlen] != cpy[mlen])
					break;
			*dst++ = ((mlen - MATCH_MIN) << (NBBY - MATCH_BITS)) |
			    (offset >> NBBY);
			*dst++ = (uint8_t)offset;
			src += mlen;
		} else {
			*dst++ = *src++;
		}
	}
	return (dst - (uint8_t *)d_start);
}

/*ARGSUSED*/
int
lzjb_decompress(void *s_start, void *d_start, size_t s_len, size_t d_len, int n)
{
	uint8_t *src = s_start;
	uint8_t *dst = d_start;
	uint8_t *d_end = (uint8_t *)d_start + d_len;
	uint8_t *cpy;
	uint8_t copymap = 0;
	int copymask = 1 << (NBBY - 1);

	while (dst < d_end) {
		if ((copymask <<= 1) == (1 << NBBY)) {
			copymask = 1;
			copymap = *src++;
		}
		if (copymap & copymask) {
			int mlen = (src[0] >> (NBBY - MATCH_BITS)) + MATCH_MIN;
			int offset = ((src[0] << NBBY) | src[1]) & OFFSET_MASK;
			src += 2;
			if ((cpy = dst - offset) < (uint8_t *)d_start)
				return (-1);
			while (--mlen >= 0 && dst < d_end)
				*dst++ = *cpy++;
		} else {
			*dst++ = *src++;
		}
	}
	return (0);
}
