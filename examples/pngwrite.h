/*
 * Minimal PNG writer.
 */

#ifndef PNG_WRITE_H
#define PNG_WRITE_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

/* Write a byte */
#define PNG_WRITE_PUT(u) fputc((u), fp)

static void png_write(FILE *fp,
                      unsigned w,
                      unsigned h,
                      const uint8_t *img,
                      bool alpha)
{
    static const unsigned t[] = {
        0,          0x1db71064, 0x3b6e20c8, 0x26d930ac, 0x76dc4190, 0x6b6b51f4,
        0x4db26158, 0x5005713c, 0xedb88320, 0xf00f9344, 0xd6d6a3e8, 0xcb61b38c,
        0x9b64c2b0, 0x86d3d2d4, 0xa00ae278, 0xbdbdf21c,
    };
    unsigned a = 1, b = 0, c, p = w * (alpha ? 4 : 3) + 1, x, y;

#define PNG_WRITE_U8A(ua, l)         \
    for (size_t i = 0; i < (l); i++) \
        PNG_WRITE_PUT((ua)[i]);

#define PNG_WRITE_U32(u)                  \
    do {                                  \
        PNG_WRITE_PUT((u) >> 24);         \
        PNG_WRITE_PUT(((u) >> 16) & 255); \
        PNG_WRITE_PUT(((u) >> 8) & 255);  \
        PNG_WRITE_PUT((u) & 255);         \
    } while (0)

#define PNG_WRITE_U8C(u)          \
    do {                          \
        PNG_WRITE_PUT(u);         \
        c ^= (u);                 \
        c = (c >> 4) ^ t[c & 15]; \
        c = (c >> 4) ^ t[c & 15]; \
    } while (0)

#define PNG_WRITE_U8AC(ua, l)        \
    for (size_t i = 0; i < (l); i++) \
    PNG_WRITE_U8C((ua)[i])

#define PNG_WRITE_U16LC(u)               \
    do {                                 \
        PNG_WRITE_U8C((u) & 255);        \
        PNG_WRITE_U8C(((u) >> 8) & 255); \
    } while (0)

#define PNG_WRITE_U32C(u)                 \
    do {                                  \
        PNG_WRITE_U8C((u) >> 24);         \
        PNG_WRITE_U8C(((u) >> 16) & 255); \
        PNG_WRITE_U8C(((u) >> 8) & 255);  \
        PNG_WRITE_U8C((u) & 255);         \
    } while (0)

#define PNG_WRITE_U8ADLER(u)   \
    do {                       \
        PNG_WRITE_U8C(u);      \
        a = (a + (u)) % 65521; \
        b = (b + a) % 65521;   \
    } while (0)

#define PNG_WRITE_BEGIN(s, l) \
    do {                      \
        PNG_WRITE_U32(l);     \
        c = ~0U;              \
        PNG_WRITE_U8AC(s, 4); \
    } while (0)

#define PNG_WRITE_END() PNG_WRITE_U32(~c)

    PNG_WRITE_U8A("\x89PNG\r\n\32\n", 8); /* Magic */
    PNG_WRITE_BEGIN("IHDR", 13);          /* IHDR chunk { */
    PNG_WRITE_U32C(w);
    PNG_WRITE_U32C(h); /*   Width & Height (8 bytes) */
    PNG_WRITE_U8C(8);
    /* Depth=8, Color=True color with/without alpha (2 bytes) */
    PNG_WRITE_U8C(alpha ? 6 : 2);
    /* Compression=Deflate, Filter=No, Interlace=No (3 bytes) */
    PNG_WRITE_U8AC("\0\0\0", 3);
    PNG_WRITE_END();                              /* } */
    PNG_WRITE_BEGIN("IDAT", 2 + h * (5 + p) + 4); /* IDAT chunk { */
    PNG_WRITE_U8AC("\x78\1", 2); /*   Deflate block begin (2 bytes) */
    /*   Each horizontal line makes a block for simplicity */
    for (y = 0; y < h; y++) {
        /* 1 for the last block, 0 for others (1 byte) */
        PNG_WRITE_U8C(y == h - 1);
        PNG_WRITE_U16LC(p);
        /* Size of block in little endian and its 1's complement (4 bytes) */
        PNG_WRITE_U16LC(~p);
        PNG_WRITE_U8ADLER(0); /*   No filter prefix (1 byte) */
        for (x = 0; x < p - 1; x++, img++)
            PNG_WRITE_U8ADLER(*img); /*   Image pixel data */
    }
    PNG_WRITE_U32C((b << 16) |
                   a); /*   Deflate block end with adler (4 bytes) */
    PNG_WRITE_END();   /* } */
    PNG_WRITE_BEGIN("IEND", 0);
    PNG_WRITE_END(); /* IEND chunk {} */
}

#undef PNG_WRITE_PUT
#endif /* PNG_WRITE_H */
