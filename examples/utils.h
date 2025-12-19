/* Common utilities */

#ifndef UTILS_H
#define UTILS_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pngwrite.h"

/* Get snapshot path from command line arguments or environment variable
 *
 * Checks for B3D_SNAPSHOT environment variable first, then looks for
 * --snapshot=PATH argument.
 */
static inline const char *get_snapshot_path(int argc, char **argv)
{
    const char *env = getenv("B3D_SNAPSHOT");
    if (env && env[0])
        return env;

    for (int i = 1; i < argc; ++i) {
        const char *flag = "--snapshot=";
        size_t len = strlen(flag);
        if (!strncmp(argv[i], flag, len))
            return argv[i] + (int) len;
    }
    return NULL;
}

/* Write PNG image to file
 *
 * Converts ARGB image data to RGBA format and writes as PNG.
 */
static inline void write_png(const char *path,
                     const uint32_t *rgba,
                     int width,
                     int height)
{
    FILE *file = fopen(path, "wb");
    if (!file)
        return;

    uint8_t *out = malloc((size_t) width * (size_t) height * 4);
    if (!out) {
        fclose(file);
        return;
    }
    size_t idx = 0;
    for (int i = 0; i < width * height; ++i) {
        uint32_t p = rgba[i];
        out[idx++] = (uint8_t) ((p >> 16) & 0xff);
        out[idx++] = (uint8_t) ((p >> 8) & 0xff);
        out[idx++] = (uint8_t) (p & 0xff);
        out[idx++] = 0xff;
    }
    png_write(file, (unsigned) width, (unsigned) height, out, true);
    free(out);
    fclose(file);
}

#endif /* UTILS_H */
