/*
    ASCII renderer inspired by the terminal output demos in externals (e.g. pingo).
    Renders a rotating cube into a tiny framebuffer and prints characters to stdout.
    Uses only the public bootleg3d API, no SDL required.
*/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "b3d.h"

static const char *palette = " .:-=+*#%@";

static void render_frame(uint32_t *pixels, b3d_depth_t *depth, int w, int h, float t)
{
    b3d_init(pixels, depth, w, h, 65.0f);
    b3d_set_camera(0.0f, 0.0f, -2.3f, 0.0f, 0.0f, 0.0f);
    b3d_clear();

    b3d_reset();
    b3d_rotate_y(t * 0.7f);
    b3d_rotate_x(t * 0.5f);
    b3d_rotate_z(t * 0.3f);

    /* Colourful cube faces for contrast */
    b3d_triangle(-0.5,-0.5,-0.5, -0.5, 0.5,-0.5,  0.5, 0.5,-0.5, 0xfcd0a1);
    b3d_triangle(-0.5,-0.5,-0.5,  0.5, 0.5,-0.5,  0.5,-0.5,-0.5, 0xb1b695);
    b3d_triangle( 0.5,-0.5,-0.5,  0.5, 0.5,-0.5,  0.5, 0.5, 0.5, 0x53917e);
    b3d_triangle( 0.5,-0.5,-0.5,  0.5, 0.5, 0.5,  0.5,-0.5, 0.5, 0x63535b);
    b3d_triangle( 0.5,-0.5, 0.5,  0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0x6d1a36);
    b3d_triangle( 0.5,-0.5, 0.5, -0.5, 0.5, 0.5, -0.5,-0.5, 0.5, 0xd4e09b);
    b3d_triangle(-0.5,-0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5,-0.5, 0xf6f4d2);
    b3d_triangle(-0.5,-0.5, 0.5, -0.5, 0.5,-0.5, -0.5,-0.5,-0.5, 0xcbdfbd);
    b3d_triangle(-0.5, 0.5,-0.5, -0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0xf19c79);
    b3d_triangle(-0.5, 0.5,-0.5,  0.5, 0.5, 0.5,  0.5, 0.5,-0.5, 0xa44a3f);
    b3d_triangle( 0.5,-0.5, 0.5, -0.5,-0.5, 0.5, -0.5,-0.5,-0.5, 0x5465ff);
    b3d_triangle( 0.5,-0.5, 0.5, -0.5,-0.5,-0.5,  0.5,-0.5,-0.5, 0x788bff);
}

static char color_to_char(uint32_t c)
{
    unsigned r = (c >> 16) & 0xff;
    unsigned g = (c >> 8) & 0xff;
    unsigned b = (c >> 0) & 0xff;
    float lum = (0.299f * r + 0.587f * g + 0.114f * b) / 255.0f;
    int idx = (int)(lum * (float)(strlen(palette) - 1) + 0.5f);
    if (idx < 0) idx = 0;
    size_t len = strlen(palette);
    if ((size_t)idx >= len) idx = (int)len - 1;
    return palette[idx];
}

int main(int argc, char **argv)
{
    (void)argc; (void)argv;
    const int width = 96;
    const int height = 64;

    uint32_t *pixels = malloc((size_t)width * (size_t)height * sizeof(pixels[0]));
    b3d_depth_t *depth = malloc((size_t)width * (size_t)height * sizeof(depth[0]));
    if (!pixels || !depth) {
        fprintf(stderr, "Allocation failed\n");
        return 1;
    }

    float t = (float)clock() / (float)CLOCKS_PER_SEC;
    render_frame(pixels, depth, width, height, t);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            uint32_t c = pixels[(size_t)y * (size_t)width + (size_t)x];
            char ch = c ? color_to_char(c) : ' ';
            fputc(ch, stdout);
        }
        fputc('\n', stdout);
    }

    free(pixels);
    free(depth);
    return 0;
}
