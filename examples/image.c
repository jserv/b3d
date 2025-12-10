/*
    Renders a cube into a .tga file.
*/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "b3d.h"

void write_tga(char * path, void * rgba, int width, int height) {
    FILE * file = fopen(path, "wb");
    if (!file) return;
    fwrite("\x00\x00\x02\x00\x00\x00\x00\x00\x00\x00\x00\x00", 12, 1, file);
    fwrite(&width, 2, 1, file);
    fwrite(&height, 2, 1, file);
    int bpc = 4 * 8 | 0x2000;
    fwrite(&bpc, 2, 1, file);
    char * s = (char *)rgba;
    for (int i = 0; i < width * height; ++i) {
        fputc(s[0], file);
        fputc(s[1], file);
        fputc(s[2], file);
        fputc(0xff, file);
        s += 4;
    }
    fclose(file);
}

int main() {
    int width = 640;
    int height = 480;
    uint32_t * pixels = malloc(width * height * sizeof(pixels[0]));
    b3d_depth_t * depth = malloc(width * height * sizeof(depth[0]));
    b3d_init(pixels, depth, width, height, 70.0);
    b3d_set_camera(0, 0, -1.5, 0, 0, 0);
    b3d_rotate_x(1.5);
    b3d_rotate_y(2.6);
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
    write_tga("assets/cube.tga", pixels, width, height);
    free(pixels);
    free(depth);
    return 0;
}
