/*
 * Microcraft - Voxel terrain demo using B3D renderer
 *
 * A Minecraft-style voxel renderer demonstrating procedural terrain generation
 * with Perlin noise, first-person camera controls, and block texturing via
 * tessellated faces with procedural color variation.
 *
 * Controls:
 *   WASD / Arrow keys - Move
 *   Mouse            - Look around
 *   Space            - Move up
 *   Shift/C          - Move down
 *   ESC              - Quit
 *
 * Supports headless snapshots with --snapshot=PATH or B3D_SNAPSHOT.
 */

#include <SDL.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "b3d-math.h"
#include "b3d.h"
#include "utils.h"

/* World dimensions */
#define MAP_WIDTH 64
#define MAP_DEPTH 64
#define MAP_HEIGHT 64

/* Terrain generation seed */
#define TERRAIN_SEED 1

/* Render distance (blocks from player) */
#define RENDER_DISTANCE 20.0f
#define RENDER_DISTANCE_SQ (RENDER_DISTANCE * RENDER_DISTANCE)

/* Water level */
#define WATER_LEVEL 28

/* Tessellation level for texture simulation (NxN grid per face) */
#define TESS_LEVEL 4
#define TESS_SIZE (1.0f / TESS_LEVEL)

/* Distance-based tessellation thresholds */
#define TESS_NEAR_DIST_SQ 64.0f /* Full tessellation within 8 blocks */
#define TESS_MID_DIST_SQ 196.0f /* Half tessellation within 14 blocks */

/* Terrain generation frequencies */
#define TERRAIN_BASE_FREQ 0.03f
#define TERRAIN_HILLS_FREQ 0.08f
#define TERRAIN_DETAIL_FREQ 0.15f
#define CAVE_FREQ 0.0625f

/* Terrain amplitude */
#define TERRAIN_BASE_AMP 20.0f
#define TERRAIN_HILLS_AMP 12.0f
#define TERRAIN_DETAIL_AMP 6.0f

/* Mouse sensitivity */
#define MOUSE_SENSITIVITY 0.003f

/* Color manipulation macros */
#define COLOR_R(c) (((c) >> 16) & 0xFF)
#define COLOR_G(c) (((c) >> 8) & 0xFF)
#define COLOR_B(c) ((c) & 0xFF)
#define COLOR_RGB(r, g, b) \
    (((uint32_t) (r) << 16) | ((uint32_t) (g) << 8) | (uint32_t) (b))
#define CLAMP_BYTE(x) ((x) < 0 ? 0 : ((x) > 255 ? 255 : (x)))

/* Block types */
typedef enum {
    BLOCK_AIR = 0,
    BLOCK_GRASS,
    BLOCK_DIRT,
    BLOCK_SAND,
    BLOCK_STONE,
    BLOCK_GRAVEL,
    BLOCK_WOOD,
    BLOCK_LEAVES,
    BLOCK_COBBLESTONE,
    BLOCK_WATER,
    BLOCK_TALL_GRASS,
    NUMBER_OF_BLOCKS,
} block_type_t;

/* Block base colors (0xRRGGBB format) - Classic Minecraft palette */
static const uint32_t block_colors[] = {
    [BLOCK_AIR] = 0x000000,
    [BLOCK_GRASS] = 0x79C05A,  /* Classic grass green (muted yellow-green) */
    [BLOCK_DIRT] = 0x8B5A2B,   /* Classic Minecraft dirt brown */
    [BLOCK_SAND] = 0xDCCFA5,   /* Classic sand (warm beige) */
    [BLOCK_STONE] = 0x7D7D7D,  /* Classic stone gray */
    [BLOCK_GRAVEL] = 0x8A8A8A, /* Classic gravel */
    [BLOCK_WOOD] = 0x6A5030,   /* Classic oak bark */
    [BLOCK_LEAVES] = 0x4CA031, /* Classic oak leaves (darker green) */
    [BLOCK_COBBLESTONE] = 0x787878,
    [BLOCK_WATER] = 0x3366CC, /* Classic water blue */
    [BLOCK_TALL_GRASS] = 0x79C05A,
};

/* Grass block face colors (top green, sides/bottom dirt) */
static const uint32_t grass_top_color = 0x79C05A;  /* Classic grass green */
static const uint32_t grass_side_color = 0x8B5A2B; /* dirt brown */

/* Wood face colors (top/bottom show rings, sides show bark) */
static const uint32_t wood_top_color = 0xB5945A;  /* Classic oak rings */
static const uint32_t wood_side_color = 0x6A5030; /* Classic oak bark */

/* World map storage */
static block_type_t world_map[MAP_WIDTH][MAP_DEPTH][MAP_HEIGHT];

/* Cached terrain heights */
static uint8_t height_cache[MAP_WIDTH][MAP_DEPTH];

/* Player state */
static float player_x, player_y, player_z;
static float player_yaw, player_pitch;

/* World time for day/night cycle */
static uint32_t world_time = 0;

/* Perlin noise hash table */
static const uint8_t perlin_hash[256] = {
    208, 34,  231, 213, 32,  248, 233, 56,  161, 78,  24,  140, 71,  48,  140,
    254, 245, 255, 247, 247, 40,  185, 248, 251, 245, 28,  124, 204, 204, 76,
    36,  1,   107, 28,  234, 163, 202, 224, 245, 128, 167, 204, 9,   92,  217,
    54,  239, 174, 173, 102, 193, 189, 190, 121, 100, 108, 167, 44,  43,  77,
    180, 204, 8,   81,  70,  223, 11,  38,  24,  254, 210, 210, 177, 32,  81,
    195, 243, 125, 8,   169, 112, 32,  97,  53,  195, 13,  203, 9,   47,  104,
    125, 117, 114, 124, 165, 203, 181, 235, 193, 206, 70,  180, 174, 0,   167,
    181, 41,  164, 30,  116, 127, 198, 245, 146, 87,  224, 149, 206, 57,  4,
    192, 210, 65,  210, 129, 240, 178, 105, 228, 108, 245, 148, 140, 40,  35,
    195, 38,  58,  65,  207, 215, 253, 65,  85,  208, 76,  62,  3,   237, 55,
    89,  232, 50,  217, 64,  244, 157, 199, 121, 252, 90,  17,  212, 203, 149,
    152, 140, 187, 234, 177, 73,  174, 193, 100, 192, 143, 97,  53,  145, 135,
    19,  103, 13,  90,  135, 151, 199, 91,  239, 247, 33,  39,  145, 101, 120,
    99,  3,   186, 86,  99,  41,  237, 203, 111, 79,  220, 135, 158, 42,  30,
    154, 120, 67,  87,  167, 135, 176, 183, 191, 253, 115, 184, 21,  233, 58,
    129, 233, 142, 39,  128, 211, 118, 137, 139, 255, 114, 20,  218, 113, 154,
    27,  127, 246, 250, 1,   8,   198, 250, 209, 92,  222, 173, 21,  88,  102,
    219,
};

/* Fast hash-based noise for texture variation
 * Returns value in range [0, 255]
 */
static inline int noise_hash(int x, int y, int z, int seed)
{
    int h = perlin_hash[(x + seed) & 0xFF];
    h = perlin_hash[(h + y) & 0xFF];
    h = perlin_hash[(h + z) & 0xFF];
    return h;
}

/* Apply noise variation to a base color
 * Returns modified color with brightness variation
 */
static uint32_t apply_noise(uint32_t base_color, int x, int y, int z, int scale)
{
    int noise = noise_hash(x, y, z, TERRAIN_SEED);
    int variation = (noise % scale) - scale / 2;

    int r = CLAMP_BYTE((int) COLOR_R(base_color) + variation);
    int g = CLAMP_BYTE((int) COLOR_G(base_color) + variation);
    int b = CLAMP_BYTE((int) COLOR_B(base_color) + variation);

    return COLOR_RGB(r, g, b);
}

/* Perlin noise functions */
static inline int perlin_noise2(int x, int y, int seed)
{
    int tmp = perlin_hash[(y + seed) & 0xFF];
    return perlin_hash[(tmp + x) & 0xFF];
}

static inline float perlin_lerp(float x, float y, float s)
{
    return x + s * s * (3.0f - 2.0f * s) * (y - x);
}

static float perlin2d(int seed, float x, float y, float freq)
{
    float xa = x * freq, ya = y * freq;
    float amp = 1.0f;
    float fin = 0.0f;
    float div = 0.0f;

    for (int i = 0; i < 4; i++) {
        div += 256.0f * amp;
        int x_int = (int) xa, y_int = (int) ya;
        float x_frac = xa - (float) x_int, y_frac = ya - (float) y_int;

        int s = perlin_noise2(x_int, y_int, seed);
        int t = perlin_noise2(x_int + 1, y_int, seed);
        int u = perlin_noise2(x_int, y_int + 1, seed);
        int v = perlin_noise2(x_int + 1, y_int + 1, seed);

        float low = perlin_lerp((float) s, (float) t, x_frac);
        float high = perlin_lerp((float) u, (float) v, x_frac);

        fin += perlin_lerp(low, high, y_frac) * amp;
        amp *= 0.5f;
        xa *= 2.0f, ya *= 2.0f;
    }
    return fin / div;
}

/* World random for structure placement */
static unsigned int world_rand_state;

static void world_srand(unsigned int seed)
{
    world_rand_state = seed;
}

static int world_rand(int max)
{
    world_rand_state = world_rand_state * 1103515245 + 12345;
    return (int) ((world_rand_state >> 16) & 0x7FFF) % max;
}

/* Check if position is valid for tree placement */
static int can_place_tree(int x, int z, int surface_y)
{
    if (x < 3 || x >= MAP_WIDTH - 3 || z < 3 || z >= MAP_DEPTH - 3)
        return 0;

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dz = -1; dz <= 1; ++dz) {
            if (dx == 0 && dz == 0)
                continue;
            int neighbor_height = height_cache[x + dx][z + dz];
            if (abs(neighbor_height - surface_y) > 2)
                return 0;
            if (neighbor_height <= WATER_LEVEL)
                return 0;
        }
    }
    return 1;
}

/* Generate a tree at position */
static void gen_tree(int x, int surface_y, int z)
{
    int trunk_height = 4 + world_rand(3);
    for (int i = 1; i <= trunk_height; ++i) {
        int ty = surface_y + i;
        if (ty >= 0 && ty < MAP_HEIGHT)
            world_map[x][z][ty] = BLOCK_WOOD;
    }

    int top = surface_y + trunk_height;

    /* Base leaves 5x2x5 */
    for (int dx = -2; dx <= 2; ++dx) {
        for (int dz = -2; dz <= 2; ++dz) {
            int lx = x + dx, lz = z + dz;
            if (lx < 0 || lx >= MAP_WIDTH || lz < 0 || lz >= MAP_DEPTH)
                continue;
            for (int dy = 0; dy < 2; ++dy) {
                int ly = top + dy;
                if (ly >= 0 && ly < MAP_HEIGHT &&
                    world_map[lx][lz][ly] == BLOCK_AIR)
                    world_map[lx][lz][ly] = BLOCK_LEAVES;
            }
        }
    }

    /* Top leaves 3x2x3 */
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dz = -1; dz <= 1; ++dz) {
            int lx = x + dx, lz = z + dz;
            if (lx < 0 || lx >= MAP_WIDTH || lz < 0 || lz >= MAP_DEPTH)
                continue;
            int ly = top + 2;
            if (ly >= 0 && ly < MAP_HEIGHT &&
                world_map[lx][lz][ly] == BLOCK_AIR)
                world_map[lx][lz][ly] = BLOCK_LEAVES;
        }
    }
}

/* Check if block type is transparent */
static inline int is_transparent(block_type_t block)
{
    return block == BLOCK_AIR || block == BLOCK_WATER ||
           block == BLOCK_LEAVES || block == BLOCK_TALL_GRASS;
}

/* Initialize the world with procedural terrain */
static void world_init(void)
{
    /* Pass 1: Generate heightmap and base terrain */
    for (int x = 0; x < MAP_WIDTH; ++x) {
        for (int z = 0; z < MAP_DEPTH; ++z) {
            float base = perlin2d(TERRAIN_SEED, (float) (x + 0xFFFFFF),
                                  (float) (z + 0xFFFFFF), 0.03f) *
                         20.0f;
            float hills = perlin2d(TERRAIN_SEED + 1, (float) (x + 0xFFFFFF),
                                   (float) (z + 0xFFFFFF), 0.08f) *
                          12.0f;
            float detail = perlin2d(TERRAIN_SEED + 2, (float) (x + 0xFFFFFF),
                                    (float) (z + 0xFFFFFF), 0.15f) *
                           6.0f;

            int height = (int) (base + hills + detail) + WATER_LEVEL - 12;
            if (height < 2)
                height = 2;
            if (height >= MAP_HEIGHT - 8)
                height = MAP_HEIGHT - 9;

            height_cache[x][z] = (uint8_t) height;

            for (int y = 0; y < MAP_HEIGHT; ++y) {
                if (y > height) {
                    if (y <= WATER_LEVEL)
                        world_map[x][z][y] = BLOCK_WATER;
                    else
                        world_map[x][z][y] = BLOCK_AIR;
                } else if (y == height) {
                    if (height <= WATER_LEVEL + 2)
                        world_map[x][z][y] = BLOCK_SAND;
                    else
                        world_map[x][z][y] = BLOCK_GRASS;
                } else if (y > height - 4) {
                    if (height <= WATER_LEVEL + 2)
                        world_map[x][z][y] = BLOCK_SAND;
                    else
                        world_map[x][z][y] = BLOCK_DIRT;
                } else {
                    world_map[x][z][y] = BLOCK_STONE;
                }
            }
            world_map[x][z][0] = BLOCK_STONE;
        }
    }

    /* Pass 2: Cave generation */
    for (int x = 0; x < MAP_WIDTH; ++x) {
        for (int z = 0; z < MAP_DEPTH; ++z) {
            float noise_point =
                perlin2d(TERRAIN_SEED + 10, (float) (x + 0xFFFFFF),
                         (float) (z + 0xFFFFFF), 0.0625f);

            if (noise_point >= 0.47f && noise_point <= 0.53f) {
                int surface = height_cache[x][z];
                int elevation =
                    (int) (perlin2d(TERRAIN_SEED + 12, (float) (x + 0xFFFFFF),
                                    (float) (z + 0xFFFFFF), 0.0625f) *
                           6.0f);
                int cave_height =
                    (int) (perlin2d(TERRAIN_SEED + 13, (float) (x + 0xFFFFFF),
                                    (float) (z + 0xFFFFFF), 0.0625f) *
                               3.0f +
                           2.0f);

                int cave_floor = surface - 8 - elevation;
                int cave_ceiling = cave_floor + cave_height;
                if (cave_floor < 2)
                    cave_floor = 2;

                for (int y = cave_floor; y < cave_ceiling && y < surface - 4;
                     ++y) {
                    if (y > 0 && y < MAP_HEIGHT &&
                        world_map[x][z][y] != BLOCK_WATER)
                        world_map[x][z][y] = BLOCK_AIR;
                }

                if (cave_floor > 1 && cave_floor < MAP_HEIGHT &&
                    world_map[x][z][cave_floor - 1] == BLOCK_STONE)
                    world_map[x][z][cave_floor - 1] = BLOCK_GRAVEL;
            }
        }
    }

    /* Pass 3: Generate trees */
    world_srand(TERRAIN_SEED * 12345);
    int tree_count = 64 + world_rand(16);
    int trees_placed = 0;
    for (int attempts = 0;
         attempts < tree_count * 3 && trees_placed < tree_count; ++attempts) {
        int x = world_rand(MAP_WIDTH);
        int z = world_rand(MAP_DEPTH);
        int height = height_cache[x][z];

        while (height > 0 && world_map[x][z][height] == BLOCK_AIR)
            height--;

        if (world_map[x][z][height] == BLOCK_GRASS && height < MAP_HEIGHT - 8 &&
            can_place_tree(x, z, height)) {
            gen_tree(x, height, z);
            trees_placed++;
        }
    }

    /* Pass 4: Plant tall grass */
    world_srand(TERRAIN_SEED * 54321);
    for (int x = 0; x < MAP_WIDTH; ++x) {
        for (int z = 0; z < MAP_DEPTH; ++z) {
            int height = height_cache[x][z];
            while (height > 0 && world_map[x][z][height] == BLOCK_AIR)
                height--;

            if (world_map[x][z][height] == BLOCK_GRASS && world_rand(2) == 0) {
                if (height + 1 < MAP_HEIGHT &&
                    world_map[x][z][height + 1] == BLOCK_AIR)
                    world_map[x][z][height + 1] = BLOCK_TALL_GRASS;
            }
        }
    }

    /* Pass 5: Update height_cache to include trees and vegetation */
    for (int x = 0; x < MAP_WIDTH; ++x) {
        for (int z = 0; z < MAP_DEPTH; ++z) {
            int max_h = height_cache[x][z];
            for (int y = max_h + 1; y < MAP_HEIGHT; ++y) {
                if (world_map[x][z][y] != BLOCK_AIR)
                    max_h = y;
            }
            height_cache[x][z] = (uint8_t) max_h;
        }
    }

    /* Set initial player position - view terrain with visible slopes */
    player_x = MAP_WIDTH / 3.0f;
    player_z = MAP_DEPTH / 3.0f;
    int center_height =
        height_cache[(int) (MAP_WIDTH / 3)][(int) (MAP_DEPTH / 3)];
    player_y =
        (float) (center_height + 3); /* Eye level slightly above ground */
    if (player_y < WATER_LEVEL + 2)
        player_y = WATER_LEVEL + 2;
    player_yaw = 0.8f;    /* Angle to see varied terrain */
    player_pitch = 0.15f; /* Slight downward look */
}

/* Face brightness multipliers (Minecraft-style directional shading)
 * Index: 0=-Z, 1=+Z, 2=+X, 3=-X, 4=+Y (top), 5=-Y (bottom)
 */
static const float face_brightness[6] = {
    0.8f, /* -Z face - 80% */
    0.8f, /* +Z face - 80% */
    0.6f, /* +X face - 60% */
    0.6f, /* -X face - 60% */
    1.0f, /* +Y top  - 100% */
    0.5f, /* -Y bottom - 50% */
};

/* Get the base color for a block face, with special handling for
 * grass (green top, dirt sides) and wood (rings top, bark sides)
 */
static uint32_t get_face_color(block_type_t block, int face)
{
    /* Face indices: 0=-Z, 1=+Z, 2=+X, 3=-X, 4=+Y (top), 5=-Y (bottom) */
    if (block == BLOCK_GRASS) {
        if (face == 4)
            return grass_top_color; /* Top is green */
        return grass_side_color;    /* Sides and bottom are dirt */
    }
    if (block == BLOCK_WOOD) {
        if (face == 4 || face == 5)
            return wood_top_color; /* Top/bottom show rings */
        return wood_side_color;    /* Sides show bark */
    }
    return block_colors[block];
}

/* Get color for grass side face with green edge at top (Minecraft style)
 * Only the topmost ~25% has green edge, rest is dirt
 */
static uint32_t get_grass_side_color(int ty, int tess)
{
    /* Minecraft: green edge only at very top of side face (~2-3 pixels of 16)
     * With TESS_LEVEL=4: only ty=3 (top row) is green
     * With tess=1 or 2: show dirt (green edge too small to see) */
    if (tess >= 3 && ty == tess - 1)
        return grass_top_color; /* Top row green only at high tessellation */
    return grass_side_color;    /* Rest is dirt */
}

/* Get noise scale for a block type (how much color variation)
 * Reduced values for cleaner Minecraft-style colors
 */
static int get_noise_scale(block_type_t block)
{
    switch (block) {
    case BLOCK_SAND:
        return 32; /* Subtle sand variation */
    case BLOCK_GRAVEL:
        return 80; /* Moderate gravel variation */
    case BLOCK_STONE:
        return 50; /* Subtle stone variation */
    case BLOCK_COBBLESTONE:
        return 60;
    case BLOCK_LEAVES:
        return 45; /* Subtle leaf variation */
    case BLOCK_GRASS:
        return 40; /* Clean grass colors */
    case BLOCK_DIRT:
        return 50; /* Subtle dirt variation */
    case BLOCK_WATER:
        return 25; /* Minimal water variation */
    case BLOCK_WOOD:
        return 45;
    default:
        return 50;
    }
}

/* Face normals for lighting */
static const float face_normals[6][3] = {
    {0.0f, 0.0f, -1.0f}, /* -Z */
    {0.0f, 0.0f, 1.0f},  /* +Z */
    {1.0f, 0.0f, 0.0f},  /* +X */
    {-1.0f, 0.0f, 0.0f}, /* -X */
    {0.0f, 1.0f, 0.0f},  /* +Y */
    {0.0f, -1.0f, 0.0f}, /* -Y */
};

/* Render a block face with distance-based tessellation
 * Near blocks get full tessellation, far blocks get simple quads
 * Vertices calculated directly in world space (no matrix stack)
 * @bx, by, bz: Block world coordinates
 * @face: Face index (0-5)
 * @block: Block type
 * @time_coef: Day/night brightness
 * @dist_sq: Squared distance from player (for LOD)
 */
static void render_face_tessellated(int bx,
                                    int by,
                                    int bz,
                                    int face,
                                    block_type_t block,
                                    float time_coef,
                                    float dist_sq)
{
    int noise_scale = get_noise_scale(block);
    float brightness = time_coef * face_brightness[face];

    /* Distance-based tessellation level */
    int tess = (dist_sq < TESS_NEAR_DIST_SQ)  ? TESS_LEVEL
               : (dist_sq < TESS_MID_DIST_SQ) ? 2
                                              : 1;
    float tess_size = 1.0f / (float) tess;

    int is_leaf = (block == BLOCK_LEAVES);
    int is_grass_side = (block == BLOCK_GRASS && face != 4 && face != 5);

    float nx = face_normals[face][0];
    float ny = face_normals[face][1];
    float nz = face_normals[face][2];

    /* Block center in world space */
    float cx = (float) bx, cy = (float) by, cz = (float) bz;

    /* Hoist base color computation outside tessellation loop */
    uint32_t face_base_color = is_grass_side ? 0 : get_face_color(block, face);

    for (int ty = 0; ty < tess; ty++) {
        /* For grass sides, color varies by row (green top, dirt bottom) */
        uint32_t row_base_color =
            is_grass_side ? get_grass_side_color(ty, tess) : face_base_color;

        for (int tx = 0; tx < tess; tx++) {
            /* Leaf transparency */
            if (is_leaf) {
                int skip = noise_hash(bx * 16 + tx, by * 16 + ty, bz + face,
                                      TERRAIN_SEED + 100);
                if (skip < 100)
                    continue;
            }

            /* UV coordinates within face */
            float u0 = tx * tess_size, v0 = ty * tess_size;
            float u1 = (tx + 1) * tess_size, v1 = (ty + 1) * tess_size;

            /* Apply noise and brightness (noise varies per quad) */
            int wx = bx * tess + tx;
            int wy = by * tess + ty;
            int wz = bz * tess + face;
            uint32_t color =
                apply_noise(row_base_color, wx, wy, wz, noise_scale);

            int r = (int) (COLOR_R(color) * brightness);
            int g = (int) (COLOR_G(color) * brightness);
            int b_col = (int) (COLOR_B(color) * brightness);
            color = COLOR_RGB(r, g, b_col);

            /* Calculate world-space vertices directly (no matrix stack) */
            b3d_tri_t tri1, tri2;
            float x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;

            switch (face) {
            case 0: /* -Z face */
                x0 = cx - 0.5f + u0;
                y0 = cy - 0.5f + v0;
                z0 = cz - 0.5f;
                x1 = cx - 0.5f + u0;
                y1 = cy - 0.5f + v1;
                z1 = cz - 0.5f;
                x2 = cx - 0.5f + u1;
                y2 = cy - 0.5f + v1;
                z2 = cz - 0.5f;
                x3 = cx - 0.5f + u1;
                y3 = cy - 0.5f + v0;
                z3 = cz - 0.5f;
                break;
            case 1: /* +Z face */
                x0 = cx - 0.5f + u1;
                y0 = cy - 0.5f + v0;
                z0 = cz + 0.5f;
                x1 = cx - 0.5f + u1;
                y1 = cy - 0.5f + v1;
                z1 = cz + 0.5f;
                x2 = cx - 0.5f + u0;
                y2 = cy - 0.5f + v1;
                z2 = cz + 0.5f;
                x3 = cx - 0.5f + u0;
                y3 = cy - 0.5f + v0;
                z3 = cz + 0.5f;
                break;
            case 2: /* +X face */
                x0 = cx + 0.5f;
                y0 = cy - 0.5f + v0;
                z0 = cz - 0.5f + u0;
                x1 = cx + 0.5f;
                y1 = cy - 0.5f + v1;
                z1 = cz - 0.5f + u0;
                x2 = cx + 0.5f;
                y2 = cy - 0.5f + v1;
                z2 = cz - 0.5f + u1;
                x3 = cx + 0.5f;
                y3 = cy - 0.5f + v0;
                z3 = cz - 0.5f + u1;
                break;
            case 3: /* -X face */
                x0 = cx - 0.5f;
                y0 = cy - 0.5f + v0;
                z0 = cz - 0.5f + u1;
                x1 = cx - 0.5f;
                y1 = cy - 0.5f + v1;
                z1 = cz - 0.5f + u1;
                x2 = cx - 0.5f;
                y2 = cy - 0.5f + v1;
                z2 = cz - 0.5f + u0;
                x3 = cx - 0.5f;
                y3 = cy - 0.5f + v0;
                z3 = cz - 0.5f + u0;
                break;
            case 4: /* +Y face (top) */
                x0 = cx - 0.5f + u0;
                y0 = cy + 0.5f;
                z0 = cz - 0.5f + v0;
                x1 = cx - 0.5f + u0;
                y1 = cy + 0.5f;
                z1 = cz - 0.5f + v1;
                x2 = cx - 0.5f + u1;
                y2 = cy + 0.5f;
                z2 = cz - 0.5f + v1;
                x3 = cx - 0.5f + u1;
                y3 = cy + 0.5f;
                z3 = cz - 0.5f + v0;
                break;
            case 5: /* -Y face (bottom) */
                x0 = cx - 0.5f + u0;
                y0 = cy - 0.5f;
                z0 = cz - 0.5f + v1;
                x1 = cx - 0.5f + u0;
                y1 = cy - 0.5f;
                z1 = cz - 0.5f + v0;
                x2 = cx - 0.5f + u1;
                y2 = cy - 0.5f;
                z2 = cz - 0.5f + v0;
                x3 = cx - 0.5f + u1;
                y3 = cy - 0.5f;
                z3 = cz - 0.5f + v1;
                break;
            default:
                continue;
            }

            tri1.v[0] = (b3d_point_t) {x0, y0, z0};
            tri1.v[1] = (b3d_point_t) {x1, y1, z1};
            tri1.v[2] = (b3d_point_t) {x2, y2, z2};

            tri2.v[0] = (b3d_point_t) {x0, y0, z0};
            tri2.v[1] = (b3d_point_t) {x2, y2, z2};
            tri2.v[2] = (b3d_point_t) {x3, y3, z3};

            b3d_triangle_lit(&tri1, nx, ny, nz, color);
            b3d_triangle_lit(&tri2, nx, ny, nz, color);
        }
    }
}

/* Render cross-shaped billboard for vegetation (tall grass)
 * Vertices calculated directly in world space
 */
static void render_cross_billboard(int bx,
                                   int by,
                                   int bz,
                                   block_type_t block,
                                   float time_coef)
{
    uint32_t base_color = block_colors[block];
    int noise_scale = get_noise_scale(block);
    uint32_t color = apply_noise(base_color, bx, by, bz, noise_scale);

    int r = (int) (COLOR_R(color) * time_coef);
    int g = (int) (COLOR_G(color) * time_coef);
    int b_val = (int) (COLOR_B(color) * time_coef);
    color = COLOR_RGB(r, g, b_val);

    float cx = (float) bx, cy = (float) by, cz = (float) bz;
    float h = 0.4f;
    b3d_tri_t tri;

    /* Diagonal 1 - front */
    tri.v[0] = (b3d_point_t) {cx - h, cy - 0.5f, cz - h};
    tri.v[1] = (b3d_point_t) {cx - h, cy + 0.5f, cz - h};
    tri.v[2] = (b3d_point_t) {cx + h, cy + 0.5f, cz + h};
    b3d_triangle_lit(&tri, 0.7f, 0.0f, 0.7f, color);

    tri.v[0] = (b3d_point_t) {cx - h, cy - 0.5f, cz - h};
    tri.v[1] = (b3d_point_t) {cx + h, cy + 0.5f, cz + h};
    tri.v[2] = (b3d_point_t) {cx + h, cy - 0.5f, cz + h};
    b3d_triangle_lit(&tri, 0.7f, 0.0f, 0.7f, color);

    /* Diagonal 1 - back */
    tri.v[0] = (b3d_point_t) {cx + h, cy - 0.5f, cz + h};
    tri.v[1] = (b3d_point_t) {cx + h, cy + 0.5f, cz + h};
    tri.v[2] = (b3d_point_t) {cx - h, cy + 0.5f, cz - h};
    b3d_triangle_lit(&tri, -0.7f, 0.0f, -0.7f, color);

    tri.v[0] = (b3d_point_t) {cx + h, cy - 0.5f, cz + h};
    tri.v[1] = (b3d_point_t) {cx - h, cy + 0.5f, cz - h};
    tri.v[2] = (b3d_point_t) {cx - h, cy - 0.5f, cz - h};
    b3d_triangle_lit(&tri, -0.7f, 0.0f, -0.7f, color);

    /* Diagonal 2 - front */
    tri.v[0] = (b3d_point_t) {cx - h, cy - 0.5f, cz + h};
    tri.v[1] = (b3d_point_t) {cx - h, cy + 0.5f, cz + h};
    tri.v[2] = (b3d_point_t) {cx + h, cy + 0.5f, cz - h};
    b3d_triangle_lit(&tri, 0.7f, 0.0f, -0.7f, color);

    tri.v[0] = (b3d_point_t) {cx - h, cy - 0.5f, cz + h};
    tri.v[1] = (b3d_point_t) {cx + h, cy + 0.5f, cz - h};
    tri.v[2] = (b3d_point_t) {cx + h, cy - 0.5f, cz - h};
    b3d_triangle_lit(&tri, 0.7f, 0.0f, -0.7f, color);

    /* Diagonal 2 - back */
    tri.v[0] = (b3d_point_t) {cx + h, cy - 0.5f, cz - h};
    tri.v[1] = (b3d_point_t) {cx + h, cy + 0.5f, cz - h};
    tri.v[2] = (b3d_point_t) {cx - h, cy + 0.5f, cz + h};
    b3d_triangle_lit(&tri, -0.7f, 0.0f, 0.7f, color);

    tri.v[0] = (b3d_point_t) {cx + h, cy - 0.5f, cz - h};
    tri.v[1] = (b3d_point_t) {cx - h, cy + 0.5f, cz + h};
    tri.v[2] = (b3d_point_t) {cx - h, cy - 0.5f, cz + h};
    b3d_triangle_lit(&tri, -0.7f, 0.0f, 0.7f, color);
}

/* Linear interpolation between two colors */
static uint32_t lerp_color(uint32_t c1, uint32_t c2, float t)
{
    int r1 = COLOR_R(c1), g1 = COLOR_G(c1), b1 = COLOR_B(c1);
    int r2 = COLOR_R(c2), g2 = COLOR_G(c2), b2 = COLOR_B(c2);
    int r = r1 + (int) ((r2 - r1) * t);
    int g = g1 + (int) ((g2 - g1) * t);
    int b = b1 + (int) ((b2 - b1) * t);
    return COLOR_RGB(r, g, b);
}

/* Apply distance fog to the framebuffer using the depth buffer.
 * Uses integer thresholds to minimize per-pixel float operations.
 * Skips unwritten pixels (cleared depth = max value) automatically.
 */
static void apply_fog(uint32_t *pixels,
                      b3d_depth_t *depth_buf,
                      int width,
                      int height,
                      uint32_t fog_color)
{
    /* Pre-computed thresholds: 90% = 58982, 99% = 64880 (of 65536)
     * Cleared/unwritten pixels have max depth (> fog_end), so they're skipped.
     */
    const b3d_depth_t fog_start = 58982;
    const b3d_depth_t fog_end = 64880;
    const float inv_range = 0.25f / (float) (fog_end - fog_start);
    int count = width * height;

    for (int i = 0; i < count; ++i) {
        b3d_depth_t d = depth_buf[i];

        /* Skip: unwritten pixels, too far, or too close for fog */
        if (d >= fog_end || d < fog_start)
            continue;

        /* Linear fog in the 90%-99% range, scaled by 0.25 */
        float fog = (float) (d - fog_start) * inv_range;
        pixels[i] = lerp_color(pixels[i], fog_color, fog);
    }
}

/* Render the world from player's perspective using B3D */
static void render_world(uint32_t *pixels, b3d_depth_t *depth_buf, int w, int h)
{
    /* Day/night cycle - keep it bright like Minecraft daytime */
    float time_coef = (float) (world_time % 102944) / 16384.0f;
    time_coef = b3d_sinf(time_coef);
    time_coef /= b3d_sqrtf(time_coef * time_coef + (1.0f / 128.0f));
    time_coef = (time_coef + 1.0f) / 2.0f;
    if (time_coef < 0.9f)
        time_coef = 0.9f; /* Bright daytime like Minecraft */

    /* Sky color - Minecraft's classic bright blue sky */
    uint32_t sky_color =
        lerp_color(0x0A0A1A, 0x7BA4FF, time_coef); /* Night to Minecraft sky */

    /* Clear and fill sky */
    b3d_clear();
    for (int i = 0; i < w * h; ++i)
        pixels[i] = sky_color;

    /* Set up camera */
    b3d_set_camera(&(b3d_camera_t) {player_x, player_y, player_z, player_yaw,
                                    player_pitch, 0});

    /* Set lighting - warm sun from above, bright ambient like Minecraft */
    b3d_set_light_direction(0.3f, -0.9f, 0.3f);
    b3d_set_ambient(0.45f);

    /* Calculate render bounds */
    int min_x = (int) (player_x - RENDER_DISTANCE);
    int max_x = (int) (player_x + RENDER_DISTANCE);
    int min_z = (int) (player_z - RENDER_DISTANCE);
    int max_z = (int) (player_z + RENDER_DISTANCE);

    if (min_x < 0)
        min_x = 0;
    if (max_x >= MAP_WIDTH)
        max_x = MAP_WIDTH - 1;
    if (min_z < 0)
        min_z = 0;
    if (max_z >= MAP_DEPTH)
        max_z = MAP_DEPTH - 1;

    /* Render visible blocks - no matrix stack, direct world coordinates */
    b3d_reset();
    for (int x = min_x; x <= max_x; ++x) {
        for (int z = min_z; z <= max_z; ++z) {
            float dx = (float) x - player_x;
            float dz = (float) z - player_z;
            float dist_sq_xz = dx * dx + dz * dz;
            if (dist_sq_xz > RENDER_DISTANCE_SQ)
                continue;

            /* Use height_cache to bound Y-loop instead of scanning all 64 */
            int max_y = height_cache[x][z] + 1;
            if (max_y > MAP_HEIGHT)
                max_y = MAP_HEIGHT;

            for (int y = 0; y < max_y; ++y) {
                block_type_t block = world_map[x][z][y];
                if (block == BLOCK_AIR)
                    continue;

                /* Calculate 3D distance squared for LOD and culling */
                float dy = (float) y - player_y;
                float dist_sq = dist_sq_xz + dy * dy;

                /* 3D distance culling - skip blocks too far in Y */
                if (dist_sq > RENDER_DISTANCE_SQ)
                    continue;

                /* Tall grass uses cross billboard */
                if (block == BLOCK_TALL_GRASS) {
                    render_cross_billboard(x, y, z, block, time_coef);
                    continue;
                }

                /* Check if block is fully occluded (all 6 neighbors solid) */
                int visible_faces = 0;
                int neighbor_solid[6];
                static const int face_offsets[6][3] = {
                    {0, 0, -1}, {0, 0, 1}, {1, 0, 0},
                    {-1, 0, 0}, {0, 1, 0}, {0, -1, 0},
                };
                for (int face = 0; face < 6; ++face) {
                    int nbx = x + face_offsets[face][0];
                    int nby = y + face_offsets[face][1];
                    int nbz = z + face_offsets[face][2];
                    neighbor_solid[face] =
                        (nbx >= 0 && nbx < MAP_WIDTH && nbz >= 0 &&
                         nbz < MAP_DEPTH && nby >= 0 && nby < MAP_HEIGHT &&
                         !is_transparent(world_map[nbx][nbz][nby]));
                    if (!neighbor_solid[face])
                        visible_faces++;
                }

                /* Skip fully occluded blocks */
                if (visible_faces == 0)
                    continue;

                /* Render visible faces with distance-based tessellation */
                for (int face = 0; face < 6; ++face) {
                    /* Skip face if neighbor is solid (already computed) */
                    if (neighbor_solid[face])
                        continue;

                    /* Don't render water faces between water blocks */
                    if (block == BLOCK_WATER) {
                        int nbx = x + face_offsets[face][0];
                        int nby = y + face_offsets[face][1];
                        int nbz = z + face_offsets[face][2];
                        if (nbx >= 0 && nbx < MAP_WIDTH && nbz >= 0 &&
                            nbz < MAP_DEPTH && nby >= 0 && nby < MAP_HEIGHT &&
                            world_map[nbx][nbz][nby] == BLOCK_WATER)
                            continue;
                    }

                    render_face_tessellated(x, y, z, face, block, time_coef,
                                            dist_sq);
                }
            }
        }
    }

    /* Apply fog */
    apply_fog(pixels, depth_buf, w, h, sky_color);

    /* Draw crosshair */
    int cx = w / 2, cy = h / 2;
    for (int i = -4; i <= 4; ++i) {
        if (i != 0) {
            pixels[(cx + i) + cy * w] = 0xFFFFFF;
            pixels[cx + (cy + i) * w] = 0xFFFFFF;
        }
    }
}

int main(int argc, char **argv)
{
    int width = 640, height = 480;
    const char *snapshot = get_snapshot_path(argc, argv);

    /* Guard against allocation overflow */
    size_t pixel_count = (size_t) width * (size_t) height;
    if (pixel_count > SIZE_MAX / sizeof(uint32_t) ||
        pixel_count > SIZE_MAX / sizeof(b3d_depth_t)) {
        fprintf(stderr, "Framebuffer dimensions too large\n");
        return 1;
    }

    uint32_t *pixels = malloc(pixel_count * sizeof(pixels[0]));
    b3d_depth_t *depth = malloc(pixel_count * sizeof(depth[0]));

    if (!pixels || !depth) {
        fprintf(stderr, "Failed to allocate framebuffer\n");
        free(pixels);
        free(depth);
        return 1;
    }

    b3d_init(pixels, depth, width, height, 90);
    world_init();

    if (snapshot) {
        render_world(pixels, depth, width, height);
        write_png(snapshot, pixels, width, height);

        free(pixels);
        free(depth);
        return 0;
    }

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        free(pixels);
        free(depth);
        return 1;
    }

    SDL_Window *window =
        SDL_CreateWindow("Microcraft", SDL_WINDOWPOS_CENTERED,
                         SDL_WINDOWPOS_CENTERED, width, height, 0);
    if (!window) {
        fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
        free(pixels);
        free(depth);
        SDL_Quit();
        return 1;
    }

    SDL_Renderer *renderer =
        SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC);
    if (!renderer) {
        fprintf(stderr, "SDL_CreateRenderer failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        free(pixels);
        free(depth);
        SDL_Quit();
        return 1;
    }

    SDL_Texture *texture =
        SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888,
                          SDL_TEXTUREACCESS_STREAMING, width, height);
    if (!texture) {
        fprintf(stderr, "SDL_CreateTexture failed: %s\n", SDL_GetError());
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        free(pixels);
        free(depth);
        SDL_Quit();
        return 1;
    }

    SDL_SetRelativeMouseMode(SDL_TRUE);

    /* Movement state */
    int move_forward = 0, move_back = 0, move_left = 0, move_right = 0;
    int move_up = 0, move_down = 0;
    float move_speed = 8.0f;
    float velocity_x = 0.0f, velocity_z = 0.0f, velocity_y = 0.0f;

    Uint32 last_time = SDL_GetTicks();

    int quit = 0;
    while (!quit) {
        Uint32 current_time = SDL_GetTicks();
        float dt = (float) (current_time - last_time) / 1000.0f;
        last_time = current_time;

        world_time += 10;

        if (dt > 0.1f)
            dt = 0.1f;

        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
            } else if (event.type == SDL_KEYDOWN || event.type == SDL_KEYUP) {
                int pressed = event.type == SDL_KEYDOWN;
                switch (event.key.keysym.scancode) {
                case SDL_SCANCODE_W:
                case SDL_SCANCODE_UP:
                    move_forward = pressed;
                    break;
                case SDL_SCANCODE_S:
                case SDL_SCANCODE_DOWN:
                    move_back = pressed;
                    break;
                case SDL_SCANCODE_A:
                case SDL_SCANCODE_LEFT:
                    move_left = pressed;
                    break;
                case SDL_SCANCODE_D:
                case SDL_SCANCODE_RIGHT:
                    move_right = pressed;
                    break;
                case SDL_SCANCODE_SPACE:
                    move_up = pressed;
                    break;
                case SDL_SCANCODE_LSHIFT:
                case SDL_SCANCODE_RSHIFT:
                case SDL_SCANCODE_C:
                    move_down = pressed;
                    break;
                case SDL_SCANCODE_ESCAPE:
                    quit = 1;
                    break;
                default:
                    break;
                }
            } else if (event.type == SDL_MOUSEMOTION) {
                /*
                 * Mouse sensitivity is NOT scaled by delta time because:
                 * - SDL mouse events report accumulated physical movement
                 * - This is inherently frame-rate independent
                 * - Multiplying by dt would make sensitivity vary with FPS
                 */
                player_yaw -= (float) event.motion.xrel * MOUSE_SENSITIVITY;
                player_pitch += (float) event.motion.yrel * MOUSE_SENSITIVITY;

                if (player_pitch < -1.5f)
                    player_pitch = -1.5f;
                if (player_pitch > 1.5f)
                    player_pitch = 1.5f;
            }
        }

        /* Calculate movement */
        float sin_yaw = b3d_sinf(player_yaw);
        float cos_yaw = b3d_cosf(player_yaw);

        float target_vx = 0.0f, target_vz = 0.0f, target_vy = 0.0f;
        if (move_forward) {
            target_vx += sin_yaw * move_speed;
            target_vz += cos_yaw * move_speed;
        }
        if (move_back) {
            target_vx -= sin_yaw * move_speed;
            target_vz -= cos_yaw * move_speed;
        }
        if (move_left) {
            target_vx -= cos_yaw * move_speed;
            target_vz += sin_yaw * move_speed;
        }
        if (move_right) {
            target_vx += cos_yaw * move_speed;
            target_vz -= sin_yaw * move_speed;
        }
        if (move_up)
            target_vy += move_speed;
        if (move_down)
            target_vy -= move_speed;

        velocity_x += (target_vx - velocity_x) * 10.0f * dt;
        velocity_z += (target_vz - velocity_z) * 10.0f * dt;
        velocity_y += (target_vy - velocity_y) * 10.0f * dt;

        player_x += velocity_x * dt;
        player_z += velocity_z * dt;
        player_y += velocity_y * dt;

        /* Keep player in bounds */
        if (player_x < 1.0f)
            player_x = 1.0f;
        if (player_x > MAP_WIDTH - 2.0f)
            player_x = MAP_WIDTH - 2.0f;
        if (player_z < 1.0f)
            player_z = 1.0f;
        if (player_z > MAP_DEPTH - 2.0f)
            player_z = MAP_DEPTH - 2.0f;
        if (player_y < 1.0f)
            player_y = 1.0f;
        if (player_y > MAP_HEIGHT - 2.0f)
            player_y = MAP_HEIGHT - 2.0f;

        /* Render */
        render_world(pixels, depth, width, height);

        SDL_RenderClear(renderer);
        SDL_UpdateTexture(texture, NULL, pixels,
                          width * (int) sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }

    free(pixels);
    free(depth);
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
