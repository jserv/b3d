/*
 * 3D Earth Globe Demo
 *
 * Renders a rotating 3D globe with country boundaries from Natural Earth data.
 *
 * Uses data from Natural Earth (https://www.naturalearthdata.com/),
 * specifically their 1:10m countries dataset.
 * Compact binary format with delta encoding:
 *   - uint32: polygon count
 *   - Per polygon:
 *     - uint16: point count
 *     - int16, int16: first point (quantized lon, lat)
 *     - Subsequent: varint-encoded deltas (ZigZag + LEB128)
 *
 * Generate data with: python3 scripts/gen-globe-data.py
 *
 * Supports headless snapshots with '--snapshot=PATH' or B3D_SNAPSHOT.
 */

#include <SDL.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "b3d-math.h"
#include "b3d.h"
#include "globe-data.h"
#include "utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Quantization scale for coordinate encoding (int16 range) */
#define QUANTIZATION_SCALE 32767.0

/* Display configuration */
#define GLOBE_WIDTH 800
#define GLOBE_HEIGHT 600

/* Rendering configuration */
#define GLOBE_RADIUS 1.5f
#define ROTATION_SPEED 0.4f /* radians per second */
#define LINE_THICKNESS 0.012f
#define CAMERA_DISTANCE 4.5f

/* Colors - classic globe style */
#define COLOR_OCEAN 0x1e5799     /* classic blue ocean */
#define COLOR_LAND 0x5a8f5a      /* muted green land */
#define COLOR_COASTLINE 0xd4c4a8 /* tan/beige coastlines */
#define COLOR_CITY 0xcc3333      /* red city markers */

/* City marker size */
#define CITY_MARKER_SIZE 0.04f

/* Vector operation macros (self-contained, no internal header deps) */
#define VEC3_LEN(x, y, z) b3d_sqrtf((x) * (x) + (y) * (y) + (z) * (z))

#define VEC3_NORMALIZE(x, y, z)                    \
    do {                                           \
        float _len = VEC3_LEN(x, y, z);            \
        if (_len > 0.0001f)                        \
            (x) /= _len, (y) /= _len, (z) /= _len; \
    } while (0)

#define VEC3_CROSS(rx, ry, rz, ax, ay, az, bx, by, bz) \
    do {                                               \
        (rx) = (ay) * (bz) - (az) * (by);              \
        (ry) = (az) * (bx) - (ax) * (bz);              \
        (rz) = (ax) * (by) - (ay) * (bx);              \
    } while (0)

#define VEC3_AVG3(rx, ry, rz, a, b, c)            \
    do {                                          \
        (rx) = ((a)[0] + (b)[0] + (c)[0]) / 3.0f; \
        (ry) = ((a)[1] + (b)[1] + (c)[1]) / 3.0f; \
        (rz) = ((a)[2] + (b)[2] + (c)[2]) / 3.0f; \
    } while (0)

/* Globe rendering parameters */
#define SPHERE_SEGMENTS 48
#define SPHERE_RINGS 32
#define MAX_SPHERE_TRIS (SPHERE_SEGMENTS * SPHERE_RINGS * 2)

/* City marker with geographic coordinates (radians) */
typedef struct {
    const char *name;
    double lat, lon;
} city_t;

/* Convert degrees to radians */
#define DEG_TO_RAD(deg) ((deg) * M_PI / 180.0)

/* 10 major world cities (Taipei first as required) */
static const city_t cities[] = {
    {"Taipei", DEG_TO_RAD(25.0330), DEG_TO_RAD(121.5654)},
    {"Tokyo", DEG_TO_RAD(35.6762), DEG_TO_RAD(139.6503)},
    {"New York", DEG_TO_RAD(40.7128), DEG_TO_RAD(-74.0060)},
    {"London", DEG_TO_RAD(51.5074), DEG_TO_RAD(-0.1278)},
    {"Paris", DEG_TO_RAD(48.8566), DEG_TO_RAD(2.3522)},
    {"Sydney", DEG_TO_RAD(-33.8688), DEG_TO_RAD(151.2093)},
    {"Cairo", DEG_TO_RAD(30.0444), DEG_TO_RAD(31.2357)},
    {"Santiago", DEG_TO_RAD(-33.4489), DEG_TO_RAD(-70.6693)},
    {"Kyiv", DEG_TO_RAD(50.4501), DEG_TO_RAD(30.5234)},
    {"Jakarta", DEG_TO_RAD(-6.2088), DEG_TO_RAD(106.8456)},
};
#define NUM_CITIES (sizeof(cities) / sizeof(cities[0]))

/* Static buffers for sphere mesh */
static b3d_tri_t sphere_tris[MAX_SPHERE_TRIS];
static float sphere_normals[MAX_SPHERE_TRIS][3];
static int sphere_tri_count = 0;

/* Cached country boundary data (pre-computed at startup) */
typedef struct {
    float (*points)[3]; /* array of 3D points */
    int point_count;
} cached_poly_t;

static cached_poly_t *country_polys = NULL;
static int country_poly_count = 0;

/* Forward declarations */
static void cleanup_countries(void);

/* Convert geographic coordinates to 3D Cartesian (unit sphere).
 * @x, @y, @z: output Cartesian coordinates (unit sphere)
 * @lon: longitude in radians
 * @lat: latitude in radians
 *
 * Coordinate system (Y-up to match sphere mesh):
 *   X = cos(lat) * cos(lon)  - toward prime meridian at equator
 *   Y = sin(lat)             - toward north pole (Y-up)
 *   Z = cos(lat) * sin(lon)  - toward 90E at equator
 */
static void geo_to_xyz(float *x, float *y, float *z, float lon, float lat)
{
    float sin_lat, cos_lat, sin_lon, cos_lon;
    b3d_sincosf(lat, &sin_lat, &cos_lat);
    b3d_sincosf(lon, &sin_lon, &cos_lon);

    *x = cos_lat * cos_lon;
    *y = sin_lat;
    *z = cos_lat * sin_lon;
}

/* Generate UV sphere mesh for ocean background.
 * Normals point outward for correct lighting from an external camera.
 */
static void generate_sphere(void)
{
    sphere_tri_count = 0;

    for (int ring = 0; ring < SPHERE_RINGS; ring++) {
        float phi0 = (float) M_PI * ring / SPHERE_RINGS;
        float phi1 = (float) M_PI * (ring + 1) / SPHERE_RINGS;

        float sin_phi0, cos_phi0, sin_phi1, cos_phi1;
        b3d_sincosf(phi0, &sin_phi0, &cos_phi0);
        b3d_sincosf(phi1, &sin_phi1, &cos_phi1);

        for (int seg = 0; seg < SPHERE_SEGMENTS; seg++) {
            float theta0 = 2.0f * (float) M_PI * seg / SPHERE_SEGMENTS;
            float theta1 = 2.0f * (float) M_PI * (seg + 1) / SPHERE_SEGMENTS;

            float sin_theta0, cos_theta0, sin_theta1, cos_theta1;
            b3d_sincosf(theta0, &sin_theta0, &cos_theta0);
            b3d_sincosf(theta1, &sin_theta1, &cos_theta1);

            /* Four corners of the quad (Y is up, X is right, Z is forward) */
            float v00[3] = {sin_phi0 * cos_theta0, cos_phi0,
                            sin_phi0 * sin_theta0};
            float v01[3] = {sin_phi0 * cos_theta1, cos_phi0,
                            sin_phi0 * sin_theta1};
            float v10[3] = {sin_phi1 * cos_theta0, cos_phi1,
                            sin_phi1 * sin_theta0};
            float v11[3] = {sin_phi1 * cos_theta1, cos_phi1,
                            sin_phi1 * sin_theta1};

            if (sphere_tri_count + 2 > MAX_SPHERE_TRIS)
                break;

            /* First triangle (upper-left of quad) */
            b3d_tri_t *tri1 = &sphere_tris[sphere_tri_count];
            tri1->v[0] = (b3d_point_t) {v00[0], v00[1], v00[2]};
            tri1->v[1] = (b3d_point_t) {v01[0], v01[1], v01[2]};
            tri1->v[2] = (b3d_point_t) {v11[0], v11[1], v11[2]};

            /* Normal for face (average of vertices for smooth sphere) */
            float nx, ny, nz;
            VEC3_AVG3(nx, ny, nz, v00, v01, v11);
            VEC3_NORMALIZE(nx, ny, nz);
            sphere_normals[sphere_tri_count][0] = nx;
            sphere_normals[sphere_tri_count][1] = ny;
            sphere_normals[sphere_tri_count][2] = nz;
            sphere_tri_count++;

            /* Second triangle (lower-right of quad) */
            b3d_tri_t *tri2 = &sphere_tris[sphere_tri_count];
            tri2->v[0] = (b3d_point_t) {v00[0], v00[1], v00[2]};
            tri2->v[1] = (b3d_point_t) {v11[0], v11[1], v11[2]};
            tri2->v[2] = (b3d_point_t) {v10[0], v10[1], v10[2]};

            VEC3_AVG3(nx, ny, nz, v00, v11, v10);
            VEC3_NORMALIZE(nx, ny, nz);
            sphere_normals[sphere_tri_count][0] = nx;
            sphere_normals[sphere_tri_count][1] = ny;
            sphere_normals[sphere_tri_count][2] = nz;
            sphere_tri_count++;
        }
    }
}

/* Render a line segment on the globe surface as a thin quad.
 * @p0, @p1: endpoints in 3D (unit sphere coordinates)
 * @color: line color
 */
static void render_globe_line(const float p0[3],
                              const float p1[3],
                              uint32_t color)
{
    /* Calculate perpendicular vector for line thickness */
    float dx = p1[0] - p0[0], dy = p1[1] - p0[1], dz = p1[2] - p0[2];

    /* Cross product with radial direction to get perpendicular */
    float cx, cy, cz;
    VEC3_CROSS(cx, cy, cz, p0[0], p0[1], p0[2], dx, dy, dz);

    float len = VEC3_LEN(cx, cy, cz);
    if (len < 0.0001f)
        return;

    /* Normalize and scale to line thickness */
    float scale = LINE_THICKNESS / len;
    cx *= scale, cy *= scale, cz *= scale;

    /* Create quad vertices (two triangles) */
    b3d_tri_t tri1 = {
        .v =
            {
                {p0[0] - cx, p0[1] - cy, p0[2] - cz},
                {p0[0] + cx, p0[1] + cy, p0[2] + cz},
                {p1[0] + cx, p1[1] + cy, p1[2] + cz},
            },
    };

    b3d_tri_t tri2 = {
        .v =
            {
                {p0[0] - cx, p0[1] - cy, p0[2] - cz},
                {p1[0] + cx, p1[1] + cy, p1[2] + cz},
                {p1[0] - cx, p1[1] - cy, p1[2] - cz},
            },
    };

    b3d_triangle(&tri1, color);
    b3d_triangle(&tri2, color);
}

/* Decode LEB128 variable-length unsigned integer with bounds checking.
 * @ptr: pointer to data pointer (advanced on return)
 * @end: pointer past the end of valid data
 * @out: output for decoded value
 * Returns true on success, false if bounds exceeded or invalid data.
 */
static bool read_varint(const uint8_t **ptr, const uint8_t *end, uint32_t *out)
{
    uint32_t result = 0;
    uint32_t shift = 0;

    while (*ptr < end) {
        uint8_t byte = *(*ptr)++;
        result |= (uint32_t) (byte & 0x7F) << shift;
        if ((byte & 0x80) == 0) {
            *out = result;
            return true;
        }
        shift += 7;
        if (shift > 28) /* prevent overflow: max 4 bytes for uint32 */
            return false;
    }
    return false; /* truncated data */
}

/* Decode ZigZag-encoded integer to signed.
 * @n: unsigned ZigZag-encoded value
 * Returns signed integer.
 */
static int32_t zigzag_decode(uint32_t n)
{
    return (int32_t) (n >> 1) ^ -((int32_t) (n & 1));
}

/* Initialize country boundary data (pre-compute geometry at startup).
 * Decodes delta-encoded varint data once, caching 3D points.
 * Returns true on success, false on allocation or decode failure.
 */
static bool init_countries(void)
{
    const uint8_t *ptr = globe_data;
    const uint8_t *end = globe_data + sizeof(globe_data);

    if (ptr + sizeof(uint32_t) > end)
        return false;

    uint32_t polygon_count;
    memcpy(&polygon_count, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    /* Sanity check: limit polygon count to prevent OOM from malformed data.
     * Natural Earth 10m has ~2000 polygons; 10000 is a generous upper bound.
     */
    if (polygon_count > 10000)
        return false;

    country_polys = calloc(polygon_count, sizeof(cached_poly_t));
    if (!country_polys)
        return false;

    float scale = GLOBE_RADIUS * 1.01f;
    int valid_count = 0;

    for (uint32_t poly = 0; poly < polygon_count; poly++) {
        if (ptr + sizeof(uint16_t) > end)
            break;

        uint16_t point_count;
        memcpy(&point_count, ptr, sizeof(uint16_t));
        ptr += sizeof(uint16_t);

        if (point_count == 0)
            continue;

        /* Sanity check: limit per-polygon points to prevent malformed data OOM.
         * Largest real polygon in Natural Earth 10m has ~13000 points (Russia).
         */
        if (point_count > 15000)
            continue; /* skip oversize, don't break entire loop */

        float (*points)[3] = malloc(point_count * sizeof(float[3]));
        if (!points)
            goto fail;

        int16_t last_lon_q = 0, last_lat_q = 0;
        bool decode_ok = true;

        for (uint16_t i = 0; i < point_count && decode_ok; i++) {
            if (i == 0) {
                if (ptr + sizeof(int16_t) * 2 > end) {
                    decode_ok = false;
                    break;
                }
                memcpy(&last_lon_q, ptr, sizeof(int16_t));
                ptr += sizeof(int16_t);
                memcpy(&last_lat_q, ptr, sizeof(int16_t));
                ptr += sizeof(int16_t);
            } else {
                uint32_t d_lon_u, d_lat_u;
                if (!read_varint(&ptr, end, &d_lon_u) ||
                    !read_varint(&ptr, end, &d_lat_u)) {
                    decode_ok = false;
                    break;
                }
                last_lon_q += (int16_t) zigzag_decode(d_lon_u);
                last_lat_q += (int16_t) zigzag_decode(d_lat_u);
            }

            double lon = last_lon_q * (M_PI / QUANTIZATION_SCALE);
            double lat = last_lat_q * ((M_PI / 2.0) / QUANTIZATION_SCALE);

            float x, y, z;
            geo_to_xyz(&x, &y, &z, lon, lat);

            points[i][0] = x * scale;
            points[i][1] = y * scale;
            points[i][2] = z * scale;
        }

        if (!decode_ok) {
            free(points);
            continue; /* skip malformed polygon */
        }

        country_polys[valid_count].points = points;
        country_polys[valid_count].point_count = point_count;
        valid_count++;
    }

    country_poly_count = valid_count;
    return true;

fail:
    /* Set count before cleanup so it frees already-allocated polygons */
    country_poly_count = valid_count;
    cleanup_countries();
    return false;
}

/* Free country boundary data allocated by init_countries. */
static void cleanup_countries(void)
{
    if (country_polys) {
        for (int i = 0; i < country_poly_count; i++)
            free(country_polys[i].points);
        free(country_polys);
        country_polys = NULL;
        country_poly_count = 0;
    }
}

/* Render pre-computed country boundaries as lines. */
static void render_country_boundaries(void)
{
    for (int i = 0; i < country_poly_count; i++) {
        cached_poly_t *poly = &country_polys[i];
        for (int j = 0; j < poly->point_count - 1; j++)
            render_globe_line(poly->points[j], poly->points[j + 1],
                              COLOR_COASTLINE);
    }
}

/* Number of segments for circular city markers */
#define CITY_MARKER_SEGMENTS 8

/* Render a city marker as a circular dot on the globe surface.
 * @city: city data with coordinates
 * @color: marker color
 */
static void render_city_marker(const city_t *city, uint32_t color)
{
    float nx, ny, nz;
    geo_to_xyz(&nx, &ny, &nz, city->lon, city->lat);

    /* Scale to globe radius (above grid to avoid z-fighting) */
    float scale = GLOBE_RADIUS * 1.05f;
    float cx = nx * scale;
    float cy = ny * scale;
    float cz = nz * scale;

    /* Calculate tangent vectors using robust reference axis.
     * Use X-axis when near poles, Y-axis otherwise.
     */
    float tx, ty, tz;
    if (b3d_fabsf(ny) > 0.9f)
        VEC3_CROSS(tx, ty, tz, nx, ny, nz, 1, 0, 0); /* cross with X-axis */
    else
        VEC3_CROSS(tx, ty, tz, nx, ny, nz, 0, 1, 0); /* cross with Y-axis */
    VEC3_NORMALIZE(tx, ty, tz);

    /* Second tangent: cross product of normal and first tangent */
    float ux, uy, uz;
    VEC3_CROSS(ux, uy, uz, nx, ny, nz, tx, ty, tz);
    VEC3_NORMALIZE(ux, uy, uz);

    /* Generate circular marker using triangle fan */
    float size = CITY_MARKER_SIZE;
    float prev_vx = cx + tx * size;
    float prev_vy = cy + ty * size;
    float prev_vz = cz + tz * size;

    for (int i = 1; i <= CITY_MARKER_SEGMENTS; i++) {
        float angle = 2.0f * (float) M_PI * i / CITY_MARKER_SEGMENTS;
        float cos_a, sin_a;
        b3d_sincosf(angle, &sin_a, &cos_a);

        float vx = cx + (tx * cos_a + ux * sin_a) * size;
        float vy = cy + (ty * cos_a + uy * sin_a) * size;
        float vz = cz + (tz * cos_a + uz * sin_a) * size;

        b3d_tri_t tri = {
            .v = {{cx, cy, cz}, {prev_vx, prev_vy, prev_vz}, {vx, vy, vz}},
        };
        b3d_triangle(&tri, color);

        prev_vx = vx, prev_vy = vy, prev_vz = vz;
    }
}

/* Render all city markers on the globe */
static void render_cities(void)
{
    for (size_t i = 0; i < NUM_CITIES; i++)
        render_city_marker(&cities[i], COLOR_CITY);
}

/* Render the complete globe scene.
 * @pixels: pixel buffer
 * @depth: depth buffer
 * @width: buffer width
 * @height: buffer height
 * @rotation: current globe rotation angle in radians
 */
static void render_scene(uint32_t *pixels,
                         b3d_depth_t *depth,
                         int width,
                         int height,
                         float rotation)
{
    b3d_init(pixels, depth, width, height, 45);
    b3d_set_camera(&(b3d_camera_t) {0, 0, -CAMERA_DISTANCE, 0, 0, 0});

    b3d_clear();

    /* Fill background with deep space color */
    for (int i = 0; i < width * height; i++)
        pixels[i] = 0x000510;

    /* Set up lighting (sun from upper-right-front) */
    b3d_set_light_direction(0.5f, 0.3f, -0.8f);
    b3d_set_ambient(0.15f);

    /* Render ocean sphere */
    b3d_reset();
    b3d_rotate_y(rotation);
    b3d_scale(GLOBE_RADIUS, GLOBE_RADIUS, GLOBE_RADIUS);

    for (int i = 0; i < sphere_tri_count; i++) {
        b3d_triangle_lit(&sphere_tris[i], sphere_normals[i][0],
                         sphere_normals[i][1], sphere_normals[i][2],
                         COLOR_OCEAN);
    }

    /* Render country boundaries */
    b3d_reset();
    b3d_rotate_y(rotation);
    render_country_boundaries();

    /* Render city markers */
    render_cities();
}

int main(int argc, char **argv)
{
    int width = GLOBE_WIDTH;
    int height = GLOBE_HEIGHT;
    const char *snapshot = get_snapshot_path(argc, argv);

    /* Allocate buffers */
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    if (!pixels) {
        fprintf(stderr, "Failed to allocate pixel buffer\n");
        return 1;
    }

    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));
    if (!depth) {
        fprintf(stderr, "Failed to allocate depth buffer\n");
        free(pixels);
        return 1;
    }

    /* Generate sphere mesh and pre-compute country boundaries */
    generate_sphere();
    if (!init_countries()) {
        fprintf(stderr, "Failed to initialize country data\n");
        free(pixels);
        free(depth);
        return 1;
    }

    if (snapshot) {
        /* Render single frame for snapshot (show Asia with Taipei) */
        render_scene(pixels, depth, width, height, 2.0f);
        write_png(snapshot, pixels, width, height);
        cleanup_countries();
        free(pixels);
        free(depth);
        return 0;
    }

    /* SDL2 setup with error checking */
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        cleanup_countries();
        free(pixels);
        free(depth);
        return 1;
    }

    SDL_Window *window =
        SDL_CreateWindow("B3D Globe", SDL_WINDOWPOS_CENTERED,
                         SDL_WINDOWPOS_CENTERED, width, height, 0);
    if (!window) {
        fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
        SDL_Quit();
        cleanup_countries();
        free(pixels);
        free(depth);
        return 1;
    }

    SDL_Renderer *renderer =
        SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC);
    if (!renderer) {
        fprintf(stderr, "SDL_CreateRenderer failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        cleanup_countries();
        free(pixels);
        free(depth);
        return 1;
    }

    SDL_Texture *texture =
        SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888,
                          SDL_TEXTUREACCESS_STREAMING, width, height);
    if (!texture) {
        fprintf(stderr, "SDL_CreateTexture failed: %s\n", SDL_GetError());
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
        cleanup_countries();
        free(pixels);
        free(depth);
        return 1;
    }

    uint32_t last_time = SDL_GetTicks();
    float rotation = 0.0f;

    int quit = 0;
    while (!quit) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT ||
                (event.type == SDL_KEYDOWN &&
                 event.key.keysym.scancode == SDL_SCANCODE_ESCAPE)) {
                quit = 1;
                break;
            }
        }
        if (quit)
            break;

        /* Calculate delta time */
        uint32_t current_time = SDL_GetTicks();
        float dt = (current_time - last_time) / 1000.0f;
        last_time = current_time;

        /* Cap delta time */
        if (dt > 0.1f)
            dt = 0.1f;

        /* Update rotation */
        rotation += ROTATION_SPEED * dt;

        /* Render */
        render_scene(pixels, depth, width, height, rotation);

        /* Display */
        SDL_RenderClear(renderer);
        SDL_UpdateTexture(texture, NULL, pixels, width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }

    /* Cleanup */
    cleanup_countries();
    free(pixels);
    free(depth);
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
