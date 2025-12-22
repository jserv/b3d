/*
 * B3D Performance Benchmark Tests
 * Measures rendering throughput and critical path performance
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/b3d.h"

/* ANSI color codes for terminal output */
#define ANSI_GREEN "\033[32m"
#define ANSI_YELLOW "\033[33m"
#define ANSI_BOLD "\033[1m"
#define ANSI_RESET "\033[0m"

/* Benchmark configuration */
#define WARMUP_ITERATIONS 100
#define BENCHMARK_DURATION_MS 1000

/* Get current time in milliseconds */
static double get_time_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1000000.0;
}

/* Benchmark result structure */
typedef struct {
    char *name; /* Dynamically allocated name (caller must free) */
    double ops_per_sec;
    double avg_time_us;
    size_t iterations;
} bench_result_t;

static void print_result(bench_result_t *result)
{
    printf("  %-36s %12.0f ops/s  %8.3f us/op\n", result->name,
           result->ops_per_sec, result->avg_time_us);
    free(result->name);
    result->name = NULL;
}

/* Helper to set up benchmark environment */
static bool setup_benchmark(int width,
                            int height,
                            uint32_t **out_pixels,
                            b3d_depth_t **out_depth)
{
    if (width <= 0 || height <= 0)
        return false;

    size_t pixel_bytes = b3d_buffer_size(width, height, sizeof(uint32_t));
    size_t depth_bytes = b3d_buffer_size(width, height, sizeof(b3d_depth_t));
    if (pixel_bytes == 0 || depth_bytes == 0)
        return false;

    uint32_t *pixels = malloc(pixel_bytes);
    b3d_depth_t *depth = malloc(depth_bytes);
    if (!pixels || !depth) {
        free(pixels);
        free(depth);
        return false;
    }

    if (!b3d_init(pixels, depth, width, height, 65.0f)) {
        free(pixels);
        free(depth);
        return false;
    }

    *out_pixels = pixels;
    *out_depth = depth;
    return true;
}

/* Cube vertices for benchmarking */
static void render_cube(float angle)
{
    b3d_reset();
    b3d_rotate_y(angle);
    b3d_rotate_x(angle * 0.7f);

    /* Front face */
    b3d_triangle(-0.5f, -0.5f, -0.5f, -0.5f, 0.5f, -0.5f, 0.5f, 0.5f, -0.5f,
                 0xfcd0a1);
    b3d_triangle(-0.5f, -0.5f, -0.5f, 0.5f, 0.5f, -0.5f, 0.5f, -0.5f, -0.5f,
                 0xb1b695);
    /* Right face */
    b3d_triangle(0.5f, -0.5f, -0.5f, 0.5f, 0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
                 0x53917e);
    b3d_triangle(0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f, 0.5f, -0.5f, 0.5f,
                 0x63535b);
    /* Back face */
    b3d_triangle(0.5f, -0.5f, 0.5f, 0.5f, 0.5f, 0.5f, -0.5f, 0.5f, 0.5f,
                 0x6d1a36);
    b3d_triangle(0.5f, -0.5f, 0.5f, -0.5f, 0.5f, 0.5f, -0.5f, -0.5f, 0.5f,
                 0xd4e09b);
    /* Left face */
    b3d_triangle(-0.5f, -0.5f, 0.5f, -0.5f, 0.5f, 0.5f, -0.5f, 0.5f, -0.5f,
                 0xf6f4d2);
    b3d_triangle(-0.5f, -0.5f, 0.5f, -0.5f, 0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                 0xcbdfbd);
    /* Top face */
    b3d_triangle(-0.5f, 0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                 0xf19c79);
    b3d_triangle(-0.5f, 0.5f, -0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, -0.5f,
                 0xa44a3f);
    /* Bottom face */
    b3d_triangle(0.5f, -0.5f, 0.5f, -0.5f, -0.5f, 0.5f, -0.5f, -0.5f, -0.5f,
                 0x5465ff);
    b3d_triangle(0.5f, -0.5f, 0.5f, -0.5f, -0.5f, -0.5f, 0.5f, -0.5f, -0.5f,
                 0x788bff);
}

/*
 * Benchmark: Triangle rendering throughput
 */
static bench_result_t bench_triangles(int width, int height)
{
    bench_result_t result = {
        .name = strdup("Triangle rendering"),
        .ops_per_sec = 0.0,
        .avg_time_us = 0.0,
        .iterations = 0,
    };

    uint32_t *pixels = NULL;
    b3d_depth_t *depth = NULL;
    if (!setup_benchmark(width, height, &pixels, &depth))
        return result;

    b3d_set_camera(0.0f, 0.0f, -3.0f, 0.0f, 0.0f, 0.0f);

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        b3d_clear();
        b3d_triangle(-0.5f, -0.5f, 0.5f, 0.5f, -0.5f, 0.5f, 0.0f, 0.5f, 0.5f,
                     0xffffff);
    }

    /* Benchmark - check time before each batch for precision */
    size_t iterations = 0;
    double start = get_time_ms();

    while (get_time_ms() - start < BENCHMARK_DURATION_MS) {
        b3d_clear();
        for (int i = 0; i < 100; i++) {
            b3d_triangle(-0.5f + (float) i * 0.001f, -0.5f, 0.5f,
                         0.5f + (float) i * 0.001f, -0.5f, 0.5f, 0.0f, 0.5f,
                         0.5f, 0xffffff);
        }
        iterations += 100;
    }

    double elapsed = get_time_ms() - start;
    free(pixels);
    free(depth);

    if (iterations == 0 || elapsed <= 0.0)
        return result;

    result.ops_per_sec = (double) iterations / (elapsed / 1000.0);
    result.avg_time_us = (elapsed * 1000.0) / (double) iterations;
    result.iterations = iterations;
    return result;
}

/*
 * Benchmark: Cube rendering (12 triangles per frame)
 */
static bench_result_t bench_cubes(int width, int height)
{
    bench_result_t result = {
        .name = strdup("Cube rendering (12 tris)"),
        .ops_per_sec = 0.0,
        .avg_time_us = 0.0,
        .iterations = 0,
    };

    uint32_t *pixels = NULL;
    b3d_depth_t *depth = NULL;
    if (!setup_benchmark(width, height, &pixels, &depth))
        return result;

    b3d_set_camera(0.0f, 0.0f, -3.0f, 0.0f, 0.0f, 0.0f);

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        b3d_clear();
        render_cube((float) i * 0.1f);
    }

    /* Benchmark - check time before each iteration for precision */
    size_t iterations = 0;
    double start = get_time_ms();

    while (get_time_ms() - start < BENCHMARK_DURATION_MS) {
        b3d_clear();
        render_cube((float) iterations * 0.1f);
        iterations++;
    }

    double elapsed = get_time_ms() - start;
    free(pixels);
    free(depth);

    if (iterations == 0 || elapsed <= 0.0)
        return result;

    result.ops_per_sec = (double) iterations / (elapsed / 1000.0);
    result.avg_time_us = (elapsed * 1000.0) / (double) iterations;
    result.iterations = iterations;
    return result;
}

/*
 * Benchmark: b3d_clear() performance
 */
static bench_result_t bench_clear(int width, int height)
{
    bench_result_t result = {
        .name = strdup("Buffer clear"),
        .ops_per_sec = 0.0,
        .avg_time_us = 0.0,
        .iterations = 0,
    };

    uint32_t *pixels = NULL;
    b3d_depth_t *depth = NULL;
    if (!setup_benchmark(width, height, &pixels, &depth))
        return result;

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++)
        b3d_clear();

    /* Benchmark - check time before each iteration for precision */
    size_t iterations = 0;
    double start = get_time_ms();

    while (get_time_ms() - start < BENCHMARK_DURATION_MS) {
        b3d_clear();
        iterations++;
    }

    double elapsed = get_time_ms() - start;
    free(pixels);
    free(depth);

    if (iterations == 0 || elapsed <= 0.0)
        return result;

    result.ops_per_sec = (double) iterations / (elapsed / 1000.0);
    result.avg_time_us = (elapsed * 1000.0) / (double) iterations;
    result.iterations = iterations;
    return result;
}

/*
 * Benchmark: Matrix operations
 */
static bench_result_t bench_matrix_ops(void)
{
    bench_result_t result = {
        .name = strdup("Matrix chain (6 ops)"),
        .ops_per_sec = 0.0,
        .avg_time_us = 0.0,
        .iterations = 0,
    };

    const int width = 64, height = 64;
    uint32_t *pixels = NULL;
    b3d_depth_t *depth = NULL;
    if (!setup_benchmark(width, height, &pixels, &depth))
        return result;

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        b3d_reset();
        b3d_translate(1.0f, 2.0f, 3.0f);
        b3d_rotate_x(0.5f), b3d_rotate_y(0.5f), b3d_rotate_z(0.5f);
        b3d_scale(2.0f, 2.0f, 2.0f);
    }

    /* Benchmark - check time before each iteration for precision */
    size_t iterations = 0;
    double start = get_time_ms();

    while (get_time_ms() - start < BENCHMARK_DURATION_MS) {
        b3d_reset();
        b3d_translate(1.0f, 2.0f, 3.0f);
        b3d_rotate_x(0.5f), b3d_rotate_y(0.5f), b3d_rotate_z(0.5f);
        b3d_scale(2.0f, 2.0f, 2.0f);
        iterations++;
    }

    double elapsed = get_time_ms() - start;
    free(pixels);
    free(depth);

    if (iterations == 0 || elapsed <= 0.0)
        return result;

    result.ops_per_sec = (double) iterations / (elapsed / 1000.0);
    result.avg_time_us = (elapsed * 1000.0) / (double) iterations;
    result.iterations = iterations;
    return result;
}

/*
 * Benchmark: Matrix stack push/pop
 */
static bench_result_t bench_matrix_stack(void)
{
    bench_result_t result = {
        .name = strdup("Matrix push/pop cycle"),
        .ops_per_sec = 0.0,
        .avg_time_us = 0.0,
        .iterations = 0,
    };

    const int width = 64, height = 64;
    uint32_t *pixels = NULL;
    b3d_depth_t *depth = NULL;
    if (!setup_benchmark(width, height, &pixels, &depth))
        return result;

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        b3d_push_matrix();
        b3d_translate(1.0f, 0.0f, 0.0f);
        b3d_pop_matrix();
    }

    /* Benchmark - check time before each iteration for precision */
    size_t iterations = 0;
    double start = get_time_ms();

    while (get_time_ms() - start < BENCHMARK_DURATION_MS) {
        b3d_push_matrix();
        b3d_translate(1.0f, 0.0f, 0.0f);
        b3d_pop_matrix();
        iterations++;
    }

    double elapsed = get_time_ms() - start;
    free(pixels);
    free(depth);

    if (iterations == 0 || elapsed <= 0.0)
        return result;

    result.ops_per_sec = (double) iterations / (elapsed / 1000.0);
    result.avg_time_us = (elapsed * 1000.0) / (double) iterations;
    result.iterations = iterations;
    return result;
}

/*
 * Benchmark: Screen projection
 */
static bench_result_t bench_to_screen(void)
{
    bench_result_t result = {
        .name = strdup("Screen projection"),
        .ops_per_sec = 0.0,
        .avg_time_us = 0.0,
        .iterations = 0,
    };

    const int width = 640, height = 480;
    uint32_t *pixels = NULL;
    b3d_depth_t *depth = NULL;
    if (!setup_benchmark(width, height, &pixels, &depth))
        return result;

    b3d_set_camera(0.0f, 0.0f, -3.0f, 0.0f, 0.0f, 0.0f);

    int sx, sy;

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++)
        b3d_to_screen(0.0f, 0.0f, 1.0f, &sx, &sy);

    /* Benchmark - check time before each batch for precision */
    size_t iterations = 0;
    double start = get_time_ms();

    while (get_time_ms() - start < BENCHMARK_DURATION_MS) {
        for (int i = 0; i < 100; i++)
            b3d_to_screen((float) i * 0.01f, 0.0f, 1.0f, &sx, &sy);
        iterations += 100;
    }

    double elapsed = get_time_ms() - start;
    free(pixels);
    free(depth);

    if (iterations == 0 || elapsed <= 0.0)
        return result;

    result.ops_per_sec = (double) iterations / (elapsed / 1000.0);
    result.avg_time_us = (elapsed * 1000.0) / (double) iterations;
    result.iterations = iterations;
    return result;
}

/*
 * Benchmark: Full frame (clear + render cube)
 */
static bench_result_t bench_full_frame(int width, int height)
{
    bench_result_t result = {
        .name = NULL,
        .ops_per_sec = 0.0,
        .avg_time_us = 0.0,
        .iterations = 0,
    };

    uint32_t *pixels = NULL;
    b3d_depth_t *depth = NULL;
    if (!setup_benchmark(width, height, &pixels, &depth)) {
        result.name = strdup("Full frame (setup failed)");
        return result;
    }

    b3d_set_camera(0.0f, 0.0f, -3.0f, 0.0f, 0.0f, 0.0f);

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        b3d_clear();
        render_cube((float) i * 0.1f);
    }

    /* Benchmark - check time before each iteration for precision */
    size_t iterations = 0;
    double start = get_time_ms();

    while (get_time_ms() - start < BENCHMARK_DURATION_MS) {
        b3d_clear();
        render_cube((float) iterations * 0.1f);
        iterations++;
    }

    double elapsed = get_time_ms() - start;
    free(pixels);
    free(depth);

    if (iterations == 0 || elapsed <= 0.0) {
        result.name = strdup("Full frame (no iterations)");
        return result;
    }

    double fps = (double) iterations / (elapsed / 1000.0);

    /* Dynamically allocate the name with resolution and FPS info */
    char name_buf[64];
    snprintf(name_buf, sizeof(name_buf), "Full frame %dx%d (%.0f FPS)", width,
             height, fps);
    result.name = strdup(name_buf);

    result.ops_per_sec = fps;
    result.avg_time_us = (elapsed * 1000.0) / (double) iterations;
    result.iterations = iterations;
    return result;
}

int main(void)
{
    printf(ANSI_BOLD "B3D Performance Benchmarks\n" ANSI_RESET);
    printf("===========================\n");
    printf("Each benchmark runs for ~1 second\n\n");

    bench_result_t results[16];
    int num_results = 0;

    printf(ANSI_BOLD "Primitive Operations:\n" ANSI_RESET);
    results[num_results++] = bench_matrix_ops();
    print_result(&results[num_results - 1]);

    results[num_results++] = bench_matrix_stack();
    print_result(&results[num_results - 1]);

    results[num_results++] = bench_to_screen();
    print_result(&results[num_results - 1]);

    printf("\n" ANSI_BOLD "Buffer Operations:\n" ANSI_RESET);
    results[num_results++] = bench_clear(320, 240);
    print_result(&results[num_results - 1]);

    results[num_results++] = bench_clear(640, 480);
    print_result(&results[num_results - 1]);

    printf("\n" ANSI_BOLD "Rendering Throughput:\n" ANSI_RESET);
    results[num_results++] = bench_triangles(320, 240);
    print_result(&results[num_results - 1]);

    results[num_results++] = bench_cubes(320, 240);
    print_result(&results[num_results - 1]);

    printf("\n" ANSI_BOLD "Frame Rate (clear + render):\n" ANSI_RESET);
    results[num_results++] = bench_full_frame(320, 240);
    print_result(&results[num_results - 1]);

    results[num_results++] = bench_full_frame(640, 480);
    print_result(&results[num_results - 1]);

    results[num_results++] = bench_full_frame(800, 600);
    print_result(&results[num_results - 1]);

    printf("\n===========================\n");
    printf(ANSI_GREEN "Completed %d benchmarks" ANSI_RESET "\n", num_results);

    return 0;
}
