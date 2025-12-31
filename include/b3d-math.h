/*
 * B3D Math Utilities
 *
 * Unified interface for float math operations. These work in both
 * fixed-point mode (default) and floating-point mode (B3D_FLOAT_POINT).
 * Use these instead of raw sinf/cosf/etc. for consistency across all
 * B3D examples.
 */

#ifndef B3D_PUBLIC_MATH_H
#define B3D_PUBLIC_MATH_H

#include <limits.h>
#include <math.h>
#include <stdint.h>

#ifndef B3D_MATH_UTILS_DEFINED
#define B3D_MATH_UTILS_DEFINED

#ifdef B3D_FLOAT_POINT

/* Floating-point mode: use standard math library */

static inline float b3d_sinf(float x)
{
    return sinf(x);
}

static inline float b3d_cosf(float x)
{
    return cosf(x);
}

static inline float b3d_tanf(float x)
{
    return tanf(x);
}

static inline float b3d_sqrtf(float x)
{
    return x <= 0.0f ? 0.0f : sqrtf(x);
}

static inline float b3d_fabsf(float x)
{
    return x < 0.0f ? -x : x;
}

/* Compute sine and cosine simultaneously (more efficient than separate calls)
 * @x: angle in radians
 * @sinp: pointer to store sine result
 * @cosp: pointer to store cosine result
 */
static inline void b3d_sincosf(float x, float *sinp, float *cosp)
{
#if defined(__GLIBC__) && defined(_GNU_SOURCE)
    sincosf(x, sinp, cosp);
#else
    *sinp = sinf(x);
    *cosp = cosf(x);
#endif
}

#else /* Fixed-point mode (default) */

/*
 * Q15.16 fixed-point implementation using lookup table (LUT).
 * Pre-computes 360 sine values at 1-degree resolution for fast lookup.
 * Uses linear interpolation for sub-degree precision (~0.1% max error).
 */

typedef int32_t b3d_math_fixed_t;

#define B3D_MATH_FP_BITS 16
#define B3D_MATH_FP_ONE (1LL << B3D_MATH_FP_BITS)

#define B3D_MATH_FLOAT_TO_FP(f) ((b3d_math_fixed_t) ((f) * B3D_MATH_FP_ONE))
#define B3D_MATH_FP_TO_FLOAT(f) ((float) (f) / B3D_MATH_FP_ONE)
#define B3D_MATH_FP_MUL(a, b) \
    ((b3d_math_fixed_t) (((int64_t) (a) * (b)) >> B3D_MATH_FP_BITS))

/* Fixed-point pi constants using rational approximation 355/113 */
#define B3D_MATH_FP_PI ((b3d_math_fixed_t) ((355LL << B3D_MATH_FP_BITS) / 113))
#define B3D_MATH_FP_2PI (B3D_MATH_FP_PI << 1)

/*
 * Sine/Cosine Lookup Table (nGL-inspired, optimized per Gemini review)
 *
 * Power-of-two table size enables bitwise masking instead of modulo.
 * Pre-computed reciprocal eliminates expensive 64-bit division.
 * Linear interpolation provides sub-degree precision (~0.1% max error).
 *
 * Memory: 513 entries * 4 bytes = ~2KB
 */
#define B3D_MATH_LUT_BITS 9
#define B3D_MATH_LUT_SIZE (1 << B3D_MATH_LUT_BITS) /* 512 entries */
#define B3D_MATH_LUT_MASK (B3D_MATH_LUT_SIZE - 1)  /* 0x1FF for fast wrap */

/* Pre-computed: (512 << 32) / B3D_MATH_FP_2PI for radians-to-index */
#define B3D_MATH_RAD_TO_IDX 5340364LL

static b3d_math_fixed_t b3d_math_sin_lut[B3D_MATH_LUT_SIZE + 1];
static int b3d_math_lut_ready = 0;

/* Initialize LUT using float sinf() - called once on first use */
static inline void b3d_math_init_lut(void)
{
    if (b3d_math_lut_ready)
        return;
    for (int i = 0; i <= B3D_MATH_LUT_SIZE; i++) {
        float rad = (float) i * (6.2831853072f / (float) B3D_MATH_LUT_SIZE);
        b3d_math_sin_lut[i] =
            (b3d_math_fixed_t) (sinf(rad) * (float) B3D_MATH_FP_ONE);
    }
    b3d_math_lut_ready = 1;
}

/* Core LUT lookup with linear interpolation */
static inline b3d_math_fixed_t b3d_math_lut_lookup(int64_t x64, int idx_offset)
{
    /* Convert radians to index using multiply+shift (no division!) */
    int64_t idx_fp = (x64 * B3D_MATH_RAD_TO_IDX) >> B3D_MATH_FP_BITS;

    /* Extract integer index and fractional part */
    uint32_t idx = (((uint32_t) (idx_fp >> B3D_MATH_FP_BITS)) + idx_offset) &
                   B3D_MATH_LUT_MASK;
    b3d_math_fixed_t frac =
        (b3d_math_fixed_t) (idx_fp & ((1 << B3D_MATH_FP_BITS) - 1));

    /* Linear interpolation */
    b3d_math_fixed_t s0 = b3d_math_sin_lut[idx];
    b3d_math_fixed_t s1 = b3d_math_sin_lut[idx + 1];
    return s0 + B3D_MATH_FP_MUL(s1 - s0, frac);
}

/* Internal: LUT sine with linear interpolation */
static inline b3d_math_fixed_t b3d_math_lut_sin(b3d_math_fixed_t x)
{
    int sign = 1;
    int64_t x64 = x;

    if (x64 < 0) {
        x64 = -x64;
        sign = -1;
    }

    b3d_math_fixed_t result = b3d_math_lut_lookup(x64, 0);
    return sign == 1 ? result : -result;
}

/* Internal: LUT cosine (sin + 90° offset) */
static inline b3d_math_fixed_t b3d_math_lut_cos(b3d_math_fixed_t x)
{
    int64_t x64 = x < 0 ? -x : x; /* cos(-x) = cos(x) */
    return b3d_math_lut_lookup(x64, B3D_MATH_LUT_SIZE / 4); /* +128 = 90° */
}

/* Primary sine using LUT (auto-initializes if needed) */
static inline b3d_math_fixed_t b3d_math_fp_sin(b3d_math_fixed_t x)
{
    if (!b3d_math_lut_ready)
        b3d_math_init_lut();
    return b3d_math_lut_sin(x);
}

/* Primary cosine using LUT (auto-initializes if needed) */
static inline b3d_math_fixed_t b3d_math_fp_cos(b3d_math_fixed_t x)
{
    if (!b3d_math_lut_ready)
        b3d_math_init_lut();
    return b3d_math_lut_cos(x);
}

/* Compute sine and cosine together using LUT */
static inline void b3d_math_fp_sincos(b3d_math_fixed_t x,
                                      b3d_math_fixed_t *sinp,
                                      b3d_math_fixed_t *cosp)
{
    if (!b3d_math_lut_ready)
        b3d_math_init_lut();

    int sin_sign = 1;
    int64_t x64 = x;
    if (x64 < 0) {
        x64 = -x64;
        sin_sign = -1;
    }

    /* Single radians-to-index conversion shared by sin and cos */
    int64_t idx_fp = (x64 * B3D_MATH_RAD_TO_IDX) >> B3D_MATH_FP_BITS;
    uint32_t sin_idx = ((uint32_t) (idx_fp >> B3D_MATH_FP_BITS)) &
                       B3D_MATH_LUT_MASK;
    b3d_math_fixed_t frac =
        (b3d_math_fixed_t) (idx_fp & ((1 << B3D_MATH_FP_BITS) - 1));

    /* Sine interpolation */
    b3d_math_fixed_t s0 = b3d_math_sin_lut[sin_idx];
    b3d_math_fixed_t s1 = b3d_math_sin_lut[sin_idx + 1];
    b3d_math_fixed_t sin_val = s0 + B3D_MATH_FP_MUL(s1 - s0, frac);

    /* Cosine: +128 offset (quarter circle) */
    uint32_t cos_idx = (sin_idx + B3D_MATH_LUT_SIZE / 4) & B3D_MATH_LUT_MASK;
    s0 = b3d_math_sin_lut[cos_idx];
    s1 = b3d_math_sin_lut[cos_idx + 1];
    b3d_math_fixed_t cos_val = s0 + B3D_MATH_FP_MUL(s1 - s0, frac);

    if (sinp)
        *sinp = (sin_sign == 1) ? sin_val : -sin_val;
    if (cosp)
        *cosp = cos_val;
}

/* Integer sqrt on Q16.16: computes floor(sqrt(a)) in fixed-point */
static inline b3d_math_fixed_t b3d_math_fp_sqrt(b3d_math_fixed_t a)
{
    if (a <= 0)
        return 0;

    uint64_t n = ((uint64_t) (uint32_t) a) << B3D_MATH_FP_BITS;
    uint64_t res = 0;
    uint64_t bit = 1ULL << 62;

    while (bit > n)
        bit >>= 2;

    while (bit != 0) {
        if (n >= res + bit) {
            n -= res + bit;
            res = (res >> 1) + bit;
        } else {
            res >>= 1;
        }
        bit >>= 2;
    }

    return res > INT32_MAX ? INT32_MAX : (b3d_math_fixed_t) res;
}

/* Public float wrappers using fixed-point internally */

static inline float b3d_sinf(float x)
{
    return B3D_MATH_FP_TO_FLOAT(b3d_math_fp_sin(B3D_MATH_FLOAT_TO_FP(x)));
}

static inline float b3d_cosf(float x)
{
    return B3D_MATH_FP_TO_FLOAT(b3d_math_fp_cos(B3D_MATH_FLOAT_TO_FP(x)));
}

static inline float b3d_tanf(float x)
{
    float c = b3d_cosf(x);
    return (c < 1e-7f && c > -1e-7f) ? 0.0f : b3d_sinf(x) / c;
}

static inline float b3d_sqrtf(float x)
{
    return B3D_MATH_FP_TO_FLOAT(b3d_math_fp_sqrt(B3D_MATH_FLOAT_TO_FP(x)));
}

static inline float b3d_fabsf(float x)
{
    return x < 0.0f ? -x : x;
}

/* Compute sine and cosine simultaneously (more efficient than separate calls)
 * @x: angle in radians
 * @sinp: pointer to store sine result
 * @cosp: pointer to store cosine result
 */
static inline void b3d_sincosf(float x, float *sinp, float *cosp)
{
    b3d_math_fixed_t angle = B3D_MATH_FLOAT_TO_FP(x);
    b3d_math_fixed_t sin_fp, cos_fp;
    b3d_math_fp_sincos(angle, &sin_fp, &cos_fp);
    *sinp = B3D_MATH_FP_TO_FLOAT(sin_fp);
    *cosp = B3D_MATH_FP_TO_FLOAT(cos_fp);
}

#endif /* B3D_FLOAT_POINT */

#endif /* B3D_MATH_UTILS_DEFINED */

#endif /* B3D_PUBLIC_MATH_H */
