/*
 * B3D is freely redistributable under the MIT License. See the file
 * "LICENSE" for information on usage and redistribution of this file.
 */

/* Internal math types and functions. */

#ifndef B3D_MATH_H
#define B3D_MATH_H

#include <limits.h>
#include <math.h>
#include <stdint.h>

/*
 * Fixed-point / scalar arithmetic
 *
 * Supports two modes:
 *   - B3D_FLOAT_POINT: Use native float (faster on FPU-equipped hardware)
 *   - Default: Q15.16 fixed-point (better for embedded without FPU)
 */

#ifndef B3D_DEPTH_16BIT
#ifndef B3D_DEPTH_32BIT
#define B3D_DEPTH_32BIT 1
#endif
#endif

#if defined(B3D_FLOAT_POINT) && defined(B3D_DEPTH_16BIT)
#error "B3D_FLOAT_POINT and B3D_DEPTH_16BIT cannot be combined"
#endif

/* Floating-point path for comparison builds */
#ifdef B3D_FLOAT_POINT

typedef float b3d_scalar_t;
#define B3D_FP_BITS 0
#define B3D_FP_ONE 1.0f
#define B3D_FP_HALF 0.5f

#define B3D_INT_TO_FP(i) ((b3d_scalar_t) (i))
#define B3D_FLOAT_TO_FP(f) ((b3d_scalar_t) (f))
#define B3D_FP_TO_INT(f) ((int) (f))
#define B3D_FP_TO_FLOAT(f) ((float) (f))

#define B3D_FP_MUL(a, b) ((a) * (b))
#define B3D_FP_DIV(a, b) ((b) == 0 ? 0 : (a) / (b))
#define B3D_FP_FLOOR(f) floorf(f)
#define B3D_FP_ADD(a, b) ((a) + (b))
#define B3D_FP_SUB(a, b) ((a) - (b))

/* Fixed-point style constants expressed in float for compatibility */
#define B3D_FP_PI 3.1415926536f
#define B3D_FP_PI_HALF 1.5707963268f
#define B3D_FP_2PI 6.2831853072f

static inline b3d_scalar_t b3d_fp_sin(b3d_scalar_t x)
{
    return sinf(x);
}

static inline b3d_scalar_t b3d_fp_cos(b3d_scalar_t x)
{
    return cosf(x);
}

static inline b3d_scalar_t b3d_fp_sqrt(b3d_scalar_t a)
{
    return a <= 0 ? 0 : sqrtf(a);
}

static inline b3d_scalar_t b3d_fp_abs(b3d_scalar_t x)
{
    return fabsf(x);
}

/* Float mode: inner loop uses same type, no conversion needed */
typedef float b3d_tri_scalar_t;
#define B3D_TRI_FP_BITS 0
#define B3D_TRI_FP_ONE 1.0f
#define B3D_TRI_FP_MUL(a, b) ((a) * (b))
#define B3D_TRI_FP_DIV(a, b) ((b) == 0 ? 0 : (a) / (b))
#define B3D_FP_TO_TRI(f) (f)
#define B3D_TRI_TO_FP(f) (f)
#define B3D_FP_MUL_TRI(fp, tri) ((fp) * (tri))
static inline b3d_tri_scalar_t b3d_tri_mul_add_sat(b3d_tri_scalar_t t,
                                                   b3d_tri_scalar_t t_step,
                                                   int skip)
{
    return t + t_step * (float) skip;
}

#else /* Fixed-point path */

/* Q15.16 format in int32_t: range ±32,768, precision 1/65536 */
typedef int32_t b3d_fixed_t;
#define B3D_FP_BITS 16
#define B3D_FP_ONE (1LL << B3D_FP_BITS) /* 1LL avoids UB */
#define B3D_FP_HALF (1LL << (B3D_FP_BITS - 1))
#define B3D_FP_EPSILON 16 /* ~0.000244 in Q15.16, for near-zero comparisons */

/* Use multiplication (not left-shift) to avoid UB with negative values */
#define B3D_INT_TO_FP(i) ((b3d_fixed_t) ((int64_t) (i) * B3D_FP_ONE))
#define B3D_FLOAT_TO_FP(f) ((b3d_fixed_t) ((f) * B3D_FP_ONE))
#define B3D_FP_TO_INT(f) ((f) >> B3D_FP_BITS)
#define B3D_FP_TO_FLOAT(f) ((float) (f) / B3D_FP_ONE)

/* int64_t intermediates required for precision.
 * Use multiplication instead of left-shift to avoid UB with negative values. */
#define B3D_FP_MUL(a, b) ((b3d_fixed_t) (((int64_t) (a) * (b)) >> B3D_FP_BITS))
#define B3D_FP_DIV(a, b) \
    ((b) == 0 ? 0 : (b3d_fixed_t) (((int64_t) (a) * B3D_FP_ONE) / (b)))
#define B3D_FP_FLOOR(f)                                  \
    ((f) & ~((1 << B3D_FP_BITS) - 1)) /* replaces floorf \
                                       */
#define B3D_FP_ADD(a, b) ((a) + (b))
#define B3D_FP_SUB(a, b) ((a) - (b))

/* Define scalar type for fixed-point */
typedef b3d_fixed_t b3d_scalar_t;

/* Fixed-point constants */
#define B3D_FP_PI ((b3d_fixed_t) ((355LL << B3D_FP_BITS) / 113))
#define B3D_FP_PI_HALF (B3D_FP_PI >> 1)
#define B3D_FP_3PI_HALF (B3D_FP_PI + B3D_FP_PI_HALF)
#define B3D_FP_2PI (B3D_FP_PI << 1)
#define B3D_FP_PI_SQ \
    ((b3d_fixed_t) (((int64_t) B3D_FP_PI * B3D_FP_PI) >> B3D_FP_BITS))

/* Platform-aware mathematical optimizations using compiler intrinsics
 * (CLZ, platform detection) for improved performance without external deps.
 *
 * Performance focus:
 * - Fast integer log2 via CLZ (count leading zeros): ~32x faster
 * - Platform-optimized sqrt (Newton-Raphson on 64-bit, bit-by-bit on 32-bit)
 * - Bit manipulation utilities
 */

/* Compiler detection and intrinsic availability */
#if defined(__GNUC__) || defined(__clang__)
#define B3D_HAS_CLZ 1
#else
#define B3D_HAS_CLZ 0
#endif

/* Platform detection for hardware 64-bit division
 *
 * Platforms WITH hardware 64-bit division (fast Newton-Raphson):
 * - x86_64 (IDIV instruction)
 * - ARM64/AArch64 (UDIV/SDIV for 64-bit)
 * - RISC-V 64-bit with M extension (DIV/DIVU)
 *
 * Platforms WITHOUT hardware 64-bit division (use bit-by-bit):
 * - ARM Cortex-M (M0/M0+/M3/M4/M7) - 32-bit only
 * - ARM Cortex-A (32-bit variants)
 * - x86 (32-bit)
 * - RISC-V 32-bit or without M extension
 */
#ifndef B3D_USE_NEWTON_SQRT
#if defined(__x86_64__) || defined(__amd64__) || defined(_M_X64) || \
    defined(__aarch64__) || defined(_M_ARM64)
/* 64-bit x86/ARM with guaranteed hardware 64-bit division */
#define B3D_USE_NEWTON_SQRT 1
#elif defined(__riscv) && __riscv_xlen == 64
/* RISC-V 64: Check for M extension (multiply/divide) */
#if defined(__riscv_mul) && defined(__riscv_div)
#define B3D_USE_NEWTON_SQRT 1
#else
/* RV64 without M extension: use bit-by-bit */
#define B3D_USE_NEWTON_SQRT 0
#endif
#else
/* 32-bit or embedded: use portable bit-by-bit method */
#define B3D_USE_NEWTON_SQRT 0
#endif
#endif

/* Fast integer log2 (floor) using count leading zeros
 *
 * For x > 0, returns floor(log2(x)).
 * For x == 0, returns -1 (undefined mathematically).
 *
 * Performance: ~1-2 cycles (CLZ instruction) vs ~50+ cycles (FPU log2)
 *
 * @param x: Input value (must be > 0 for defined result)
 * @return: floor(log2(x)), or -1 if x == 0
 */
static inline int b3d_fast_log2_u32(uint32_t x)
{
#if B3D_HAS_CLZ
    return x == 0 ? -1 : 31 - __builtin_clz(x);
#else
    /* Fallback: de Bruijn sequence method (still fast) */
    static const int debruijn32[32] = {
        0, 9,  1,  10, 13, 21, 2,  29, 11, 14, 16, 18, 22, 25, 3, 30,
        8, 12, 20, 28, 15, 17, 24, 7,  19, 27, 23, 6,  26, 5,  4, 31,
    };
    if (x == 0)
        return -1;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return debruijn32[(uint32_t) (x * 0x07C4ACDDU) >> 27];
#endif
}

static inline int b3d_fast_log2_u64(uint64_t x)
{
#if B3D_HAS_CLZ
    /* GCC/Clang handle __builtin_clzll efficiently on both 32-bit and 64-bit */
    return x == 0 ? -1 : 63 - __builtin_clzll(x);
#else
    /* Fall back to 32-bit version for portability */
    uint32_t hi = (uint32_t) (x >> 32);
    if (hi != 0)
        return 32 + b3d_fast_log2_u32(hi);
    return b3d_fast_log2_u32((uint32_t) x);
#endif
}

/* Bit manipulation utilities */

/* Sine/Cosine Lookup Table
 *
 * Power-of-two table size enables bitwise masking instead of modulo.
 * Pre-computed reciprocal eliminates expensive 64-bit division.
 * Linear interpolation provides sub-degree precision (~0.1% max error).
 *
 * Memory: 513 entries * 4 bytes = ~2KB
 * Performance: multiply+shift+mask vs. division+modulo
 */
#define B3D_TRIG_LUT_BITS 9
#define B3D_TRIG_LUT_SIZE (1 << B3D_TRIG_LUT_BITS) /* 512 entries */
#define B3D_TRIG_LUT_MASK (B3D_TRIG_LUT_SIZE - 1)  /* 0x1FF for fast wrap */

/* Pre-computed: (512 << 32) / B3D_FP_2PI for radians-to-index conversion.
 * This replaces division with multiply+shift: index = (x * RECIP) >> 32
 */
#define B3D_RAD_TO_IDX_RECIP 5340364LL

static b3d_fixed_t b3d_sin_lut[B3D_TRIG_LUT_SIZE + 1]; /* +1 for interp wrap */
static int b3d_trig_lut_ready = 0;

/* Initialize trig LUT - call from b3d_init() or on first use */
static inline void b3d_init_trig_lut(void)
{
    if (b3d_trig_lut_ready)
        return;
    for (int i = 0; i <= B3D_TRIG_LUT_SIZE; i++) {
        float rad = (float) i * (6.2831853072f / (float) B3D_TRIG_LUT_SIZE);
        b3d_sin_lut[i] = (b3d_fixed_t) (sinf(rad) * (float) B3D_FP_ONE);
    }
    b3d_trig_lut_ready = 1;
}

/* Core LUT lookup with linear interpolation (no init check - caller ensures).
 * @x: angle in Q15.16 fixed-point radians (must be non-negative)
 * @idx_offset: 0 for sine, 128 for cosine (quarter circle)
 * Returns: result in Q15.16 fixed-point
 */
static inline b3d_fixed_t b3d_lut_lookup(int64_t x64, int idx_offset)
{
    /* Convert radians to index using multiply+shift (no division!)
     * index_fp = x * (512 / 2π) in Q16.16 format
     */
    int64_t idx_fp = (x64 * B3D_RAD_TO_IDX_RECIP) >> B3D_FP_BITS;

    /* Extract integer index and fractional part */
    uint32_t idx =
        (((uint32_t) (idx_fp >> B3D_FP_BITS)) + idx_offset) & B3D_TRIG_LUT_MASK;
    b3d_fixed_t frac = (b3d_fixed_t) (idx_fp & ((1 << B3D_FP_BITS) - 1));

    /* Linear interpolation between adjacent table entries */
    b3d_fixed_t s0 = b3d_sin_lut[idx];
    b3d_fixed_t s1 = b3d_sin_lut[idx + 1];
    return s0 + B3D_FP_MUL(s1 - s0, frac);
}

/* Fast sine using power-of-two LUT with linear interpolation.
 * @x: angle in Q15.16 fixed-point radians
 * Returns: sin(x) in Q15.16 fixed-point
 */
static inline b3d_fixed_t b3d_lut_sin(b3d_fixed_t x)
{
    int sign = 1;
    int64_t x64 = x;

    if (x64 < 0) {
        x64 = -x64;
        sign = -1;
    }

    b3d_fixed_t result = b3d_lut_lookup(x64, 0);
    return sign == 1 ? result : -result;
}

/* Fast cosine using LUT: cos(x) = sin(x + π/2), offset by 128 in 512-entry
 * table
 * @x: angle in Q15.16 fixed-point radians
 * Returns: cos(x) in Q15.16 fixed-point
 */
static inline b3d_fixed_t b3d_lut_cos(b3d_fixed_t x)
{
    int64_t x64 = x < 0 ? -x : x;                      /* cos(-x) = cos(x) */
    return b3d_lut_lookup(x64, B3D_TRIG_LUT_SIZE / 4); /* +128 = 90° offset */
}

/* Compute sine and cosine together (shares conversion work) */
static inline void b3d_lut_sincos(b3d_fixed_t x,
                                  b3d_fixed_t *sinp,
                                  b3d_fixed_t *cosp)
{
    int sin_sign = 1;
    int64_t x64 = x;

    if (x64 < 0) {
        x64 = -x64;
        sin_sign = -1;
    }

    /* Single radians-to-index conversion shared by sin and cos */
    int64_t idx_fp = (x64 * B3D_RAD_TO_IDX_RECIP) >> B3D_FP_BITS;
    uint32_t sin_idx = ((uint32_t) (idx_fp >> B3D_FP_BITS)) & B3D_TRIG_LUT_MASK;
    b3d_fixed_t frac = (b3d_fixed_t) (idx_fp & ((1 << B3D_FP_BITS) - 1));

    /* Sine interpolation */
    b3d_fixed_t s0 = b3d_sin_lut[sin_idx];
    b3d_fixed_t s1 = b3d_sin_lut[sin_idx + 1];
    b3d_fixed_t sin_val = s0 + B3D_FP_MUL(s1 - s0, frac);

    /* Cosine: +128 offset (quarter circle) */
    uint32_t cos_idx = (sin_idx + B3D_TRIG_LUT_SIZE / 4) & B3D_TRIG_LUT_MASK;
    s0 = b3d_sin_lut[cos_idx];
    s1 = b3d_sin_lut[cos_idx + 1];
    b3d_fixed_t cos_val = s0 + B3D_FP_MUL(s1 - s0, frac);

    if (sinp)
        *sinp = (sin_sign == 1) ? sin_val : -sin_val;
    if (cosp)
        *cosp = cos_val;
}

/* Primary sine - initializes LUT on first use if needed.
 * @x: angle in Q15.16 fixed-point radians
 * Returns: sin(x) in Q15.16 fixed-point
 */
static inline b3d_fixed_t b3d_fp_sin(b3d_fixed_t x)
{
    if (!b3d_trig_lut_ready)
        b3d_init_trig_lut();
    return b3d_lut_sin(x);
}

/* Primary cosine - initializes LUT on first use if needed.
 * @x: angle in Q15.16 fixed-point radians
 * Returns: cos(x) in Q15.16 fixed-point
 */
static inline b3d_fixed_t b3d_fp_cos(b3d_fixed_t x)
{
    if (!b3d_trig_lut_ready)
        b3d_init_trig_lut();
    return b3d_lut_cos(x);
}

/* Compute sine and cosine together using LUT (faster than separate calls).
 * @x:    angle in Q15.16 fixed-point radians
 * @sinp: pointer to store sine result (may be NULL)
 * @cosp: pointer to store cosine result (may be NULL)
 */
static inline void b3d_fp_sincos(b3d_fixed_t x,
                                 b3d_fixed_t *sinp,
                                 b3d_fixed_t *cosp)
{
    if (!b3d_trig_lut_ready)
        b3d_init_trig_lut();
    b3d_lut_sincos(x, sinp, cosp);
}

/* Fixed-point square root (Q15.16 format)
 *
 * Platform-optimized implementation:
 * - x86_64/ARM64: Newton-Raphson (3x faster, hardware 64-bit division)
 * - Cortex-M/32-bit: Bit-by-bit algorithm (optimal, no 64-bit division)
 *
 * Platform detection is automatic via B3D_USE_NEWTON_SQRT defined above.
 * Manual override: gcc -DB3D_USE_NEWTON_SQRT=0 (force bit-by-bit)
 *
 * @a:    Input in Q15.16 fixed-point
 * Returns: floor(sqrt(a)) in Q15.16 fixed-point
 */
static inline int32_t b3d_fast_sqrt_fp16(int32_t a)
{
    if (a <= 0)
        return 0;

#if !B3D_USE_NEWTON_SQRT
    /* 32-bit/embedded: bit-by-bit method (no 64-bit division) */
    uint64_t n = ((uint64_t) (uint32_t) a) << 16;
    uint64_t res = 0;

    /* Optimization: Use CLZ to find initial bit position directly
     * This eliminates the setup loop (saves ~16 iterations)
     */
    int log2_n = b3d_fast_log2_u64(n);
    uint64_t bit = 1ULL << ((log2_n + 1) & ~1); /* Round to even power */

    while (bit != 0) {
        if (n >= res + bit) {
            n -= res + bit;
            res = (res >> 1) + bit;
        } else {
            res >>= 1;
        }
        bit >>= 2;
    }

    return res > INT32_MAX ? INT32_MAX : (int32_t) res;
#else  /* B3D_USE_NEWTON_SQRT */
    /* 64-bit: Newton-Raphson with hardware division */
    uint64_t n = ((uint64_t) (uint32_t) a) << 16;

    /* Fast initial guess using log2 */
    int log2_n = b3d_fast_log2_u64(n);
    uint64_t g = 1ULL << (log2_n >> 1);

    /* Newton-Raphson: g = (g + n/g) / 2
     * Uses hardware 64-bit division (fast on x86_64/ARM64)
     * Reduced from 5 to 4 iterations (Codex recommendation)
     */
    g = (g + n / g) >> 1;
    g = (g + n / g) >> 1;
    g = (g + n / g) >> 1;
    g = (g + n / g) >> 1;

    /* Floor correction: Ensure floor(sqrt(n)) correctness
     * Newton-Raphson can be off by ±1-2 ULPs due to truncation
     * This guarantees exact floor() result (Codex recommendation)
     */
    while ((g + 1) * (g + 1) <= n)
        g++; /* Round up if next integer is still valid */
    while (g * g > n)
        g--; /* Round down if current is too large */

    /* Clamp to 32-bit range */
    return g > INT32_MAX ? INT32_MAX : (int32_t) g;
#endif /* B3D_USE_NEWTON_SQRT */
}

#if B3D_USE_NEWTON_SQRT
/* x86_64/ARM64: Newton-Raphson with hardware 64-bit division */
static inline b3d_fixed_t b3d_fp_sqrt(b3d_fixed_t a)
{
    return b3d_fast_sqrt_fp16(a);
}

#else /* 32-bit/Cortex-M: Bit-by-bit (no 64-bit division) */

static inline b3d_fixed_t b3d_fp_sqrt(b3d_fixed_t a)
{
    if (a <= 0)
        return 0;

    /* Scale by 2^16 so sqrt preserves fixed-point fraction bits */
    uint64_t n = ((uint64_t) (uint32_t) a) << B3D_FP_BITS;
    uint64_t res = 0;
    uint64_t bit = 1ULL << 62; /* Largest power-of-four starting point */

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

    return res > INT32_MAX ? INT32_MAX : (b3d_fixed_t) res;
}

#endif /* B3D_FP_SQRT implementation */

/* Fixed-point absolute value - guards against INT32_MIN overflow */
static inline b3d_fixed_t b3d_fp_abs(b3d_fixed_t x)
{
    if (x == INT32_MIN)
        return INT32_MAX; /* -INT32_MIN overflows; clamp to max */
    return x < 0 ? -x : x;
}

/* Higher-precision inner loop format (Q12.20).
 *
 * Screen coordinates are bounded (typically 0-1920), so we can trade
 * integer range for fractional precision. Q12.20 provides 16x more
 * precision than Q15.16 for smoother edge interpolation.
 */
typedef int32_t b3d_tri_scalar_t;
#define B3D_TRI_FP_BITS 20
#define B3D_TRI_FP_ONE (1LL << B3D_TRI_FP_BITS)
#define B3D_TRI_FP_MUL(a, b) \
    ((b3d_tri_scalar_t) (((int64_t) (a) * (b)) >> B3D_TRI_FP_BITS))
#define B3D_TRI_FP_DIV(a, b) \
    ((b) == 0 ? 0 : (b3d_tri_scalar_t) (((int64_t) (a) * B3D_TRI_FP_ONE) / (b)))

/* Convert Q15.16 to Q12.20 with saturation.
 * Q12.20 has 12 integer bits, so values >= 2048 or <= -2048 saturate.
 * Uses inline function for overflow protection.
 */
static inline b3d_tri_scalar_t b3d_fp_to_tri(b3d_fixed_t f)
{
    int64_t wide = (int64_t) f << (B3D_TRI_FP_BITS - B3D_FP_BITS);
    if (wide > INT32_MAX)
        return INT32_MAX;
    if (wide < INT32_MIN)
        return INT32_MIN;
    return (b3d_tri_scalar_t) wide;
}
#define B3D_FP_TO_TRI(f) b3d_fp_to_tri(f)

/* Convert Q12.20 to Q15.16: shift right 4 bits */
#define B3D_TRI_TO_FP(f) \
    ((b3d_scalar_t) ((f) >> (B3D_TRI_FP_BITS - B3D_FP_BITS)))
/* Mixed multiply: Q15.16 * Q12.20 → Q15.16 (shift by TRI bits) */
#define B3D_FP_MUL_TRI(fp, tri) \
    ((b3d_scalar_t) (((int64_t) (fp) * (tri)) >> B3D_TRI_FP_BITS))

/* Saturating multiply-add for skip accumulation: t += t_step * skip.
 * Prevents overflow when skip is large (off-screen geometry).
 */
static inline b3d_tri_scalar_t b3d_tri_mul_add_sat(b3d_tri_scalar_t t,
                                                   b3d_tri_scalar_t t_step,
                                                   int skip)
{
    int64_t delta = (int64_t) t_step * skip;
    int64_t result = (int64_t) t + delta;
    if (result > INT32_MAX)
        return INT32_MAX;
    if (result < INT32_MIN)
        return INT32_MIN;
    return (b3d_tri_scalar_t) result;
}

#endif /* B3D_FLOAT_POINT */

/*
 * Float-based math wrappers for public API and examples.
 * These provide a unified interface regardless of internal fixed/float mode.
 * Define guard to prevent redefinition if b3d.h is included later.
 */
#ifndef B3D_MATH_UTILS_DEFINED
#define B3D_MATH_UTILS_DEFINED

static inline float b3d_sinf(float x)
{
#ifdef B3D_FLOAT_POINT
    return sinf(x);
#else
    return B3D_FP_TO_FLOAT(b3d_fp_sin(B3D_FLOAT_TO_FP(x)));
#endif
}

static inline float b3d_cosf(float x)
{
#ifdef B3D_FLOAT_POINT
    return cosf(x);
#else
    return B3D_FP_TO_FLOAT(b3d_fp_cos(B3D_FLOAT_TO_FP(x)));
#endif
}

static inline float b3d_sqrtf(float x)
{
#ifdef B3D_FLOAT_POINT
    return x <= 0.0f ? 0.0f : sqrtf(x);
#else
    return B3D_FP_TO_FLOAT(b3d_fp_sqrt(B3D_FLOAT_TO_FP(x)));
#endif
}

static inline float b3d_fabsf(float x)
{
#ifdef B3D_FLOAT_POINT
    return fabsf(x);
#else
    return x < 0.0f ? -x : x;
#endif
}

static inline float b3d_tanf(float x)
{
#ifdef B3D_FLOAT_POINT
    return tanf(x);
#else
    float c = b3d_cosf(x);
    return (c < 1e-7f && c > -1e-7f) ? 0.0f : b3d_sinf(x) / c;
#endif
}

/* Compute sine and cosine simultaneously (more efficient than separate calls)
 */
static inline void b3d_sincosf(float x, float *sinp, float *cosp)
{
#ifdef B3D_FLOAT_POINT
#if defined(__GLIBC__) && defined(_GNU_SOURCE)
    sincosf(x, sinp, cosp);
#else
    *sinp = sinf(x);
    *cosp = cosf(x);
#endif
#else
    b3d_scalar_t angle = B3D_FLOAT_TO_FP(x);
    b3d_fixed_t sin_fp, cos_fp;
    b3d_fp_sincos(angle, &sin_fp, &cos_fp);
    *sinp = B3D_FP_TO_FLOAT(sin_fp);
    *cosp = B3D_FP_TO_FLOAT(cos_fp);
#endif
}

#endif /* B3D_MATH_UTILS_DEFINED */

/*
 * Constants
 */
#define B3D_NEAR_DISTANCE 0.1f
#define B3D_FAR_DISTANCE 100.0f
#define B3D_EPSILON 1e-8f /* Near-zero threshold for division guards */
#define B3D_DEGEN_THRESHOLD 0.0001f /* Degenerate trig/scanline threshold */
#define B3D_CULL_THRESHOLD 0.01f /* Back-face culling dot product threshold */
#define B3D_DEPTH_FAR 1e30f      /* Depth buffer clear value (far plane) */
#define B3D_PI 3.1415926536f     /* Pi constant for angle conversions */
#define B3D_CLIP_BUFFER_SIZE 32  /* Maximum triangles in clipping buffer */

/*
 * Depth buffer type (shared with public header)
 */
#ifndef B3D_DEPTH_T_DEFINED
#ifdef B3D_FLOAT_POINT
typedef float b3d_depth_t;
#elif defined(B3D_DEPTH_16BIT)
typedef uint16_t b3d_depth_t;
#else
typedef int32_t b3d_depth_t;
#endif
#define B3D_DEPTH_T_DEFINED
#endif

#if defined(B3D_FLOAT_POINT)
#define B3D_DEPTH_CLEAR B3D_DEPTH_FAR
static inline b3d_depth_t b3d_depth_from_float(float d)
{
    return d;
}
static inline float b3d_depth_to_float(b3d_depth_t d)
{
    return d;
}
#elif defined(B3D_DEPTH_16BIT)
#define B3D_DEPTH_CLEAR 0xFFFF
static inline b3d_depth_t b3d_depth_from_float(float d)
{
    if (d < 0.0f)
        d = 0.0f;
    if (d > 1.0f)
        d = 1.0f;
    return (b3d_depth_t) (d * 65535.0f + 0.5f);
}
static inline float b3d_depth_to_float(b3d_depth_t d)
{
    return (float) d * (1.0f / 65535.0f);
}
#else
#define B3D_DEPTH_CLEAR INT32_MAX
static inline b3d_depth_t b3d_depth_from_float(float d)
{
    return B3D_FLOAT_TO_FP(d);
}
static inline float b3d_depth_to_float(b3d_depth_t d)
{
    return B3D_FP_TO_FLOAT(d);
}
#endif

/*
 * Internal vector/matrix types
 */
typedef struct {
    float x, y, z, w;
} b3d_vec_t;

typedef struct {
    float m[4][4];
} b3d_mat_t;

typedef struct {
    b3d_vec_t p[3];
} b3d_triangle_t;

/*
 * Vector and matrix operations - generated from src/math.dsl
 * Provides vector/matrix functions (always float-based since b3d_vec_t uses
 * float)
 */
#include "math-gen.inc"

/* Note: Most math functions are generated from DSL (src/math.dsl) */

/*
 * Geometry operations
 */
static inline int b3d_clip_against_plane(b3d_vec_t plane,
                                         b3d_vec_t norm,
                                         b3d_triangle_t in,
                                         b3d_triangle_t out[2])
{
    norm = b3d_vec_norm(norm);
    float plane_d = b3d_vec_dot(norm, plane);
    b3d_vec_t *inside[3];
    int inside_count = 0;
    b3d_vec_t *outside[3];
    int outside_count = 0;
    float d0 = b3d_vec_dot(in.p[0], norm) - plane_d;
    float d1 = b3d_vec_dot(in.p[1], norm) - plane_d;
    float d2 = b3d_vec_dot(in.p[2], norm) - plane_d;
    if (d0 >= 0)
        inside[inside_count++] = &in.p[0];
    else
        outside[outside_count++] = &in.p[0];
    if (d1 >= 0)
        inside[inside_count++] = &in.p[1];
    else
        outside[outside_count++] = &in.p[1];
    if (d2 >= 0)
        inside[inside_count++] = &in.p[2];
    else
        outside[outside_count++] = &in.p[2];
    if (inside_count == 3) {
        out[0] = in;
        return 1;
    } else if (inside_count == 1 && outside_count == 2) {
        out[0].p[0] = *inside[0];
        out[0].p[1] =
            b3d_intersect_plane(norm, plane_d, *inside[0], *outside[0]);
        out[0].p[2] =
            b3d_intersect_plane(norm, plane_d, *inside[0], *outside[1]);
        return 1;
    } else if (inside_count == 2 && outside_count == 1) {
        out[0].p[0] = *inside[0];
        out[0].p[1] = *inside[1];
        out[0].p[2] =
            b3d_intersect_plane(norm, plane_d, *inside[0], *outside[0]);
        out[1].p[0] = *inside[1];
        out[1].p[1] = out[0].p[2];
        out[1].p[2] =
            b3d_intersect_plane(norm, plane_d, *inside[1], *outside[0]);
        return 2;
    }
    return 0;
}

#endif /* B3D_MATH_H */
