#undef NDEBUG
#include "pycrypto_common.h"

/**
 * Double-precision multiplication
 */
#if defined(HAVE_UINT128) & 0

#define DP_MULT(a,b,ol,oh) do { \
    __uint128_t pr;             \
    pr = (__uint128_t)(a)*(b);  \
    ol = (__uint128_t)pr;       \
    oh = pr >> 64;              \
    } while (0)

#elif defined(_MSC_VER) && defined(_WIN64)

#include <windows.h>
#define DP_MULT(a,b,ol,oh) do { ol = UnsignedMultiply128(a,b,&oh); } while (0)

#else

uint64_t static inline dp_mult_128_32(uint64_t a, uint64_t b, uint64_t *oh)
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("-O3")))
#endif
;

#if 0
uint64_t static inline dp_mult_128_32(uint64_t a, uint64_t b, uint64_t *oh)
{
    uint32_t al = (uint32_t) a;
    uint32_t ah = a >> 32;
    uint32_t bl = (uint32_t) b;
    uint32_t bh = b >> 32;

    uint64_t sum0, sum1a, sum1b, sum2, sum3;

    sum0 = (uint64_t)al*bl;
    sum1a = (uint64_t)al*bh;
    sum1b = (uint64_t)ah*bl;
    sum2 = (uint64_t)ah*bh;

    sum1a += sum0 >> 32;
    sum1b += sum1a;
    sum3 = sum1b < sum1a;
    sum2 += sum1b >> 32;
    sum3 += sum2 >> 32;

    *oh = (sum3 << 32) + (uint32_t)sum2;
    return (sum1b << 32) + (uint32_t)sum0;
}

#else

#include <x86intrin.h>
#pragma GCC target "sse4.1"
uint64_t static inline dp_mult_128_32(uint64_t x, uint64_t y, uint64_t *oh)
{
    uint64_t lo;

    __m128i r0 = _mm_set1_epi32((uint32_t)(x >> 32));                           // { xH, xH, xH, xH }
    __m128i r1 = _mm_set1_epi32((uint32_t)x);                                   // { xL, xL, xL, xL }
    __m128i r2 = _mm_shuffle_epi32(_mm_cvtsi64_si128(y), _MM_SHUFFLE(2,1,2,0)); // { 0,  yH, 0,  yL }
    __m128i r3 = _mm_mul_epu32(r0, r2);
    __m128i r4 = _mm_mul_epu32(r1, r2);
    lo = (uint32_t) _mm_extract_epi32(r4, 0);      // lowest 32-bit word (TODO: replace with 2x _mm_extract_epi16)
    __m128i r5 = _mm_unpackhi_epi64(r4, _mm_setzero_si128());
    __m128i r6 = _mm_shuffle_epi32(r5, _MM_SHUFFLE(3,1,3,0));  // { 0, H(xL*yH), 0, L(xL*yH) }
    __m128i r7 = _mm_add_epi64(r3, r6); // Step i
    __m128i r8 = _mm_move_epi64(r4);
    __m128i r9 = _mm_shuffle_epi32(r8, _MM_SHUFFLE(3,3,3,1));  // { 0, 0, 0, H(xL*yL) }
    __m128i r10 = _mm_add_epi64(r7, r9); // Step ii
    lo += ((uint64_t)_mm_extract_epi32(r10, 0)) << 32;      // lowest 32-bit word 
    __m128i r11 = _mm_move_epi64(r10);
    __m128i r12 = _mm_shuffle_epi32(r11, _MM_SHUFFLE(3,1,3,3));  // { 0, H, 0, 0 }
    __m128i r13 = _mm_add_epi64(r10, r12); // Step iii
    *oh = _mm_extract_epi64(r13, 1);

    return lo;
}

#endif

#define DP_MULT(a,b,ol,oh) do { ol = dp_mult_128_32(a,b,&oh); } while (0)

#endif

size_t square_w(uint64_t *t, const uint64_t *a, size_t words)
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("-O3")))
#endif
;

uint64_t addmul128(uint64_t * RESTRICT t, const uint64_t * RESTRICT a, uint64_t b0, uint64_t b1, size_t words)
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("-O3")))
#endif
;
