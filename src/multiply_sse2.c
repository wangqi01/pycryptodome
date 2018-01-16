/*
Copyright (c) 2017, Helder Eijs <helderijs@gmail.com>
All rights reserved. 

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met: 

 * Redistributions of source code must retain the above copyright notice, 
   this list of conditions and the following disclaimer. 
 * Redistributions in binary form must reproduce the above copyright 
   notice, this list of conditions and the following disclaimer in the 
   documentation and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND ANY 
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY 
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY 
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
DAMAGE. 
*/

#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "multiply.h"

#if defined(__GNUC__)
// Note: this pragma does not work with clang
#pragma GCC target("sse2")
#endif


/**
 * Add a 64-bit value x to y/sum_mid/sum_hi
 */
#if defined(_WIN64) && (_MSC_VER>=1900)

#include <intrin.h>
#define ADD192(y, x) do {           \
    unsigned char c = 0;            \
    c = _addcarry_u64(c, x, y, &y); \
    c = _addcarry_u64(c, 0, sum_mid, &sum_mid); \
    _addcarry_u64(c, 0, sum_hi, &sum_hi);   \
    } while (0)

#else

#define ADD192(y, x) do {       \
    uint64_t c;                 \
    y += x;                     \
    sum_mid += (c = (y < x));   \
    sum_hi += sum_mid < c;      \
    } while (0)

#endif

#if defined(_MSC_VER) && !defined(_WIN64) && (_MSC_VER<1900)

__m128i static inline _mm_set_epi64x(uint64_t e1, uint64_t e0)
{
	union{
        uint64_t u64;
        uint32_t u32[2];
    } u1, u0;
    u1.u64 = e1;
    u0.u64 = e0;
    
	return _mm_set_epi32(u1.u32[1], u1.u32[0], u0.u32[1], u0.u32[0]);
}

__m128i static inline _mm_set1_epi64x(uint64_t e0)
{
	union{
        uint64_t u64;
        uint32_t u32[2];
    } u1, u0;
    u1.u64 = e0;
    u0.u64 = e0;
    
	return _mm_set_epi32(u1.u32[1], u1.u32[0], u0.u32[1], u0.u32[0]);
}
#endif

void static inline print128(const char *str, __m128i x)
{
    uint64_t hi, lo;

    _mm_storeh_pd((double*)&hi, _mm_castsi128_pd(x));
    _mm_storel_pd((double*)&lo, _mm_castsi128_pd(x));

    printf("%s: %016lX %016lX\n", str, hi, lo);
}

int static inline get_carry_bits(__m128i x, __m128i y, __m128i sum, __m128i *all_ones)
{
    __m128i r0, r1, r2, r3, r4;
    int carries;

    // Carry check can be done efficiently with a 64-bit compare
    // Unfortunately, that requires SSE4.1

    r0 = _mm_xor_si128(sum, *all_ones);  // ~sum
    r1 = _mm_or_si128(x, y);            // x | y
    r2 = _mm_and_si128(x, y);           // x & y
    r3 = _mm_and_si128(r1, r0);         // (x | y) & ~sum
    r4 = _mm_or_si128(r2, r3);          // (x & y) | ((x & y) & ~sum)

    carries = _mm_movemask_pd(_mm_castsi128_pd(r4));

    return carries;
}

__m128i static inline mult64_128(uint64_t x, uint64_t y)
{
    __m128i r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13, r14, r15;

    r0 = _mm_set1_epi32((uint32_t)(x >> 32));                           // { xH, xH, xH, xH }
    r1 = _mm_set1_epi32((uint32_t)x);                                   // { xL, xL, xL, xL }
    r2 = _mm_shuffle_epi32(_mm_castpd_si128(_mm_set_sd(*(double*)&y)), _MM_SHUFFLE(2,1,2,0));   // { 0, yH, 0, yL }
    r3 = _mm_mul_epu32(r0, r2);                         // { xH*yH, xH*yL }
    r4 = _mm_mul_epu32(r1, r2);                         // { xL*yH, xL*yL }
    r5 = _mm_unpackhi_epi64(r4, _mm_setzero_si128());   // { xL*yH, 0 }
    r6 = _mm_shuffle_epi32(r5, _MM_SHUFFLE(3,1,3,0));   // { 0, H(xL*yH), 0, L(xL*yH) }
    r7 = _mm_add_epi64(r3, r6); // Step i
    r8 = _mm_move_epi64(r4);                            // { 0, xL*yL }
    r9 = _mm_shuffle_epi32(r8, _MM_SHUFFLE(3,3,3,1));   // { 0, 0, 0, H(xL*yL) }
    r10 = _mm_add_epi64(r7, r9); // Step ii
    r11 = _mm_move_epi64(r10);                          // { 0, L(r10) }
    r12 = _mm_shuffle_epi32(r11, _MM_SHUFFLE(3,1,3,3)); // { 0, H, 0, 0 }
    r13 = _mm_add_epi64(r10, r12); // Step iii

    //_mm_store_ss((float*)&lo, _mm_castsi128_ps(r4));
    //_mm_store_ss((float*)&lo+1, _mm_castsi128_ps(r10));
    //_mm_storeh_pi((__m64*)oh, _mm_castsi128_ps(r13));
    
    r14 = _mm_unpacklo_epi32(r4, r10);                  // { *, *, r10[31:0], r4[31:0] }
    r15 = _mm_castpd_si128(
            _mm_move_sd(
                _mm_castsi128_pd(r13),
                _mm_castsi128_pd(r14)
            )
          );
 
    return r15;
}

uint64_t addmul128(uint64_t * RESTRICT t, const uint64_t * RESTRICT a, uint64_t b0, uint64_t b1, size_t words)
{
    uint64_t sum_low, sum_mid,
             sum_hi /** < 6 **/;
    uint64_t pr_low, pr_high, aim1;
    __m128i sum_ml;
    size_t i;
    __m128i all_ones;

    if (words == 0) {
        return 0;
    }

    /** LSW **/
    DP_MULT(a[0], b0, sum_low, sum_mid);
    sum_hi = 0;
    
    ADD192(t[0], sum_low);

    sum_low = sum_mid;
    sum_mid = sum_hi;
    sum_hi = 0;

    aim1 = a[0];
    all_ones = _mm_set_epi64x(-1, -1); 
    sum_ml = _mm_set_epi64x(sum_mid, sum_low);
    for (i=1; i<words;i++) {
        int c1, c2, c3;
        int carries;
        __m128i r20, r22;
        __m128i r40, r42;
        __m128i r50, r52;
        __m128i r60, r62;

        
        r20 = mult64_128(aim1, b1);

        // --------
        // sum_low += pr_low; c1 = sum_low < pr_low;
        // sum_mid += pr_high; c2 = sum_mid < pr_high;
        
        r22 = _mm_add_epi64(r20, sum_ml);  // sum (two 64 bit words)

        carries = get_carry_bits(r20, sum_ml, r22, &all_ones);
        c1 = carries & 1;
        c2 = carries >> 1;

        sum_ml = r22;
    
        // --------

        r40 = mult64_128(a[i], b0);
        
        // --------
        // sum_low += pr_low; c1 += sum_low < pr_low;
        // sum_mid += pr_high; c2 += sum_mid < pr_high;

        r42 = _mm_add_epi64(r40, sum_ml);  // sum (two 64 bit words)

        carries = get_carry_bits(r40, sum_ml, r42, &all_ones);
        c1 += carries & 1;
        c2 += carries >> 1;

        sum_ml = r42;
   
        // --------
        // sum_low += t[i]; c3 = sum_low < t[i];
        // sum_mid += c1; c2 += sum_mid < c1;
        // t[i] = sum_low;
        
        r50 = _mm_set_epi64x((uint64_t)c1, t[i]); 
        r52 = _mm_add_epi64(r50, sum_ml);  // sum (two 64 bit words)

        carries = get_carry_bits(r50, sum_ml, r52, &all_ones);
        c3 = carries & 1;
        c2 += carries >> 1;
    
        _mm_storel_pd((double*)&t[i], _mm_castsi128_pd(r52));

        sum_ml = r52;

        // --------
        // sum_mid += c3; c2 += sum_mid < c3;
        // sum_low = sum_mid;
        // sum_mid = c2;
    
        r60 = _mm_set1_epi64x((uint64_t)c3);
        r62 = _mm_add_epi64(r60, sum_ml);  // sum (trash lower half)
        carries = get_carry_bits(r60, sum_ml, r62, &all_ones);
        c2 += carries >> 1;

        // sum_ml = { c2, H(r62) }
        
        sum_ml = _mm_set1_epi64x((uint64_t)c2);
        sum_ml = _mm_castps_si128(
                    _mm_shuffle_ps( _mm_castsi128_ps(r62),
                                    _mm_castsi128_ps(sum_ml),
                                    _MM_SHUFFLE(1,0,3,2)
                    )
                );
        
        // --------
        
        aim1 = a[i];
    }
    _mm_storel_pd((double*)&sum_low, _mm_castsi128_pd(sum_ml));
    _mm_storeh_pd((double*)&sum_mid, _mm_castsi128_pd(sum_ml));

    /** MSW **/
    DP_MULT(a[i-1], b1, pr_low, pr_high);
    ADD192(sum_low, pr_low);
    sum_mid += pr_high;
    sum_hi += sum_mid < pr_high;
    ADD192(t[i], sum_low);

    i++;
    
    sum_low = sum_mid;
    sum_mid = sum_hi;
    sum_hi = 0;
 
    /** Extend carry indefinetly **/
    for (; sum_low || sum_mid; i++) {
        ADD192(t[i], sum_low);

        sum_low = sum_mid;
        sum_mid = sum_hi;
        sum_hi = 0;
    }
    
    return i;
}


size_t square_w(uint64_t *t, const uint64_t *a, size_t words)
{
    size_t i, j;
    uint64_t carry;

    if (words == 0) {
        return 0;
    }

    memset(t, 0, 2*sizeof(uint64_t)*words);

    /** Compute all mix-products without doubling **/
    for (i=0; i<words; i++) {
        carry = 0;
        
        for (j=i+1; j<words; j++) {
            uint64_t sum_lo, sum_hi;

            DP_MULT(a[j], a[i], sum_lo, sum_hi);

            sum_lo += carry;
            sum_hi += sum_lo < carry;

            t[i+j] += sum_lo;
            carry = sum_hi + (t[i+j] < sum_lo);
        }

        /** Propagate carry **/
        for (j=i+words; carry>0; j++) {
            t[j] += (uint64_t)carry;
            carry = t[j] < carry;
        }
    }

    /** Double mix-products and add squares **/
    carry = 0;
    for (i=0, j=0; i<words; i++, j+=2) {
        uint64_t sum_lo, sum_hi, tmp, tmp2;

        DP_MULT(a[i], a[i], sum_lo, sum_hi);

        sum_lo += carry;
        sum_hi += sum_lo < carry;

        sum_hi += (tmp = ((t[j+1] << 1) + (t[j] >> 63)));
        carry = (t[j+1] >> 63) + (sum_hi < tmp);

        sum_lo += (tmp = (t[j] << 1));
        sum_hi += (tmp2 = (sum_lo < tmp));
        carry += sum_hi < tmp2;
 
        t[j] = sum_lo;
        t[j+1] = sum_hi;
    }
    assert(carry == 0);

    return 2*words;
}

