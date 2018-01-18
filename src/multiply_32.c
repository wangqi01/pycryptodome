#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "multiply.h"

void print128(const char* str, __m128i x)
{
    uint32_t *y = (uint32_t*)&x;

    printf("%s = %08lX %08X %08X %08X\n", str, y[3], y[2], y[1], y[0]);
}

void addmul32(uint32_t* t, const uint32_t *a, uint32_t b, size_t words)
{
    uint32_t carry = 0;
    int i = 0;

#if 0
    for (; i<(words & ~1); i++) {
        uint64_t prod;
        uint32_t prodl, prodh;

        prod = (uint64_t)a[i]*b;
        prodl = (uint32_t)prod;
        prodh = (uint32_t)(prod >> 32);

        prodl += carry; prodh += prodl < carry;
        t[i] += prodl; prodh += t[i] < prodl;
        carry = prodh;
    }
#else

    // assume a[] is 64-bit aligned?
    
    __m128i r0;
    __m128i r1;
   
    //printf("b = %08X\n", b); 
    r0 = _mm_set1_epi32(b);    // { b, b, b, b }
    //print128("r0", r0);
    r1 = _mm_cvtsi32_si128(carry);  // { 0, 0, 0, carry }
    for (i=0; i<(words & ~1); i+=2) {
        __m128i r10, r11, r12, r13, r14, r15, r16, r17;
        uint64_t prod;
        uint32_t prodl, prodh;
        uint32_t new_ti, new_tip1;
        uint32_t int_carry, new_carry;
       
#if 0 
        // Expectation 
        prod = (uint64_t)a[i]*b;
        prodl = (uint32_t)prod;
        prodh = (uint32_t)(prod >> 32);

        prodl += carry; prodh += prodl < carry;
        new_ti = prodl + t[i]; prodh += new_ti < prodl;
        int_carry = prodh;
        
        prod = (uint64_t)a[i+1]*b;
        prodl = (uint32_t)prod;
        prodh = (uint32_t)(prod >> 32);

        prodl += int_carry; prodh += prodl < int_carry;
        new_tip1 = prodl + t[i+1]; prodh += new_tip1 < prodl;
        new_carry = prodh;
        // End
#endif

        //printf("a[%d] = %08X\n", i, a[i]); 
        //printf("a[%d] = %08X\n", i+1, a[i+1]); 
        //printf("t[%d] = %08X\n", i, t[i]); 
        //printf("t[%d] = %08X\n", i+1, t[i+1]); 

        //print128("r1", r1);
        
        r10 = _mm_shuffle_epi32(
                _mm_castpd_si128(
                    _mm_set_sd(*(double*)&a[i])
                ),
             _MM_SHUFFLE(2,1,2,0));     // { 0, a[i+1], 0, a[i] }
        //print128("r10", r10);
        r11 = _mm_mul_epu32(r0,  r10);  // { a[i+1]*b,  a[i]*b  }
        //print128("r11", r11);
        
        r12 = _mm_shuffle_epi32(
                _mm_castpd_si128(
                    _mm_set_sd(*(double*)&t[i])
                ),
             _MM_SHUFFLE(2,1,2,0));     // { 0, t[i+1], 0, t[i] }
        //print128("r12", r12);
        r13 = _mm_add_epi64(r12, r1);   // { t[i+1],  t[i]+carry }
        //print128("r13", r13);
        r14 = _mm_add_epi64(r11, r13);  // { a[i+1]*b+t[i+1],  a[i]*b+t[i]+carry }
        //print128("r14", r14);
        r15 = _mm_shuffle_epi32(
                _mm_move_epi64(r14),    // { 0, a[i]*b+t[i]+carry }
                _MM_SHUFFLE(2,1,2,2)
              );                        // { 0, H(a[i]*b+t[i]+carry), 0, 0 }
        //print128("r15", r15);
        r16 = _mm_add_epi64(r14, r15);  // { next_carry, new t[i+1], *, new t[i] }
        //print128("r16", r16);
        r17 = _mm_shuffle_epi32(r16, _MM_SHUFFLE(2,0,1,3));
                                        // { new t[i+1], new t[i], *, new carry }

        _mm_storeh_pd((double*)&t[i],
                      _mm_castsi128_pd(r17)); // Store upper 64 bit word (also t[i+1])

        //print128("r17", r17);
        r1 = _mm_castps_si128(_mm_move_ss(
                _mm_castsi128_ps(r1),
                _mm_castsi128_ps(r17)
                ));
        //print128("new r1", r1);
        
        //printf("new t[%d] = %08X\n", i, t[i]); 
        //printf("new t[%d] = %08X\n", i+1, t[i+1]); 
    
        //printf("new carry = %08X, exp = %08X\n", carry, new_carry); 

#if 0
        assert(carry == new_carry);
        assert(new_ti == t[i]);
        assert(new_tip1 == t[i+1]);
#endif

        //assert(i<30);
        //assert(0);
    }
    carry = _mm_cvtsi128_si32(r1);

    //assert(0);
#endif
    
    for (; i<words; i++) {
        uint64_t prod;
        uint32_t prodl, prodh;

        prod = (uint64_t)a[i]*b;
        prodl = (uint32_t)prod;
        prodh = (uint32_t)(prod >> 32);

        prodl += carry; prodh += prodl < carry;
        t[i] += prodl; prodh += t[i] < prodl;
        carry = prodh;
    }

    for (;carry; i++) {
        t[i] += carry; carry = t[i] < carry;
    }
}


uint64_t addmul128(uint64_t * RESTRICT t, const uint64_t * RESTRICT a, uint64_t b0, uint64_t b1, size_t words)
{
#if 1
    uint32_t b0l, b0h, b1l, b1h;

    b0l = (uint32_t)b0;
    b0h = (uint32_t)(b0 >> 32);
    b1l = (uint32_t)b1;
    b1h = (uint32_t)(b1 >> 32);

    addmul32((uint32_t*)t+0, (uint32_t*)a, b0l, 2*words);
    addmul32((uint32_t*)t+1, (uint32_t*)a, b0h, 2*words);
    addmul32((uint32_t*)t+2, (uint32_t*)a, b1l, 2*words);
    addmul32((uint32_t*)t+3, (uint32_t*)a, b1h, 2*words);
#else
    addmul64(t+0, a, b0, words);
    addmul64(t+1, a, b1, words);
#endif
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

