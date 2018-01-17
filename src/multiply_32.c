#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "multiply.h"

void inline static addmul32(uint32_t* t, const uint32_t *a, uint32_t b, size_t words)
{
    uint32_t carry = 0;
    int i;

    for (i=0; i<words; i++) {
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
    uint32_t b0l, b0h, b1l, b1h;

    b0l = (uint32_t)b0;
    b0h = (uint32_t)(b0 >> 32);
    b1l = (uint32_t)b1;
    b1h = (uint32_t)(b1 >> 32);

    addmul32((uint32_t*)t+0, (uint32_t*)a, b0l, 2*words);
    addmul32((uint32_t*)t+1, (uint32_t*)a, b0h, 2*words);
    addmul32((uint32_t*)t+2, (uint32_t*)a, b1l, 2*words);
    addmul32((uint32_t*)t+3, (uint32_t*)a, b1h, 2*words);
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

