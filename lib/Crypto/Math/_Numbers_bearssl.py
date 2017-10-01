from Crypto.Util._raw_api import (load_pycryptodome_raw_lib,
                                  create_string_buffer,
                                  get_raw_buffer)

from Crypto.Util.number import long_to_bytes, bytes_to_long

# For now, only support CFFI, not ctypes
from Crypto.Util._raw_api import ffi

_bearssl_defs = """
uint32_t br_i31_iszero(const uint32_t *x);
uint32_t br_i31_add(uint32_t *a, const uint32_t *b, uint32_t ctl);
uint32_t br_i31_sub(uint32_t *a, const uint32_t *b, uint32_t ctl);
uint32_t br_i31_bit_length(uint32_t *x, size_t xlen);
void br_i31_decode(uint32_t *x, const void *src, size_t len);
uint32_t br_i31_decode_mod(uint32_t *x,
	const void *src, size_t len, const uint32_t *m);
void br_i31_rshift(uint32_t *x, int count);
void br_i31_reduce(uint32_t *x, const uint32_t *a, const uint32_t *m);
void br_i31_decode_reduce(uint32_t *x,
	const void *src, size_t len, const uint32_t *m);
void br_i31_muladd_small(uint32_t *x, uint32_t z, const uint32_t *m);
void br_i31_encode(void *dst, size_t len, const uint32_t *x);
uint32_t br_i31_ninv31(uint32_t x);
void br_i31_montymul(uint32_t *d, const uint32_t *x, const uint32_t *y,
	const uint32_t *m, uint32_t m0i);
void br_i31_to_monty(uint32_t *x, const uint32_t *m);
void br_i31_from_monty(uint32_t *x, const uint32_t *m, uint32_t m0i);
void br_i31_modpow(uint32_t *x, const unsigned char *e, size_t elen,
	const uint32_t *m, uint32_t m0i, uint32_t *t1, uint32_t *t2);
uint32_t br_i31_modpow_opt(uint32_t *x, const unsigned char *e, size_t elen,
	const uint32_t *m, uint32_t m0i, uint32_t *tmp, size_t twlen);
void br_i31_mulacc(uint32_t *d, const uint32_t *a, const uint32_t *b);
"""

_raw_bearssl = load_pycryptodome_raw_lib("Crypto.Util._bearssl", _bearssl_defs)

# With i31, the number is represented with 31-bit limbs via an array of 32-bit words.
#
#
#
# The first word of the array is the "announced bit length" (i.e. the # of meaningful bits):
# - the upper 27 bits are the numer of full limbs
# - the lower 5 bits are the residual bits in the last limb - if 0 there is no last limb
#
# Example of total array size for:
# - a 0 bit number:     1 word
# - a 1 bit number:     2 words
# - a 32 bit number:    2 words
# - a 33 bit number:    3 words
#
# For an array X, the idiom (X[0] + 31) >> 5 is used to compute the number of limbs.
# The total array length is clearly (X[0] + 31) >> 5 + 1.

_MAX_BITSIZE = 4096
_MAX_SIZE = 1 << _MAX_BITSIZE - 1
_MAX_NR_WORDS = _MAX_BITSIZE // 31 + 2
_MAX_BYTESIZE = _MAX_BITSIZE // 8 + 1

class Integer(object):

    def __init__(self, value):
        """Initialize the integer to the given value."""

        if isinstance(value, float):
            raise ValueError("A floating point type is not a natural number")

        if value < 0 or value > _MAX_SIZE:
            raise ValueError("Integer not in valid range")

        self._i31 = ffi.new("uint32_t[]", _MAX_NR_WORDS)

        buf = long_to_bytes(value)
        _raw_bearssl.br_i31_decode(self._i31, buf, len(buf))

    # Conversions
    def __int__(self):
        return bytes_to_long(self.to_bytes())

    def __str__(self):
        return str(int(self))

    def __repr__(self):
        return "Integer(%s)" % str(self)

    def to_bytes(self, block_size=0):
        """Convert the number into a byte string.

        This method encodes the number in network order and prepends
        as many zero bytes as required. It only works for non-negative
        values.

        :Parameters:
          block_size : integer
            The exact size the output byte string must have.
            If zero, the string has the minimal length.
        :Returns:
          A byte string.
        :Raise ValueError:
          If the value is negative or if ``block_size`` is
          provided and the length of the byte string would exceed it.
        """

        buf = create_string_buffer(self.size_in_bytes())
        _raw_bearssl.br_i31_encode(buf, len(buf), self._i31)
        return get_raw_buffer(buf)

    @classmethod
    def from_bytes(cls, byte_string):
        """Convert a byte string into a number.

        :Parameters:
          byte_string : byte string
            The input number, encoded in network order.
            It can only be non-negative.
        :Return:
          The ``Integer`` object carrying the same value as the input.
        """

        self = cls.__new__(cls)
        self._i31 = ffi.new("uint32_t[]", _MAX_NR_WORDS)
        _raw_bearssl.br_i31_decode(self._i31, byte_string, len(byte_string))
        return self

    # Relations
    def _apply_and_return(self, func, term):
        if not isinstance(term, Integer):
            term = Integer(term)
        return func(self._mpz_p, term._mpz_p)

    def __eq__(self, term):
        return int(self) == int(term)

    def __ne__(self, term):
        return int(self) != int(term)

    def __lt__(self, term):
        return int(self) < int(term)

    def __le__(self, term):
        return int(self) <= int(term)

    def __gt__(self, term):
        return int(self) > int(term)

    def __ge__(self, term):
        return int(self) >= int(term)

    def __nonzero__(self):
        return _raw_bearssl.br_i31_iszero(self._i31) == 0

    def is_negative(self):
        raise NotImplementedError("Not implemented yet")

    # Arithmetic operations
    def __add__(self, term):
        result = cls.__new__(cls)
        result._i31 = ffi.new("uint32_t[]", _MAX_NR_WORDS)
        ffi.memmove(result._i31, self._i31, MAX_NR_WORDS * 4)
        _raw_bearssl.br_i31_add(result._i31, term._i31, 1);
        return result

    def __sub__(self, term):
        raise NotImplementedError("Not implemented yet")

    def __mul__(self, term):
        raise NotImplementedError("Not implemented yet")

    def __floordiv__(self, divisor):
        raise NotImplementedError("Not implemented yet")

    def __mod__(self, divisor):
        raise NotImplementedError("Not implemented yet")

    def inplace_pow(self, exponent, modulus=None):
        raise NotImplementedError("Not implemented yet")

    def __pow__(self, exponent, modulus=None):
        raise NotImplementedError("Not implemented yet")

    def __abs__(self):
        raise NotImplementedError("Not implemented yet")

    def sqrt(self):
        """Return the largest Integer that does not
        exceed the square root"""

        raise NotImplementedError("Not implemented yet")

    def __iadd__(self, term):
        raise NotImplementedError("Not implemented yet")

    def __isub__(self, term):
        raise NotImplementedError("Not implemented yet")

    def __imul__(self, term):
        raise NotImplementedError("Not implemented yet")

    def __imod__(self, divisor):
        raise NotImplementedError("Not implemented yet")

    # Boolean/bit operations
    def __and__(self, term):
        raise NotImplementedError("Not implemented yet")

    def __or__(self, term):
        raise NotImplementedError("Not implemented yet")

    def __rshift__(self, pos):
        result = cls.__new__(cls)
        result._i31 = ffi.new("uint32_t[]", _MAX_NR_WORDS)
        ffi.memmove(result._i31, self._i31, MAX_NR_WORDS * 4)
        result >>= pos
        return result

    def __irshift__(self, pos):
        _raw_bearssl.br_i31_rshift(self._i31, pos)
        return self

    def __lshift__(self, pos):
        raise NotImplementedError("Not implemented yet")

    def __ilshift__(self, pos):
        raise NotImplementedError("Not implemented yet")

    def get_bit(self, n):
        """Return True if the n-th bit is set to 1.
        Bit 0 is the least significant."""

        raise NotImplementedError("Not implemented yet")

    # Extra
    def is_odd(self):
        raise NotImplementedError("Not implemented yet")

    def is_even(self):
        raise NotImplementedError("Not implemented yet")

    def size_in_bits(self):
        """Return the minimum number of bits that can encode the number."""

        if self._i31[0] == 0:
            return 1
        limbs = (self._i31[0] >> 5) + 1
        res = _raw_bearssl.br_i31_bit_length(self.i31 + 1, limbs)
        return res

    def size_in_bytes(self):
        """Return the minimum number of bytes that can encode the number."""
        return (self.size_in_bits() - 1) // 8 + 1

    def is_perfect_square(self):
        raise NotImplementedError("Not implemented yet")

    def fail_if_divisible_by(self, small_prime):
        """Raise an exception if the small prime is a divisor."""
        raise NotImplementedError("Not implemented yet")

    def multiply_accumulate(self, a, b):
        """Increment the number by the product of a and b."""
        raise NotImplementedError("Not implemented yet")

    def set(self, source):
        """Set the Integer to have the given value"""
        buf = long_to_bytes(value)
        _raw_bearssl.br_i31_decode(self._i31, buf, len(buf))

    def inplace_inverse(self, modulus):
        """Compute the inverse of this number in the ring of
        modulo integers.

        Raise an exception if no inverse exists.
        """

        raise NotImplementedError("Not implemented yet")

    def inverse(self, modulus):
        raise NotImplementedError("Not implemented yet")

    def gcd(self, term):
        """Compute the greatest common denominator between this
        number and another term."""

        raise NotImplementedError("Not implemented yet")

    def lcm(self, term):
        """Compute the least common multiplier between this
        number and another term."""

        raise NotImplementedError("Not implemented yet")

    @staticmethod
    def jacobi_symbol(a, n):
        """Compute the Jacobi symbol"""
        raise NotImplementedError("Not implemented yet")

    # Clean-up
    def __del__(self):

        try:
            self._i31 = None
        except AttributeError:
            pass
