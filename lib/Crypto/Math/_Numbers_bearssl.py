from Crypto.Util._raw_api import (load_pycryptodome_raw_lib,
                                  create_string_buffer,
                                  get_raw_buffer)

from Crypto.Util.number import long_to_bytes, bytes_to_long

# For now, only support CFFI, not ctypes
from Crypto.Util._raw_api import ffi

_bearssl_defs = """
    void br_i31_decode(uint32_t *x, const void *src, size_t len);
    void br_i31_encode(void *dst, size_t len, const uint32_t *x);
"""

_raw_bearssl = load_pycryptodome_raw_lib("Crypto.Util._bearssl", _bearssl_defs)


class Integer(object):

    def __init__(self, value):
        """Initialize the integer to the given value."""

        if isinstance(value, float):
            raise ValueError("A floating point type is not a natural number")

        if value < 0 or value >= 2**4096:
            raise ValueError("Integer not in valid range")

        # We only support integer up to 4096 bits (133+1 31-bit words)
        self._i31 = ffi.new("uint32_t[]", 134)

        big_endian = long_to_bytes(value)
        _raw_bearssl.br_i31_decode(self._i31, big_endian, len(big_endian))

    # Conversions
    def __int__(self):

        # Enough room for 4096 bits
        # BearSSL will zero fill if too big
        buf = create_string_buffer(512)
        _raw_bearssl.br_i31_encode(buf, len(buf), self._i31)
        return bytes_to_long(get_raw_buffer(buf))

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

        raise NotImplementedError("Not implemented yet")

    @staticmethod
    def from_bytes(byte_string):
        """Convert a byte string into a number.

        :Parameters:
          byte_string : byte string
            The input number, encoded in network order.
            It can only be non-negative.
        :Return:
          The ``Integer`` object carrying the same value as the input.
        """

        raise NotImplementedError("Not implemented yet")

    # Relations
    def _apply_and_return(self, func, term):
        if not isinstance(term, Integer):
            term = Integer(term)
        return func(self._mpz_p, term._mpz_p)

    def __eq__(self, term):
        if not isinstance(term, (Integer, int, long)):
            return False
        raise NotImplementedError("Not implemented yet")

    def __ne__(self, term):
        if not isinstance(term, (Integer, int, long)):
            return True
        raise NotImplementedError("Not implemented yet")

    def __lt__(self, term):
        raise NotImplementedError("Not implemented yet")

    def __le__(self, term):
        raise NotImplementedError("Not implemented yet")

    def __gt__(self, term):
        raise NotImplementedError("Not implemented yet")

    def __ge__(self, term):
        raise NotImplementedError("Not implemented yet")

    def __nonzero__(self):
        raise NotImplementedError("Not implemented yet")

    def is_negative(self):
        raise NotImplementedError("Not implemented yet")

    # Arithmetic operations
    def __add__(self, term):
        raise NotImplementedError("Not implemented yet")

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
        raise NotImplementedError("Not implemented yet")

    def __irshift__(self, pos):
        raise NotImplementedError("Not implemented yet")

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
        raise NotImplementedError("Not implemented yet")

    def size_in_bytes(self):
        """Return the minimum number of bytes that can encode the number."""
        raise NotImplementedError("Not implemented yet")

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
        raise NotImplementedError("Not implemented yet")

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
