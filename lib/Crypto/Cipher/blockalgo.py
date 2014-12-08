# -*- coding: utf-8 -*-
#
#  Cipher/blockalgo.py
#
# ===================================================================
# The contents of this file are dedicated to the public domain.  To
# the extent that dedication to the public domain is not available,
# everyone is granted a worldwide, perpetual, royalty-free,
# non-exclusive license to exercise all rights associated with the
# contents of this file for any purpose whatsoever.
# No rights are reserved.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ===================================================================
"""Module with definitions common to all block ciphers."""

from Crypto.Util.py3compat import *

from binascii import unhexlify

from Crypto.Util import Counter
from Crypto.Util.strxor import strxor
from Crypto.Util.number import long_to_bytes, bytes_to_long
import Crypto.Util.Counter
from Crypto.Hash import CMAC
from Crypto.Hash.CMAC import _SmoothMAC

from Crypto.Util import _galois

#: *Electronic Code Book (ECB)*.
#: This is the simplest encryption mode. Each of the plaintext blocks
#: is directly encrypted into a ciphertext block, independently of
#: any other block. This mode exposes frequency of symbols
#: in your plaintext. Other modes (e.g. *CBC*) should be used instead.
#:
#: See `NIST SP800-38A`_ , Section 6.1 .
#:
#: .. _`NIST SP800-38A` : http://csrc.nist.gov/publications/nistpubs/800-38a/sp800-38a.pdf
MODE_ECB = 1

#: *Cipher-Block Chaining (CBC)*. Each of the ciphertext blocks depends
#: on the current and all previous plaintext blocks. An Initialization Vector
#: (*IV*) is required.
#:
#: The *IV* is a data block to be transmitted to the receiver.
#: The *IV* can be made public, but it must be authenticated by the receiver
#: and it should be picked randomly.
#:
#: See `NIST SP800-38A`_ , Section 6.2 .
#:
#: .. _`NIST SP800-38A` : http://csrc.nist.gov/publications/nistpubs/800-38a/sp800-38a.pdf
MODE_CBC = 2

#: *Cipher FeedBack (CFB)*. This mode is similar to CBC, but it transforms
#: the underlying block cipher into a stream cipher. Plaintext and ciphertext
#: are processed in *segments* of **s** bits. The mode is therefore sometimes
#: labelled **s**-bit CFB. An Initialization Vector (*IV*) is required.
#:
#: When encrypting, each ciphertext segment contributes to the encryption of
#: the next plaintext segment.
#:
#: This *IV* is a data block to be transmitted to the receiver.
#: The *IV* can be made public, but it should be picked randomly.
#: Reusing the same *IV* for encryptions done with the same key lead to
#: catastrophic cryptographic failures.
#:
#: See `NIST SP800-38A`_ , Section 6.3 .
#:
#: .. _`NIST SP800-38A` : http://csrc.nist.gov/publications/nistpubs/800-38a/sp800-38a.pdf
MODE_CFB = 3

#: This mode should not be used.
MODE_PGP = 4

#: *Output FeedBack (OFB)*. This mode is very similar to CBC, but it
#: transforms the underlying block cipher into a stream cipher.
#: The keystream is the iterated block encryption of an
#: Initialization Vector (*IV*).
#:
#: The *IV* is a data block to be transmitted to the receiver.
#: The *IV* can be made public, but it should be picked randomly.
#:
#: Reusing the same *IV* for encryptions done with the same key lead to
#: catastrophic cryptograhic failures.
#:
#: See `NIST SP800-38A`_ , Section 6.4 .
#:
#: .. _`NIST SP800-38A` : http://csrc.nist.gov/publications/nistpubs/800-38a/sp800-38a.pdf
MODE_OFB = 5

#: *CounTeR (CTR)*. This mode is very similar to ECB, in that
#: encryption of one block is done independently of all other blocks.
#: Unlike ECB, the block *position* contributes to the encryption and no
#: information leaks about symbol frequency.
#:
#: Each message block is associated to a *counter* which must be unique
#: across all messages that get encrypted with the same key (not just within
#: the same message). The counter is as big as the block size.
#:
#: Counters can be generated in several ways. The most straightword one is
#: to choose an *initial counter block* (which can be made public, similarly
#: to the *IV* for the other modes) and increment its lowest **m** bits by
#: one (modulo *2^m*) for each block. In most cases, **m** is chosen to be half
#: the block size.
#:
#: Reusing the same *initial counter block* for encryptions done with the same
#: key lead to catastrophic cryptograhic failures.
#:
#: See `NIST SP800-38A`_ , Section 6.5 (for the mode) and Appendix B (for how
#: to manage the *initial counter block*).
#:
#: .. _`NIST SP800-38A` : http://csrc.nist.gov/publications/nistpubs/800-38a/sp800-38a.pdf
MODE_CTR = 6

#: *EAX*. This is an Authenticated Encryption with Associated Data
#: (`AEAD`_) mode. It provides both confidentiality and authenticity.
#:
#: The header of the message may be left in the clear, if needed, and it will
#: still be subject to authentication.
#:
#: The decryption step tells the receiver if the message comes from a source
#: that really knowns the secret key.
#: Additionally, decryption detects if any part of the message - including the
#: header - has been modified or corrupted.
#:
#: This mode requires a nonce. The nonce shall never repeat for two
#: different messages encrypted with the same key, but it does not need to
#: be random.
#
#: This mode is only available for ciphers that operate on 64 or
#: 128 bits blocks.
#:
#: There are no official standards defining EAX. The implementation is based on
#: `a proposal`__ that was presented to NIST.
#:
#: .. _AEAD: http://blog.cryptographyengineering.com/2012/05/how-to-choose-authenticated-encryption.html
#: .. __: http://csrc.nist.gov/groups/ST/toolkit/BCM/documents/proposedmodes/eax/eax-spec.pdf
MODE_EAX = 9

#: *Galois/Counter Mode (GCM)*. This is an Authenticated Encryption with
#: Associated Data (`AEAD`_) mode. It provides both confidentiality and
#: authenticity.
#: The header of the message may be left in the clear, if needed, and it will
#: still be subject to authentication. The decryption step tells the receiver
#: if the message comes from a source that really knowns the secret key.
#: Additionally, decryption detects if any part of the message - including the
#: header - has been modified or corrupted.
#:
#: This mode requires a nonce. The nonce shall never repeat for two
#: different messages encrypted with the same key, but it does not need to
#: be random.
#:
#: This mode is only available for ciphers that operate on 128 bits blocks
#: (e.g. AES but not TDES).
#:
#: See `NIST SP800-38D`_ .
#:
#: .. _`NIST SP800-38D`: http://csrc.nist.gov/publications/nistpubs/800-38D/SP-800-38D.pdf
#: .. _AEAD: http://blog.cryptographyengineering.com/2012/05/how-to-choose-authenticated-encryption.html
MODE_GCM = 11


def _getParameter(name, index, args, kwargs, default=None):
    """Find a parameter in tuple and dictionary arguments a function receives"""

    param = kwargs.get(name)
    if len(args) > index:
        if param:
            raise TypeError("Parameter '%s' is specified twice" % name)
        param = args[index]
    return param or default


class _GHASH(_SmoothMAC):
    """GHASH function defined in NIST SP 800-38D, Algorithm 2.

    If X_1, X_2, .. X_m are the blocks of input data, the function
    computes:

       X_1*H^{m} + X_2*H^{m-1} + ... + X_m*H

    in the Galois field GF(2^256) using the reducing polynomial
    (x^128 + x^7 + x^2 + x + 1).
    """

    def __init__(self, hash_subkey, block_size):
        _SmoothMAC.__init__(self, block_size, None, 0)
        self._hash_subkey = _galois.ghash_expand(hash_subkey)
        self._last_y = bchr(0) * 16
        self._mac = _galois.ghash

    def copy(self):
        clone = _GHASH(self._hash_subkey, self._bs)
        _SmoothMAC._deep_copy(self, clone)
        clone._last_y = self._last_y
        return clone

    def _update(self, block_data):
        self._last_y = _galois.ghash(block_data, self._last_y,
                                     self._hash_subkey)

    def _digest(self, left_data):
        return self._last_y


class BlockAlgo:
    """Class modelling an abstract block cipher."""

    def __init__(self, factory, key, *args, **kwargs):
        self.mode = _getParameter('mode', 0, args, kwargs, default=MODE_ECB)
        self.block_size = factory.block_size
        self._factory = factory
        self._tag = None

        if self.mode == MODE_EAX:
            self._start_eax(factory, key, *args, **kwargs)
        elif self.mode == MODE_GCM:
            self._start_gcm(factory, key, *args, **kwargs)
        else:
            self._cipher = factory.new(key, *args, **kwargs)
            self.IV = self._cipher.IV

    def _start_gcm(self, factory, key, *args, **kwargs):

        if self.block_size != 16:
            raise TypeError("GCM mode is only available for ciphers that operate on 128 bits blocks")

        self.nonce = _getParameter('nonce', 1, args, kwargs)
        if not self.nonce:
            raise TypeError("MODE_GCM requires a nonce")

        self._mac_len = kwargs.get('mac_len', 16)
        if not (self._mac_len and 4 <= self._mac_len <= 16):
            raise ValueError("Parameter 'mac_len' must not be larger than 16 bytes")

        # Allowed transitions after initialization
        self._next = [self.update, self.encrypt, self.decrypt,
                      self.digest, self.verify]

        self._done_assoc_data = False

        # Length of the ciphertext or plaintext
        self._msg_len = 0

        # Step 1 in SP800-38D, Algorithm 4 (encryption) - Compute H
        # See also Algorithm 5 (decryption)
        hash_subkey = factory.new(key, MODE_ECB).encrypt(bchr(0) * 16)

        # Step 2 - Compute J0 (integer, not byte string!)
        if len(self.nonce) == 12:
            self._j0 = bytes_to_long(self.nonce + b("\x00\x00\x00\x01"))
        else:
            fill = (16 - (len(self.nonce) % 16)) % 16 + 8
            ghash_in = (self.nonce +
                        bchr(0) * fill +
                        long_to_bytes(8 * len(self.nonce), 8))

            mac = _GHASH(hash_subkey, factory.block_size)
            mac.update(ghash_in)
            self._j0 = bytes_to_long(mac.digest())

        # Step 3 - Prepare GCTR cipher for encryption/decryption
        ctr = Counter.new(128, initial_value=self._j0 + 1,
                          allow_wraparound=True)
        self._cipher = self._factory.new(key, MODE_CTR, counter=ctr)

        # Step 5 - Bootstrat GHASH
        self._cipherMAC = _GHASH(hash_subkey, factory.block_size)

        # Step 6 - Prepare GCTR cipher for GMAC
        ctr = Counter.new(128, initial_value=self._j0, allow_wraparound=True)
        self._tag_cipher = self._factory.new(key, MODE_CTR, counter=ctr)

    def _start_eax(self, factory, key, *args, **kwargs):

        self.nonce = _getParameter('nonce', 1, args, kwargs)
        if not self.nonce:
            raise TypeError("MODE_EAX requires a nonce")

        # Allowed transitions after initialization
        self._next = [self.update, self.encrypt, self.decrypt,
                      self.digest, self.verify]

        self._mac_len = kwargs.get('mac_len', self.block_size)
        if not (self._mac_len and 4 <= self._mac_len <= self.block_size):
            raise ValueError("Parameter 'mac_len' must not be larger than %d"
                             % self.block_size)

        self._omac = [
                CMAC.new(key, bchr(0) * (self.block_size - 1) + bchr(i),
                         ciphermod=factory)
                for i in xrange(0, 3)
                ]

        # Compute MAC of nonce
        self._omac[0].update(self.nonce)

        self._cipherMAC = self._omac[1]

        # MAC of the nonce is also the initial counter for CTR encryption
        counter_int = bytes_to_long(self._omac[0].digest())
        counter_obj = Crypto.Util.Counter.new(
                        self.block_size * 8,
                        initial_value=counter_int,
                        allow_wraparound=True)
        self._cipher = factory.new(key, MODE_CTR, counter=counter_obj)

    def update(self, assoc_data):
        """Protect associated data

        When using an AEAD mode like CCM, EAX or GCM and
        if there is any associated data, the caller has to invoke
        this function one or more times, before using
        ``decrypt`` or ``encrypt``.

        By *associated data* it is meant any data (e.g. packet headers) that
        will not be encrypted and will be transmitted in the clear.
        However, the receiver is still able to detect any modification to it.
        In CCM and GCM, the *associated data* is also called
        *additional authenticated data* (AAD).
        In EAX, the *associated data* is called *header*.

        If there is no associated data, this method must not be called.

        The caller may split associated data in segments of any size, and
        invoke this method multiple times, each time with the next segment.

        :Parameters:
          assoc_data : byte string
            A piece of associated data. There are no restrictions on its size.
        """

        if self.mode not in (MODE_EAX, MODE_GCM):
            raise TypeError("update() not supported by this mode of operation")

        if self.update not in self._next:
            raise TypeError("update() can only be called immediately after initialization")

        self._next = [self.update, self.encrypt, self.decrypt,
                      self.digest, self.verify]

        return self._cipherMAC.update(assoc_data)

    def encrypt(self, plaintext):
        """Encrypt data with the key and the parameters set at initialization.

        A cipher object is stateful: once you have encrypted a message
        you cannot encrypt (or decrypt) another message using the same
        object.

        For all other modes, the data to encrypt can be broken up in two or
        more pieces and `encrypt` can be called multiple times.

        That is, the statement:

            >>> c.encrypt(a) + c.encrypt(b)

        is equivalent to:

             >>> c.encrypt(a+b)

        That also means that you cannot reuse an object for encrypting
        or decrypting other data with the same key.

        This function does not add any padding to the plaintext.

         - For `MODE_ECB` and `MODE_CBC`, *plaintext* length (in bytes) must be
           a multiple of *block_size*.

         - For `MODE_CFB`, *plaintext* length (in bytes) must be a multiple
           of *segment_size*/8.

         - For `MODE_OFB`, `MODE_CTR` and all AEAD modes
           *plaintext* can be of any length.

        :Parameters:
          plaintext : byte string
            The piece of data to encrypt.
        :Return:
            the encrypted data, as a byte string. It is as long as
            *plaintext*.
        """

        if self.mode in (MODE_EAX, MODE_GCM):
            if self.encrypt not in self._next:
                raise TypeError("encrypt() can only be called after initialization or an update()")
            self._next = [self.encrypt, self.digest]

        ct = self._cipher.encrypt(plaintext)

        if self.mode == MODE_EAX:
            self._omac[2].update(ct)

        if self.mode == MODE_GCM:
            if not self._done_assoc_data:
                self._cipherMAC.zero_pad()
                self._done_assoc_data = True
            self._cipherMAC.update(ct)
            self._msg_len += len(plaintext)

        return ct

    def decrypt(self, ciphertext):
        """Decrypt data with the key and the parameters set at initialization.

        A cipher object is stateful: once you have decrypted a message
        you cannot decrypt (or encrypt) another message with the same
        object.

        For all other modes, the data to decrypt can be broken up in two or
        more pieces and `decrypt` can be called multiple times.

        That is, the statement:

            >>> c.decrypt(a) + c.decrypt(b)

        is equivalent to:

             >>> c.decrypt(a+b)

        That also means that you cannot reuse an object for encrypting
        or decrypting other data with the same key.

        This function does not remove any padding from the plaintext.

         - For `MODE_ECB` and `MODE_CBC`, *ciphertext* length (in bytes) must
           be a multiple of *block_size*.

         - For `MODE_CFB`, *ciphertext* length (in bytes) must be a multiple
           of *segment_size*/8.

         - For `MODE_OFB`, `MODE_CTR` and all AEAD modes
           *ciphertext* can be of any length.

        :Parameters:
          ciphertext : byte string
            The piece of data to decrypt.

        :Return: the decrypted data (byte string).
        """

        if self.mode in (MODE_EAX, MODE_GCM):

            if self.decrypt not in self._next:
                raise TypeError("decrypt() can only be called after initialization or an update()")
            self._next = [self.decrypt, self.verify]

            if self.mode == MODE_GCM:
                if not self._done_assoc_data:
                    self._cipherMAC.zero_pad()
                    self._done_assoc_data = True

                self._cipherMAC.update(ciphertext)
                self._msg_len += len(ciphertext)

            if self.mode == MODE_EAX:
                self._omac[2].update(ciphertext)

        pt = self._cipher.decrypt(ciphertext)

        return pt

    def digest(self):
        """Compute the *binary* MAC tag in an AEAD mode.

        When using an AEAD mode like CCM or EAX, the caller invokes
        this function at the very end.

        This method returns the MAC that shall be sent to the receiver,
        together with the ciphertext.

        :Return: the MAC, as a byte string.
        """

        if self.mode not in (MODE_EAX, MODE_GCM):
            raise TypeError("digest() not supported by this mode of operation")

        if self.digest not in self._next:
            raise TypeError("digest() cannot be called when decrypting or validating a message")
        self._next = [self.digest]

        return self._compute_mac()

    def _compute_mac(self):
        """Compute MAC without any FSM checks."""

        if self._tag:
            return self._tag

        if self.mode == MODE_GCM:

            # Step 5 in NIST SP 800-38D, Algorithm 4 - Compute S
            self._cipherMAC.zero_pad()
            auth_len = self._cipherMAC.data_signed_so_far() - self._msg_len
            for tlen in (auth_len, self._msg_len):
                self._cipherMAC.update(long_to_bytes(8 * tlen, 8))
            s_tag = self._cipherMAC.digest()

            # Step 6 - Compute T
            self._tag = self._tag_cipher.encrypt(s_tag)[:self._mac_len]

        if self.mode == MODE_EAX:
            tag = bchr(0) * self.block_size
            for i in xrange(3):
                tag = strxor(tag, self._omac[i].digest())
            self._tag = tag[:self._mac_len]

        return self._tag

    def hexdigest(self):
        """Compute the *printable* MAC tag in an AEAD mode.

        This method is like `digest`.

        :Return: the MAC, as a hexadecimal string.
        """
        return "".join(["%02x" % bord(x) for x in self.digest()])

    def verify(self, mac_tag):
        """Validate the *binary* MAC tag in an AEAD mode.

        When using an AEAD mode like CCM or EAX, the caller invokes
        this function at the very end.

        This method checks if the decrypted message is indeed valid
        (that is, if the key is correct) and it has not been
        tampered with while in transit.

        :Parameters:
          mac_tag : byte string
            This is the *binary* MAC, as received from the sender.
        :Raises ValueError:
            if the MAC does not match. The message has been tampered with
            or the key is incorrect.
        """

        if self.mode not in (MODE_EAX, MODE_GCM):
            raise TypeError("verify() not supported by this mode of operation")

        if self.verify not in self._next:
            raise TypeError("verify() cannot be called when encrypting a message")
        self._next = [self.verify]

        res = 0
        # Constant-time comparison
        for x, y in zip(self._compute_mac(), mac_tag):
            res |= bord(x) ^ bord(y)
        if res or len(mac_tag) != self._mac_len:
            raise ValueError("MAC check failed")

    def hexverify(self, hex_mac_tag):
        """Validate the *printable* MAC tag in an AEAD mode.

        This method is like `verify`.

        :Parameters:
          hex_mac_tag : string
            This is the *printable* MAC, as received from the sender.
        :Raises ValueError:
            if the MAC does not match. The message has been tampered with
            or the key is incorrect.
        """

        self.verify(unhexlify(hex_mac_tag))

    def encrypt_and_digest(self, plaintext):
        """Perform encrypt() and digest() in one step.

        :Parameters:
          plaintext : byte string
            The piece of data to encrypt.
        :Return:
            a tuple with two byte strings:

            - the encrypted data
            - the MAC
        """

        return self.encrypt(plaintext), self.digest()

    def decrypt_and_verify(self, ciphertext, mac_tag):
        """Perform decrypt() and verify() in one step.

        :Parameters:
          ciphertext : byte string
            The piece of data to decrypt.
          mac_tag : byte string
            This is the *binary* MAC, as received from the sender.

        :Return: the decrypted data (byte string).
        :Raises ValueError:
            if the MAC does not match. The message has been tampered with
            or the key is incorrect.
        """

        pt = self.decrypt(ciphertext)
        self.verify(mac_tag)
        return pt
