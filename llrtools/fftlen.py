"""FFT length adjustment routines."""
# MIT License
#
# Copyright (c) 2015, 2020, 2024 Alexander Jones
#
# Based on formulas from LLRTools by Thomas Ritschel.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Imports
import math


class FFTLengthK(object):
    def __init__(self, k: int) -> None:
        self.k = k

    def __int__(self) -> int:
        return self.k

    @property
    def log2k(self) -> float:
        return math.log(self.k, 2)

    def _adjust(self, fftlen: int) -> float:
        return (self.log2k + self.log2k * (fftlen / 2.0))

    def n_max(self, fftlen: int, mersenne: int) -> float:
        return mersenne - self._adjust(fftlen)

    def n_max_zero_padded(self, fftlen: int, mersenne: int) -> float:
        return (mersenne + 0.3*fftlen)/2.0

    def mersenne(self, fftlen: int, n_max: int) -> float:
        return n_max + self._adjust(fftlen)
