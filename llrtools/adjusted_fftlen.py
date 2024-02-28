#!/usr/bin/python3
"""Script to print out max n values for FFT lengths adjusted for a given k."""
# MIT License
#
# Copyright (c) 2018-2020 Alexander Jones
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
import sys

from fftlen import FFTLengthK


max_data = []


def load_max_file() -> None:
    """Load the max length data file."""
    with open('maxlen.txt', 'rt') as max_file:
        for line in max_file:
            nums = line.split()
            max_data.append(tuple(map(int, nums)))


def nmax(fftlen: int, mersenne_n: int, fftlen_k: FFTLengthK) -> float:
    """Calculate the adjusted maximum n for the given k and FFT length."""
    if fftlen_k.k < 2**20:
        return fftlen_k.n_max(fftlen, mersenne_n)
    else:
        return fftlen_k.n_max_zero_padded(fftlen, mersenne_n)


def list_fftlens(fftlen_k: FFTLengthK) -> None:
    """List FFT lengths for the given k."""
    print("FFT lengths for \033[3mk\033[m={}".format(fftlen_k.k))
    for fftlen, mersenne_n in max_data:
        print("{:>8} {:>9}".format(fftlen, int(nmax(fftlen, mersenne_n, fftlen_k))))


def main() -> None:
    k = int(sys.argv[1])
    fftlen_k = FFTLengthK(k)
    load_max_file()
    list_fftlens(fftlen_k)


if __name__ == '__main__':
    main()
