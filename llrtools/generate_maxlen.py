#!/usr/bin/python3
"""Script to generate maxlen.txt and LLR test input for use with LLRTools."""
# MIT License
#
# Copyright (c) 2020 Alexander Jones
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

import argparse
import math
import re
import subprocess
import sys
import tempfile
from typing import Optional, Tuple

import fftlen as fftlen_lib

# Globals
VERBOSE = False

temp_dir_name = ''


# Processing code

fftlen_re = re.compile(rb'FFT length (\d+[kmKM]?)')


def parse_formatted_fftlen(fftlen: str) -> int:
    """Parse a formatted FFT length."""
    if fftlen.endswith(('m', 'M')):
        mega = int(fftlen[:-1])
        return mega * 2**20
    elif fftlen.endswith(('k', 'K')):
        kilo = int(fftlen[:-1])
        return kilo * 2**10
    else:
        return int(fftlen)


def get_fftlen_from_test(k: int, n: int) -> Optional[int]:
    """Get an FFT length from parsing the output of an LLR test."""
    try:
        proc = subprocess.run(['./llr', f'-w{temp_dir_name}', '-d', f'-q{k}*2^{n}-1'], capture_output=True, timeout=2)
        output = proc.stdout
    except subprocess.TimeoutExpired as ex:
        output = ex.output

    match = fftlen_re.search(output)
    if match:
        return parse_formatted_fftlen(str(match[1], 'utf-8'))
    else:
        return None


# Adapted from https://en.wikipedia.org/wiki/Binary_search_algorithm#Procedure_for_finding_the_leftmost_element.
def binary_search(k: int, fftlen: int, start: int, finish: int) -> Tuple[int, int]:
    """Find the smallest n with the next FFT length."""
    left = start
    right = finish
    old_n = 0
    while left < right:
        n = int((left + right) / 2)
        n_fft = get_fftlen_from_test(k, n)
        while n_fft is None:
            if old_n == n:
                return n, right
            else:
                old_n = n
            n += 1
            n_fft = get_fftlen_from_test(k, n)
        if VERBOSE:
            print(f'n={n}, FFT={n_fft}')
        if n_fft < fftlen:
            left = n + 1
        else:
            right = n
    return left, left


def loop(k: int, fftlen_max: Optional[int], n_max: Optional[int], n_min: Optional[int]) -> None:
    """Loop until the maximum FFT length is reached."""
    with open('testinput.new.txt', 'wt') as testinput, open('maxlen.new.txt', 'wt') as maxlen:
        print('1000000000000:M:1:2:258', file=testinput)
        start_n = n_min or fftlen_lib.n_max(32, 700)
        last_fftlen = get_fftlen_from_test(k, n_min) if n_min else 32
        next_n = start_n
        next_fftlen = last_fftlen
        while (fftlen_max and last_fftlen <= fftlen_max) or (n_max and start_n <= n_max):
            while next_fftlen == last_fftlen:
                next_n = int(next_n * 1.02)
                next_fftlen = get_fftlen_from_test(k, next_n) or last_fftlen
                if VERBOSE:
                    print(f'n={next_n}, FFT={next_fftlen}')
            left_n, test_n = binary_search(k, next_fftlen, start_n, next_n)
            if last_fftlen == 32:
                first_n = test_n - 1
                while get_fftlen_from_test(k, first_n) is None:
                    first_n -= 1
                print(f'{k} {first_n}', file=testinput)
            print(f'{k} {test_n}', file=testinput)
            mersenne_left_n = int(fftlen_lib.mersenne(last_fftlen, left_n))
            print(f'{last_fftlen:>8} {mersenne_left_n:>9}', file=maxlen)
            print(f'FFT={last_fftlen} done.')
            last_fftlen = next_fftlen
            start_n = next_n


def start() -> None:
    """Main starting function."""
    parser = argparse.ArgumentParser(description='Generate maxlen.txt and LLR test input for use with LLRTools.')
    parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
    parser.add_argument('-k', type=int, default=100005, help='set testing k')
    parser.add_argument('-m', type=int, help='set minimum n')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-n', type=int, help='set maximum n')
    group.add_argument('-f', '--fftlen', type=str, help='set maximum FFT length')
    args = parser.parse_args()
    global VERBOSE
    VERBOSE = args.verbose
    fftlen_max = parse_formatted_fftlen(args.fftlen) if args.fftlen else None
    n_max = args.n
    n_min = args.m
    k = args.k
    fftlen_lib.change_k(k)
    with tempfile.TemporaryDirectory() as tmpdirname:
        global temp_dir_name
        temp_dir_name = tmpdirname
        loop(k, fftlen_max, n_max, n_min)


if __name__ == '__main__':
    start()
