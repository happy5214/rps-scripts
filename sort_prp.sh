#!/bin/sh

head -n 1 "$1"
tail -q -n +2 "$@" | sort -t ' ' -k 1,1n -k 2,2n
