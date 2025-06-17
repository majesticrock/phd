#!/bin/bash

rm -rf build/
mkdir build

cmake . -B build/

grep "INCLUDE_DIR" build/CMakeCache.txt

