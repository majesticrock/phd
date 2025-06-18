#!/bin/bash

rm -rf build/
mkdir build

cmake . -B build/
make -C build/

./build/test_boost_cmake

grep "INCLUDE_DIR" build/CMakeCache.txt

