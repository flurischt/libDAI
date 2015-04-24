#!/bin/bash
# jenkins uses this script to build and test
BUILD_DIR=./testBuild/

rm -rf $BUILD_DIR
mkdir $BUILD_DIR
cd $BUILD_DIR
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
cd examples/
./example_recommendation -dataset u1 -test -cpufreq 2300000000 > output.txt
cat output.txt
