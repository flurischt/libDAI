#!/bin/bash
# jenkins uses this script to build and test
BUILD_DIR=./testBuild/
FLOAT_BUILD_DIR=./testBuildFloat/

# compile double
rm -rf $BUILD_DIR
mkdir $BUILD_DIR
cd $BUILD_DIR
cmake -DCMAKE_BUILD_TYPE=Release -DWITH_PERF=TRUE -DWITH_RECOMMENDER_BOOST=TRUE ..
make -j8 VERBOSE=1
cd examples/
cp ../../examples/new_datasets/ml100kfull.* .
cp ../../examples/new_datasets/u1/32MoviesUID1.* .
cp ../../examples/new_datasets/u1/512MoviesUID1.* .

# compile float
cd ../..
rm -rf $FLOAT_BUILD_DIR
mkdir $FLOAT_BUILD_DIR
cd $FLOAT_BUILD_DIR
cmake -DCMAKE_BUILD_TYPE=Release -DWITH_PERF=TRUE -DWITH_RECOMMENDER_BOOST=TRUE -DWITH_SINGLE_PRECISION=TRUE ..
make -j8 VERBOSE=1
cd examples/
cp ../../examples/new_datasets/ml100kfull.* .
cp ../../examples/new_datasets/u1/32MoviesUID1.* .
cp ../../examples/new_datasets/u1/512MoviesUID1.* .

# test double
cd ../../$BUILD_DIR/examples/
./example_recommendation -dataset u1 -numMeasurements 7 -test -cpufreq 2300000000 > /dev/null
./example_recommendation -dataset ml100kfull -numMeasurements 7 -cpufreq 2300000000 > output_ml100kfull.txt
cat output_ml100kfull.txt
./example_recommendation -dataset 32MoviesUID1 -numMeasurements 7 -cpufreq 2300000000 > output_32.txt
cat output_32.txt
./example_recommendation -dataset 512MoviesUID1 -numMeasurements 7 -cpufreq 2300000000 > output_512.txt
cat output_512.txt

# test float
cd ../../$FLOAT_BUILD_DIR/examples/
./example_recommendation -dataset ml100kfull -numMeasurements 7 -cpufreq 2300000000 > output_ml100kfull_float.txt
cat output_ml100kfull_float.txt
