#!/bin/bash

#
# measures the runtime for a list of commits. will create REPO/results directory with a textfile for each commit and 
# input-file.
#

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
COMMITS=$(cat $DIR/old_measurements.txt  | grep '^commit' | sed -e 's/commit //')

BUILD_DIR=$DIR/../../build/
REPO_BASE=$DIR/../../
# result and testfiles need to be copied somewhere else. going back in git history will remove them from $DIR/testfiles
RESULT_DIR=$REPO_BASE/results/
TEST_FILES=$REPO_BASE/input_files/
rm -rf $RESULT_DIR
rm -rf $TEST_FILES
mkdir -p $RESULT_DIR
mkdir -p $TEST_FILES
cp $DIR/testfiles/* $TEST_FILES

for c in $COMMITS; do
    # empty the build dir
    rm -rf $BUILD_DIR
    mkdir $BUILD_DIR
    cd $BUILD_DIR
    
    git checkout $c

    cmake -DCMAKE_BUILD_TYPE=Release -DWITH_PERF=TRUE -DWITH_RECOMMENDER_BOOST=TRUE ..
    make -j8
    cd examples/

    cp $TEST_FILES/* .
    ./example_recommendation -dataset 32MoviesUID1 -numMeasurements 7 > $RESULT_DIR/${c}_32MoviesUID1.txt
    ./example_recommendation -dataset 512MoviesUID1 -numMeasurements 7 > $RESULT_DIR/${c}_512MoviesUID1.txt
    ./example_recommendation -dataset ml100kfull -numMeasurements 7 > $RESULT_DIR/${c}_ml100kfull.txt
done
