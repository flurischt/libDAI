#!/bin/bash

#
# merge the measurement results into three files for plotting
#

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
COMMITS=$(cat $DIR/old_measurements.txt  | grep '^commit' | sed -e 's/commit //')

BUILD_DIR=$DIR/../../build/
REPO_BASE=$DIR/../../
RESULT_DIR=$REPO_BASE/results/

for c in $COMMITS; do
    git log $c | head -n5 | grep -v '^$' >> measurement_32.txt
    git log $c | head -n5 | grep -v '^$' >> measurement_512.txt
    git log $c | head -n5 | grep -v '^$' >> measurement_ml100kfull.txt
    tail -n2 $RESULT_DIR/${c}_32MoviesUID1.txt >> measurement_32.txt
    tail -n2 $RESULT_DIR/${c}_512MoviesUID1.txt >> measurement_512.txt
    tail -n2 $RESULT_DIR/${c}_ml100kfull.txt >> measurement_ml100kfull.txt
    echo "" >> measurement_32.txt
    echo "" >> measurement_512.txt
    echo "" >> measurement_ml100kfull.txt
done
