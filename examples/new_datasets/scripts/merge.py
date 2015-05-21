#!/usr/bin/env python

"""
creates a ml20.base and ml20.test file.
run split.py first. 
for each user the first 20% of their ratings are moved to ml20.test
the rest is written to the base file
"""

splitted_files = 'split_per_user/'
import sys

if len(sys.argv) < 2:
    print 'USAGE: python merge.py OUTPUT_FILE_NAME'
    sys.exit(1)

from os import listdir
from os.path import isfile, join
onlyfiles = [ f for f in listdir(splitted_files) if isfile(join(splitted_files,f)) ]
files = sorted(onlyfiles, key=lambda x: int(x.split('_')[1].split('.')[0]))

base_file = file(sys.argv[1] + '.base', 'w')
test_file = file(sys.argv[1] + '.test', 'w')

for fname in files:
    with open(join(splitted_files, fname), 'r') as f:
        lines = f.readlines()
        num_test_lines = (len(lines) -1) / 5 # 20 percent tests
        for l in lines[1:num_test_lines+1]:
            test_file.write(l)
        for l in lines[num_test_lines+1:]:
            base_file.write(l)

base_file.close()
test_file.close()

