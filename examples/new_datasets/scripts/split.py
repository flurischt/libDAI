#!/usr/bin/env python
import os
import sys

"""
preprocess step for merge.py
"""

if len(sys.argv) < 2:
    print 'USAGE: python split.py INPUT_FILE_NAME'
    sys.exit(1)

in_file = file(sys.argv[1], 'r')

header = None
last_user = 0
out_file = None

try:
    os.mkdir('split_per_user')
except OSError:
    pass

# split the ratings.csv file into a file for each user
for idx, line in enumerate(in_file.readlines()):
    if idx == 0:
        header = line
        continue
    current_user = int(line.split('\t')[0])
    if last_user != current_user:
        if last_user != 0:
            out_file.close()
        out_file = file('split_per_user/user_%d.csv' % current_user, 'w')
        out_file.write(header)
        last_user = current_user
    out_file.write(line)
out_file.close()

