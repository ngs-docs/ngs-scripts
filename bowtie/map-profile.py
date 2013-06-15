#! /usr/bin/env python
import csv
import sys

map_file = sys.argv[1]

mismatch_count = [0] * 200

r = csv.reader(open(map_file), delimiter='\t')
for n, row in enumerate(r):
    if n % 100000 == 0:
        print >>sys.stderr, '...', n

    mismatches = row[7]

    mismatch_list = mismatches.split(',')
    for mismatch in mismatch_list:
        if not mismatch:
            continue
        mismatch = mismatch.split(':')
        mismatch = mismatch[0]
        mismatch = int(mismatch)
        mismatch_count[mismatch] += 1

for position, count in enumerate(mismatch_count):
    print position, count
