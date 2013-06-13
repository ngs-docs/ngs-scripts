#! /usr/bin/env python
import sys
import blastparser

import screed

if len(sys.argv) != 5:
   print>>sys.stderr, "Usage: calc-blast-cover.py b.seqs a.x.b matchlen a.seqs"
   sys.exit(-1)

MIN_SCORE=200
MIN_QUERY_LEN = int(sys.argv[3])

query_seqs = set([ record.name for record in screed.open(sys.argv[4]) \
                       if len(record.sequence) >= MIN_QUERY_LEN ])

covs = {}
for n, record in enumerate(screed.open(sys.argv[1])):
    if n % 1000 == 0:
        sys.stdout.write('+')
        sys.stdout.flush()

    covs[record.name] = [0] * len(record.sequence)

for n, record in enumerate(blastparser.parse_fp(open(sys.argv[2]))):
    if n % 100 == 0:
        sys.stdout.write('.')
        sys.stdout.flush()

    if record.query_name not in query_seqs:
        continue

    for hit in record.hits:
        for match in hit.matches:
            if match.score < MIN_SCORE:
                continue

            cov = covs.get(hit.subject_name)
            if not cov:
                continue

            start = min(match.subject_start, match.subject_end) - 1
            end = max(match.subject_start, match.subject_end)
            for i in range(start, end):
                cov[i] = 1

print ''


coved = 0
total = 0
for name in covs:
    coved += sum(covs[name])
    total += len(covs[name])
    f = sum(covs[name]) / float(len(covs[name]))
    #print name, sum(covs[name]), len(covs[name]), f

print coved, total, coved / float(total), sys.argv[1], sys.argv[2], MIN_QUERY_LEN
