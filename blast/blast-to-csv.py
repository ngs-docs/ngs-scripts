import sys
import csv
import blastparser
 
# open the output file for reading
fp = open(sys.argv[1])
 
# send output as comma-separated values to stdout
output = csv.writer(sys.stdout)
 
# parse BLAST records
for record in blastparser.parse_fp(fp):
    for hit in record:
        for match in hit.matches:
            # output each match as a separate row
            row = [record.query_name, hit.subject_name, match.score,
                   match.expect]
            output.writerow(row)
