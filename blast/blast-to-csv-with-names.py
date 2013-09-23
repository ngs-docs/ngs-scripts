import sys
import csv
import blastparser
import screed

def load_names(filename):
    d = {}
    for record in screed.open(filename):
        if record.name.startswith('gi|'):
           ident = record.name.split('|', 2)[2]
        else:
           ident = record.name
        d[ident] = record.description
    return d
 
# open the output file for reading
query_seqs = sys.argv[1]
against_seqs = sys.argv[2]

fp = open(sys.argv[3])

print >>sys.stderr, "reading query seq names from", query_seqs
query_db = load_names(query_seqs)
print >>sys.stderr, "reading against seq names from", against_seqs
against_db = load_names(against_seqs)
 
# send output as comma-separated values to stdout
output = csv.writer(sys.stdout)
 
# parse BLAST records
print >>sys.stderr, 'parsing BLAST output'
for record in blastparser.parse_fp(fp):
    query_name = record.query_name
    for hit in record:
        for match in hit.matches:
            if query_name.startswith('gi|'):
               query_name= query_name.split('|', 2)[2]
#               print >>sys.stderr, '**', query_descr
            query_descr = query_db.get(query_name, "")
            against_descr = against_db.get(hit.subject_name, "")
            
            # output each match as a separate row
            row = [query_name, query_descr,
                   hit.subject_name, against_descr, match.score,
                   match.expect]
            output.writerow(row)
