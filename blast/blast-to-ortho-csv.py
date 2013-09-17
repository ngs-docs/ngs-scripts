import sys
import csv
import blastparser
import screed

def collect_best_hits(filename):
    d = {}
    for n, record in enumerate(blastparser.parse_fp(open(filename))):
        if n % 25000 == 0:
            print '...', filename, n
        best_score = None
        for hit in record.hits:
            for match in hit.matches:
                query = record.query_name
                if query.startswith('gi'):
                    query = query.split('|', 2)[2]
                subject = hit.subject_name

                score = match.score

                # only keep the best set of scores for any query
                if best_score and best_score > score:
                    continue
                best_score = score

                x = d.get(query, [])
                x.append((subject, score))
                d[query] = x

            if best_score and best_score != score:
                break
    return d

def parse_ncbi_query(name):
    name = name.split('|')[2:]
    name = '|'.join(name)
    return name

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

ab = sys.argv[3]
ba = sys.argv[4]

print >>sys.stderr, "reading query seq names from", query_seqs
query_db = load_names(query_seqs)
print >>sys.stderr, "reading against seq names from", against_seqs
against_db = load_names(against_seqs)
 
# send output as comma-separated values to stdout
output = csv.writer(sys.stdout)
 
# parse BLAST records
print >>sys.stderr, 'parsing BLAST output', ab
ab_dict = collect_best_hits(ab)
print >>sys.stderr, 'parsing BLAST output', ba
ba_dict = collect_best_hits(ba)

print ab_dict.items()[:5]
print ba_dict.items()[:5]

print 'calculating reciprocal best hits'
dd = {}
ee = {}
for k in ab_dict:
    v = map(lambda x: x[0], ab_dict[k])

    for k2 in v:
        v2 = map(lambda x: x[0], ba_dict.get(k2, []))

        if k in v2:
            dd[k] = k2
            ee[k2] = k

for k in dd:
    v = dd[k]

    query_descr = query_db.get(k, "")
    against_descr = against_db.get(v, "")

    # output each match as a separate row
    row = [k, query_descr, v, against_descr]
    output.writerow(row)
