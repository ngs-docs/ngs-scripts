#! /usr/bin/env python
import argparse
import load_gff
import screed

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('genome')
    parser.add_argument('gff_files', nargs='+')

    args = parser.parse_args()
    
    genomefile = args.genome
    gff_d = load_gff.load(args.gff_files)

    chr_operons = load_gff.make_operons(gff_d)
    chr_intergenic = load_gff.make_intergenic(genomefile, chr_operons)

    seqs = {}
    for record in screed.open(genomefile):
        name = record.name
        if name.startswith('gi|'):
            name = name.split('|')[3]
        seqs[name] = record.sequence

    for chr in chr_intergenic:
        seq = seqs[chr]
        intergenic = chr_intergenic[chr]
        for ig in intergenic:
            (start, name1, stop, name2) = ig
            if stop - start < 20:
                continue
            print '>ig:%s:%s\n%s' % (name1, name2, seq[start:stop])

if __name__ == '__main__':
    main()
