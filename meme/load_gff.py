import sys
import screed

def load(filenames):
    d = {}
    for filename in filenames:
        for line in open(filename):
            if line.startswith('#'):
                continue

            chr, origin, typ, start, end, _, orient, _, info = line.split('\t')
            if typ != 'CDS':
                continue

            info = info.strip().split(';')
            info = [ x.split('=') for x in info ]
            info = dict(info)

            start, end = int(start), int(end)

            ident = info['Name']
            assert ident not in d, ident
            d[ident] = (chr, start, end, orient, info)

    return d

def make_operons(gff_d, overlap=20):
    chrs = {}
    for ident in gff_d:
        (chr, start, end, orient, info) = gff_d[ident]
        x = chrs.get(chr, [])
        x.append((start, end, orient, info))
        chrs[chr] = x

    chr_operons = {}

    # turn genes into operons, by chr
    for chr in chrs:
        x = chrs[chr]
        x.sort()

        idx = 0
        
        this_op = []
        op_start, op_end, op_orient, op_info = x[idx]
        this_op.append((op_start, op_end, op_orient, op_info))

        if len(x) > 1:

            idx += 1
            next_op_start, next_op_end, next_op_orient, next_op_info = \
                x[idx]

        operons = {}
        while idx < len(x) - 1:
            while next_op_start < op_end + overlap and \
                    op_orient == next_op_orient:
                this_op.append((next_op_start, next_op_end, next_op_orient,
                                next_op_info))
                idx += 1
                if idx == len(x) - 1: break

                op_start, op_end, op_orient, op_info = \
                    next_op_start, next_op_end, next_op_orient, next_op_info
                next_op_start, next_op_end, next_op_orient, next_op_info = \
                    x[idx]

            operon_name = ";".join([ op[3]['Name'] for op in this_op ])
            operons[operon_name] = this_op

            op_start, op_end, op_orient, op_info = \
                next_op_start, next_op_end, next_op_orient, next_op_info
            this_op = [(op_start, op_end, op_orient, op_info)]
            idx += 1
            if idx >= len(x) - 1: break
            next_op_start, next_op_end, next_op_orient, next_op_info = \
                x[idx]

        operon_name = ";".join([ op[3]['Name'] for op in this_op ])
        operons[operon_name] = this_op

        chr_operons[chr] = operons

    return chr_operons

def get_operon_first_name(operon):
    return operon[0][3]['Name']

def get_operon_last_name(operon):
    return operon[-1][3]['Name']

def get_operon_start(operon):
    return operon[0][0]

def get_operon_stop(operon):
    return operon[-1][1]

def make_intergenic(dbfile, chr_operons):
    seqs = {}
    for record in screed.open(dbfile):
        name = record.name
        if name.startswith('gi|'):
            name = name.split('|')[3]
        seqs[name] = len(record.sequence)

    chr_intergenic = {}
    for chr in chr_operons:
        operons = list(sorted(chr_operons[chr].values(), key=get_operon_start))

        operon = operons[0]
        last_name0 = get_operon_first_name(operon)
        last_name1 = get_operon_last_name(operon)
        last_start = get_operon_start(operon)
        last_end = get_operon_stop(operon)

        intergenic = [(0, '-', last_start, last_name0)]

        for idx in range(1, len(operons) - 1):
            operon = operons[idx]

            name0 = get_operon_first_name(operon)
            name1 = get_operon_last_name(operon)
            start = get_operon_start(operon)
            end = get_operon_stop(operon)

            if start - last_end > 0:
                intergenic.append((last_end, last_name1, start, name0))

            last_name0 = name0
            last_name1 = name1
            last_start = start
            last_end = end

        if seqs[chr] - last_end > 0:
            intergenic.append((last_end, last_name1, seqs[chr], '-'))

        chr_intergenic[chr] = intergenic

    return chr_intergenic

if __name__ == '__main__':
    d = load(sys.argv[2:])
    print 'loaded %d CDS' % len(d)
    
    chr_operons = make_operons(d)
    for chr in chr_operons:
        print 'chr', chr, 'has', len(chr_operons[chr]), 'operons'

    chr_intergenic = make_intergenic(sys.argv[1], chr_operons)
