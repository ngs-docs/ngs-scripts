#! /usr/bin/env python
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('memefile')

    args = parser.parse_args()
    
    data = open(args.memefile).read()
    
    n = 0
    loc = data.find('\nMOTIF ')
    while loc != -1:
        n += 1
        table_loc = data.find('\nSequence name', loc)
        assert table_loc != -1
        table_loc = data.find('\nig:', table_loc + 1)
        assert table_loc != -1

        end_table = data.find('\n----', table_loc)
        assert end_table != -1

        table = data[table_loc + 1:end_table]
        table = table.splitlines()

        sites = []
        for line in table:
            line = line.split()
            sites.append(line[4])

        filename = 'motif%d.sites' % n
        print 'writing', filename
        fp = open(filename, 'w')
        fp.write("\n".join(sites))
        fp.close()

        loc = data.find('\nMOTIF ', loc + 1)

    print n

if __name__ == '__main__':
    main()
