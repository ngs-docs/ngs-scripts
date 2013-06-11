#! /usr/bin/env python
"""
Yet Another BLAST parser for NCBI BLAST output.

Goals:

 - nice introspection
 - nice Pythonic accessibility
 - maintainability in the face of changing NCBI BLAST output formats

Sample usage: ::

   for record in parse_file('blast_output.txt'):
      print '-', record.query_name, record.database.name
      for hit in record.hits:
         print '--', hit.subject_name, hit.subject_length
         print '  ', hit.total_score, hit.total_expect
         for submatch in hit:
            print submatch.expect, submatch.bits
            
            print submatch.query_sequence
            print submatch.alignment
            print submatch.subject_sequence

Here, 'submatch' is a BlastObjectSubmatch; 'hit' is a BlastSubjectHits;
'record' is a BlastQuery; and 'record.database' is a BlastDatabase.  See
the docstrings below for attributes available on these objects.

Author: C. Titus Brown <titus@caltech.edu>
"""

__version__ = 0.2

__all__ = ['BlastParser', 'parse_fp', 'parse_file', 'parse_string',
           'open_shelf']
__docformat__ = 'restructuredtext'

import math
from cStringIO import StringIO
import parse_blast

###

class BlastSubjectSubmatch(object):
    """
    BlastSubjectSubmatch.
    
    A specific submatch (score/alignment) of a query sequence to a
    subject sequence.

    Attributes:

     - expect
     - frame1
     - frame2
     - score
     - query_start
     - query_end
     - subject_start
     - subject_end
     - query_sequence
     - subject_sequence

    Usage: ::

        print submatch_obj.expect
        
    (etc.)
     
    """
#    __slots__ = ['expect', 'frame1', 'frame2',
#                 'query_start', 'query_end', 'query_sequence',
#                 'subject_start', 'subject_end', 'subject_sequence', 'identity']
    
    def __init__(self, expect, frame1, frame2,
                 q_start, q_end, q_seq, s_start, s_end, s_seq, identity, score):
        self.expect = math.pow(10, -expect)
        self.frame1 = frame1
        self.frame2 = frame2
        self.query_start = q_start
        self.query_end = q_end
        self.query_sequence = q_seq

        self.subject_start = s_start
        self.subject_end = s_end
        self.subject_sequence = s_seq
        self.score = score

    def __repr__(self):
        return "<BlastSubjectSubmatch(expect=%g, query %d-%d, subject %d-%d))>"\
               % (self.expect, self.query_start, self.query_end,
                self.subject_start, self.subject_end)

class BlastSubjectHits(object):
    """
    BlastSubjectHits.

    A list of all of the matches between a query sequence and a subject
    sequence.

    Attributes:
     * subject_name -- name of subject sequence.
     * matches -- list of BlastSubjectSubmatch objects.

    Usage: ::

        print hits_object.subject_name
        for match in hits_object:
           print match
    """
#    __slots__ = ['subject_name', 'matches' ]
    def __init__(self, subject_name, matches):
        self.subject_name = str(subject_name)
        self.matches = matches

    def __getitem__(self, i):
        return self.matches[i]

    def __len__(self):
        return len(self.matches)

    def __repr__(self):
        seqname = build_short_sequence_name(self.subject_name)
        return "<BlastSubjectHits(%s, %d matches)>" % (seqname, len(self))

class BlastQuery(object):
    """
    A BLAST query (single sequence against database) containing all results.
    
    Attributes:

      * query_name -- name of query sequence (following 'Query=').
      * hits -- a list of BlastSubjectHits, containing each match + alignment.
      
    Usage: ::

        print query_object.query_name
        for hits_object in query_object:
           print hits_object.subject_name
    """
#    __slots__ = ['query_name', 'hits' ]
    def __init__(self, query_name, hits):
        self.query_name = query_name
        self.hits = list(hits)

    def __repr__(self):
        query_short = build_short_sequence_name(self.query_name)
        return "<BlastQuery(%s (%d hits))>" % (query_short, len(self.hits))

    def __len__(self):
        return len(self.hits)

    def __getitem__(self, i):
        return self.hits[i]

class _BlastShelf(object):
    def __init__(self, filename, mode='r'):
        from shelve import BsdDbShelf
        from bsddb import btopen

        _db = btopen(filename, 'r')
        self.db = BsdDbShelf(_db)

    def __iter__(self):
        db = self.db
        last_k, _ = db.last()
        k, v = db.first()
        while k != last_k:
            yield k, v
            k, v = db.next()
        yield k, v

def open_shelf(filename, mode='r'):
    from shelve import BsdDbShelf
    from bsddb import btopen

    return _BlastShelf(filename, mode)

def parse_file(filename):
    """
    Parse records from a given file; 'filename' is the path to the file.
    """
    b = BlastParser()
    for record in b.parse_file(filename):
        yield record

def parse_fp(fp, **kw):
    """
    Parse records out of the given file handle.
    """
    b = BlastParser()
    
    for record in b.parse_fp(fp, **kw):
        yield record

def parse_string(s):
    """
    Parse records out of a string buffer.
    """
    fp = StringIO(s)
    b = BlastParser()

    for record in b.parse_fp(fp):
        yield record

class _PygrBlastHitParser(parse_blast.BlastHitParser):
    def generate_intervals(self):
        yield self.query_id, self.subject_id, \
              BlastSubjectSubmatch(self.e_value,
                                   None,
                                   None,
                                   self.query_start,
                                   self.query_end,
                                   self.query_seq,
                                   self.subject_start,
                                   self.subject_end,
                                   self.subject_seq,
                                   self.identity_percent,
                                   self.blast_score)
    
class BlastParser(object):
    """
    BlastParser objects coordinate the use of pyparsing parsers to
    parse complete BLAST records.

    Attributes:

      * blast_record -- an individual BLAST record; returns BlastQuery object.
      * blast_output -- list of records; returns list of BlastQuery objects.

    Methods:

      * reset() -- clear the blast parser of persistent information.
      * parse_string(s)
      * parse_file(filename)
      * parse_fp(fp)
    """
    def __init__(self):
        self.p = _PygrBlastHitParser()

    def parse_file(self, filename):
        fp = open(filename)
        for record in self.parse_fp(fp):
            yield record

    def parse_fp(self, fp):

        subjects = []
        matches = []

        cur_query = None
        cur_subject = None
        
        for query_id, subject_id, submatch in self.p.parse_file(fp):
            if cur_subject != subject_id or cur_query != query_id:
                if matches:
                    assert cur_subject
                    subject_hits = BlastSubjectHits(cur_subject, matches)
                    subjects.append(subject_hits)
                    matches = []

                cur_subject = subject_id
                
            if cur_query != query_id:
                if cur_query:
                    assert subjects, cur_query
                    yield BlastQuery(cur_query, subjects)
                    subjects = []

                cur_query = query_id

            matches.append(submatch)

        if matches:
            subjects.append(BlastSubjectHits(cur_subject, matches))
            
        if subjects:
            yield BlastQuery(cur_query, subjects)


def build_short_sequence_name(name, max_len=20):
    if len(name) < max_len:
        return name

    name_l = name.split()
    if len(name_l) > 1:
        return build_short_sequence_name(name_l[0], max_len)

    name = name_l[0]
    if len(name) > max_len:
        name = name[:max_len-3] + '...'
    return name

#####

if __name__ == '__main__':
    import sys
    from shelve import BsdDbShelf
    from bsddb import btopen
    from optparse import OptionParser

    ### read command line parameters

    parser = OptionParser()
    parser.add_option('-z', '--zlib-compressed', action="store_true",
                      dest="zlib_compressed",
                      help="read gzipped BLAST output file")

    parser.add_option('-n', '--ignore-empty-hits', action="store_true",
                      dest="ignore_empty_hits",
                      help="ignore BLAST hits with no results")

    (options, args) = parser.parse_args()

    (blast_file, output_file) = args

    ### open blast file, open/create database r/w

    if options.zlib_compressed:
        import gzip
        blast_fp = gzip.open(blast_file)
    else:
        blast_fp = open(blast_file)

    _db = btopen(output_file, 'c')
    db = BsdDbShelf(_db)

    ### go!

    for n, record in enumerate(parse_fp(blast_fp,
                                        ignore_no_hits=options.ignore_empty_hits)):
        if n % 100 == 0:
            print '...', n

        if options.ignore_empty_hits and not record:
            continue

        name = record.query_name
        db[name] = record

    print 'read %d records total' % (n + 1,)
