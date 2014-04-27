#! /usr/bin/python3

import argparse
import csv
import random
import re
import sys
import math
from collections import defaultdict

# ================
# Argument Parsing
# ================

def parse():
    parser = argparse.ArgumentParser(prog='smof')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

    subparsers = parser.add_subparsers()

    # FSTAT
    fstat_parser = subparsers.add_parser(
        'fstat',
        help="Provides total counts for file sequence characters")
    fstat_parser.set_defaults(func=fstat)

    # UNMASK
    unmask_parser = subparsers.add_parser(
        'unmask',
        help="Converts all letters to uppercase")
    unmask_parser.add_argument(
        '-x', '--to-x',
        help="Convert lower case letters to X",
        action='store_true',
        default=False)
    unmask_parser.set_defaults(func=unmask)

    # HSTAT
    hstat_parser = subparsers.add_parser(
        'hstat',
        help="Extract info from headers (also seq length)")
    hstat_parser.add_argument(
        'fields',
        help="Header fields to write to csv",
        nargs='+')
    hstat_parser.add_argument(
        '--length',
        help="Report length of each sequence",
        action='store_true',
        default=False)
    hstat_parser.set_defaults(func=hstat)

    # IDSEARCH
    idsearch_parser = subparsers.add_parser(
        'idsearch',
        help='Find sequences by field/value pair')
    idsearch_parser.add_argument(
        'field',
        help="Header field (e.g. 'gi' or 'locus')")
    idsearch_parser.add_argument(
        'value',
        help="Header field value (e.g. '15237703' or 'AT5G64430')")
    idsearch_parser.set_defaults(func=idsearch)

    # PRETTYPRINT
    prettyprint_parser = subparsers.add_parser(
        'prettyprint',
        help='Prints fasta file in neat columns')
    prettyprint_parser.add_argument(
        'cwidth',
        help='Output column width',
        type=int,
        default=60)
    prettyprint_parser.set_defaults(func=prettyprint)

    # QSTAT
    qstat_parser = subparsers.add_parser(
        'qstat',
        help="Gathers statistics on each sequence")
    qstat_parser.add_argument(
        '-f', '--fields',
        help="Header fields by which each sequence is identified in output csv file",
        nargs='+',
        default=[])
    qstat_parser.add_argument(
        '-m', '--countmasked',
        help="Count the number of masked (lowcase) characters",
        action='store_true',
        default=False)
    qstat_parser.add_argument(
        '-i', '--ignorecase',
        help="Ignores case when counting characters",
        action='store_true',
        default=False)
    qstat_parser.add_argument(
        '-s', '--start-offset',
        help='Number of letters to ignore at beginning (default=0)',
        type=int,
        default=0
    )
    qstat_parser.add_argument(
        '-e', '--end-offset',
        help='Number of letters to ignore at end (default=0)',
        type=int,
        default=0
    )
    qstat_parser.set_defaults(func=qstat)

    # RSEARCH
    rsearch_parser = subparsers.add_parser(
        'rsearch',
        help="Mass extract of sequences from regex header search")
    rsearch_parser.add_argument(
        '-p', '--pattern',
        help="Perl regex pattern containing one set of parentheses")
    rsearch_parser.add_argument(
        '-v', '--values',
        help="The acceptable values for the parenthesized group",
        nargs='+')
    rsearch_parser.add_argument(
        '-f', '--valuefile',
        help="Read values from file")
    rsearch_parser.set_defaults(func=rsearch)

    # SAMPLE
    sample_parser = subparsers.add_parser(
        'sample',
        help="Randomly select entries from fasta file")
    sample_parser.add_argument(
        'n',
        help="Sample size",
        type=int,
        default=1)
    sample_parser.set_defaults(func=sample)

    # SORT
    sort_parser = subparsers.add_parser(
        'sort',
        help="Sort sequences by given fields")
    sort_parser.add_argument(
        'fields',
        help="Header fields by which to sort sequences",
        nargs='+')
    sort_parser.set_defaults(func=sort)

    # SEARCH
    search_parser = subparsers.add_parser(
        'search',
        help='Search for pattern')
    search_parser.add_argument(
        'pattern',
        help='Perl regular expression search pattern')
    search_parser.add_argument(
        '-i', '--invert',
        help="Drop all not matching the pattern",
        action='store_false',
        default=True
    )
    search_parser.add_argument(
        '-s', '--seq',
        help='Search for pattern in the sequence',
        action='store_true',
        default=False
    )
    search_parser.set_defaults(func=search)

    # SPLIT
    split_parser = subparsers.add_parser(
        'split',
        help='Split a multifasta file into k smaller filers'
    )
    split_parser.add_argument(
        '-n', '--nfiles',
        help='Number of output files'
    )
    split_parser.add_argument(
        '-p', '--prefix',
        help='Prefix for output files',
        default='xxx'
    )
    split_parser.set_defaults(func=split)

    # SUBSEQ
    subseq_parser = subparsers.add_parser(
        'subseq',
        help="Extract subsequence from each entry")
    subseq_parser.add_argument(
        'bounds',
        help="from and to values (indexed from 1)",
        nargs=2,
        type=int)
    subseq_parser.add_argument(
        '-r', '--revcomp',
        help='Take the reverse complement if bounds[0] > bounds[1]',
        action='store_true',
        default=False)
    subseq_parser.set_defaults(func=subseq)

    # FSUBSEQ
    fsubseq_parser = subparsers.add_parser(
        'fsubseq',
        help="Mass extraction of subsequences from first fasta entry")
    fsubseq_parser.add_argument(
        'file',
        help="File containing bounds for subsequence extraction")
    fsubseq_parser.add_argument(
        '-r', '--revcomp',
        help='Take the reverse complement if bounds[0] > bounds[1]',
        action='store_true',
        default=False)
    fsubseq_parser.add_argument(
        '-p', '--pattern-index',
        help='Index of regex pattern in each row',
        type=int,
        default=-1)
    fsubseq_parser.set_defaults(func=fsubseq)

    # FASTA2CSV
    fasta2csv_parser = subparsers.add_parser(
        'fasta2csv',
        help="Converts a fasta file to a csv, two columns (header|sequence)")
    fasta2csv_parser.add_argument(
        '-d', '--delimiter',
        help="Set delimiter (',' by default)",
        default=',')
    fasta2csv_parser.add_argument(
        '-r', '--header',
        help='Write header (default=False)',
        action='store_true',
        default=False)
    fasta2csv_parser.add_argument(
        '-f', '--fields',
        help='Extract given fields from the header',
        nargs='+')
    fasta2csv_parser.set_defaults(func=fasta2csv)

    # PERM
    perm_parser = subparsers.add_parser(
        'perm',
        help="Randomly order sequence by words of length w"
    )
    perm_parser.add_argument(
        '-w', '--word-size',
        help='Size of each word (default=1)',
        type=int,
        default=1
    )
    perm_parser.add_argument(
        '-s', '--start-offset',
        help='Number of letters to ignore at beginning (default=0)',
        type=int,
        default=0
    )
    perm_parser.add_argument(
        '-e', '--end-offset',
        help='Number of letters to ignore at end (default=0)',
        type=int,
        default=0
    )
    perm_parser.add_argument(
        '-f', '--field',
        help="Header field (e.g. 'gi' or 'locus')"
    )
    perm_parser.set_defaults(func=perm)

    # Simplifyheader
    simplifyheader_parser = subparsers.add_parser(
        'simplifyheader',
        help="Reduce header to given fields"
    )
    simplifyheader_parser.add_argument(
        'fields',
        help="Fields to retain",
        nargs='+'
    )
    simplifyheader_parser.set_defaults(func=simplifyheader)

    # REVERSE
    reverse_parser = subparsers.add_parser(
        'reverse',
        help="Reverse each sequence")
    reverse_parser.set_defaults(func=reverse)

    # TRANSLATE
    translate_parser = subparsers.add_parser(
        'translate',
        help="Translates DNA (wrapper for EMBOSS transeq)"
    )
    translate_parser.add_argument(
        '-f', '--frame',
        help="See EMBOSS transeq help",
        default=1
    )
    translate_parser.add_argument(
        '-t', '--table',
        help="See EMBOSS transeq help",
        default=0
    )
    translate_parser.add_argument(
        '-r', '--regions',
        help="See EMBOSS transeq help"
    )
    translate_parser.add_argument(
        '-m', '--trim',
        help="See EMBOSS transeq help",
        action='store_true',
        default=False
    )
    translate_parser.add_argument(
        '-c', '--clean',
        help="See EMBOSS transeq help",
        action='store_true',
        default=False
    )

    if(len(sys.argv) == 1):
        parser.print_help()
        raise SystemExit

    args = parser.parse_args()

    return(args)


# =================
# CLASS DEFINITIONS
# =================

class FSeqGenerator:
    def __init__(self, fh=sys.stdin):
        self.fh = fh

    def next(self):
        seq_list = []
        header = ''
        for line in self.fh:
            line = line.strip()
            if ">" in line:
                if(seq_list):
                    yield FSeq(header, ''.join(seq_list))
                seq_list = []
                header = line.split('>')[1]
            else:
                seq_list.append(line)
        yield FSeq(header, ''.join(seq_list))

class FSeq:
    # The translator for taking reverse complements
    revcomp_translator = str.maketrans('acgtACGT', 'tgcaTGCA')

    def __init__(self, header, seq):
        self.seq = seq
        self.header = header
        self.headerfields = {}

    def parse_header(self):
        '''
        Parses headers of the format:
           '>field_1|value_1|field_2|value_2|...|field_k|value_k| description (optional)
        For example:
            >locus|AT3G01015|taxon|3702|gi|186509637|gb|NP_186749.2| TPX2 targeting protein
        '''
        fields = self.header.split('|')
        # If there is no delimitation with '|', set whatever there is to field
        # 'header'
        if(len(fields) == 1):
            self.headerfields['header'] = fields[0]
        else:
            if(len(fields) % 2  == 1 and fields[-1]):
                self.headerfields['desc'] = fields.pop().strip()
            for i in range(0, len(fields), 2):
                self.headerfields[fields[i].strip()] = fields[i+1].strip()

    def has_field(self, field):
        if(not self.headerfields):
            self.parse_header()
        if field in self.headerfields:
            return True
        else:
            return False

    def field_is(self, field, value):
        value = str(value)
        if(not self.headerfields):
            self.parse_header()
        try:
            if(self.headerfields[field] == value):
                return True
            else:
                return False
        except:
            print("Header lacks field {}".format(field), file=sys.stderr)
            print(self.header, file=sys.stderr)

    def getvalue(self, field, quiet=False):
        if(not self.headerfields):
            self.parse_header()
        try:
            return(self.headerfields[field])
        except:
            if(not quiet):
                print("Header lacks field '{}'".format(field), file=sys.stderr)
                print(self.header, file=sys.stderr)
            raise SystemExit

    def countmasked(self):
        '''
        Masked sequences are lowercase, unmasked are uppercase. This function
        simply counts the number of lowercase characters.
        '''
        # ASCII values: (97-122 for a-z), (65-90 for A-Z)
        nmasked = sum(c > 90 for c in self.seq.encode('ascii'))
        return(nmasked)

    def subseq(self, a, b):
        try:
            return(self.seq[a:b])
        except IndexError as e:
            print("Cannot extract subsequence ({}, {}) from".format(a,b))
            print('>' + self.header)
        except:
            print("Unknown error in extraction of subsequence ({}, {}) from".format(a,b))
            print('>' + self.header)

    def print(self, column_width=80):
        print('>' + self.header)
        for i in range(0, len(self.seq), column_width):
            print(self.seq[i:i + column_width])
    def get_pretty_string(self, column_width=80):
        out = ['>' + self.header]
        for i in range(0, len(self.seq), column_width):
            out.append(self.seq[i:i + column_width])
        outstr = '\n'.join(out)
        return(outstr)

    @classmethod
    def getrevcomp(cls, seq):
        return(seq[::-1].translate(FSeq.revcomp_translator))

class Stat:
    def __init__(self, fields=[]):
        self.counts = defaultdict(int)
        self.fields = {}
        self.length = 0

    def update_counts(self, seq, ignoreCase=False):
        s = seq.seq
        if(ignoreCase):
            s = s.upper()
        for c in s:
            self.counts[c] += 1
        self.length += len(s)

    def update_fields(self, seq, fields):
        for field in fields:
            self.fields[field] = seq.getvalue(field)

    def getdict(self, charset=None, getlength=False):
        if(not charset):
            charset = self.counts.keys()
        out = {}
        for k in charset:
            out[k] = self.counts[k]
        for k, v in self.fields.items():
            out[k] = v
        if(getlength):
            out['length'] = self.length
        return(out)

    @classmethod
    def getcharset(cls, statlist):
        charset = set()
        for s in statlist:
            charset.update(s.counts.keys())
        return(charset)


# =================
# UTILITY FUNCTIONS
# =================

def csvrowGenerator(filename):
    try:
        f = open(filename, 'r')
        dialect = csv.Sniffer().sniff(f.read(1024))
        f.seek(0)
        reader = csv.reader(f, dialect)
    except IOError as e:
        print(e)
        sys.exit(0)
    for row in reader:
        yield row


# ====================
# ONE-BY-ONE FUNCTIONS
# ====================

def fasta2csv(args):
    w = csv.writer(sys.stdout, delimiter=args.delimiter)
    out = []
    if(args.header):
        if(args.fields):
            w.writerow(args.fields + ['seq'])
        else:
            w.writerow(['header', 'seq'])
    for seq in FSeqGenerator().next():
        if(args.fields):
            elements = [seq.getvalue(field) for field in args.fields]
        else:
            elements = [seq.header]
        w.writerow(tuple(elements + [seq.seq]))

def fstat(args):
    stats = Stat()
    nseqs = 0
    nchars = 0
    for seq in FSeqGenerator().next():
        nseqs += 1
        nchars += len(seq.seq)
        stats.update_counts(seq)
    for k, v in sorted(stats.getdict().items(), key=lambda x: x[1], reverse=True):
        print("{}: {} {}".format(k, v, round(v/nchars, 4)))
    print("nseqs: {}".format(nseqs))
    print("nchars: {}".format(nchars))
    print("mean length: {}".format(round(nchars/nseqs, 4)))

def hstat(args):
    ''' Writes chosen header and seq length data to csv '''
    fieldnames = list(args.fields)
    if(args.length):
        fieldnames.append('length')
    w = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    w.writeheader()
    for seq in FSeqGenerator().next():
        row = {}
        for field in args.fields:
            row[field] = seq.getvalue(field)
        if(args.length):
            row['length'] = len(seq.seq)
        w.writerow(row)

def idsearch(args):
    ''' Print entries whose headers contain a field with a given value '''
    # TODO make mass search
    for seq in FSeqGenerator().next():
        if(seq.field_is(args.field, args.value)):
            seq.print()

def perm(args):
    w = args.word_size
    start = args.start_offset
    end = args.end_offset
    for seq in FSeqGenerator().next():
        s = seq.seq
        L = len(s)
        prefix = s[0:start]
        suffix = s[L-end:L]
        rseq = s[start:L-end]
        M = len(rseq)
        words = list(rseq[i:i+w] for i in range(0, M - w + 1, w))
        words.append(rseq[(M - M % w):M])
        random.shuffle(words)
        out = ''.join(prefix + ''.join(words) + suffix)
        if(args.field):
            header='|'.join((args.field, seq.getvalue(args.field), 'start', str(start), 'end', str(end), 'word_size', str(w)))
        else:
            header=seq.header
        FSeq(header, out).print()

def prettyprint(args):
    ''' Print each sequence with even columns of length cwidth '''
    for seq in FSeqGenerator().next():
        seq.print(args.cwidth)

def qstat(args):
    class Result:
        def __init__(self, seq, args):
            self.stats = self._getstats(seq, args)
            self.other = self._getother(seq, args)
            try:
                self.fields = {x:seq.getvalue(x) for x in args.fields}
            except:
                self.fields = {}

        def _getstats(self, seq, args):
            seqstats = Stat()
            seqstats.update_counts(seq, ignoreCase=args.ignorecase)
            return(seqstats)

        def _getother(self, seq, args):
            others = {}
            if(args.countmasked):
                others['masked'] = seq.countmasked()
            others['length'] = len(seq.seq)
            return(others)

        def getdict(self, charset):
            d = {k:v for k,v in self.fields.items()}
            for k,v in self.other.items():
                d[k] = v
            for k, v in self.stats.getdict(charset=charset, getlength=False).items():
                d[k] = v
            return(d)

        def getfields(self, charset):
            names = sorted(self.fields.keys()) + \
                    sorted(self.other.keys()) + \
                    sorted(charset)
            return(names)

    results = []
    # Fill a list with Stat objects
    for seq in FSeqGenerator().next():
        # Trim the beginning and end of the sequence, as specified
        trunc_seq = seq.seq[args.start_offset:len(seq.seq) - args.end_offset]
        seq = FSeq(seq.header, trunc_seq)
        results.append(Result(seq, args))
    # Set of all unique characters seen in the fasta file
    charset = Stat.getcharset([result.stats for result in results])
    # Rownames for the csv file
    fieldnames = results[0].getfields(charset)
    w = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    w.writeheader()
    [ w.writerow(r.getdict(charset)) for r in results ]

def reverse(args):
    ''' Reverse each sequence '''
    for seq in FSeqGenerator().next():
        rseq = FSeq(seq.header, seq.seq[::-1])
        rseq.print()

def rsearch(args):
    '''
    Extracts sequences matching a certain pattern
    '''
    values = set()
    if(args.valuefile):
        with open(args.valuefile) as f:
            filevalues = [line.rstrip('\n') for line in f.readlines()]
            values.update(filevalues)
    if(args.values):
        values.update(args.values)
    pat = re.compile(args.pattern)
    for seq in FSeqGenerator().next():
        m = re.search(pat, seq.header)
        # If no match
        if(not m):
            continue
        try:
            # If at least one group defined
            match = m.group(1)
        except:
            # If no groups defined
            match = m.group(0)
        if(match in values):
            seq.print()

def search(args):
    ''' Print entries whose headers contain a given pattern '''
    # TODO add highlighting
    prog = re.compile(args.pattern)
    for seq in FSeqGenerator().next():
        if(args.seq):
            searchseq = seq.seq
        else:
            searchseq = seq.header
        hasmatch = prog.search(searchseq) is None
        if(hasmatch and not args.invert):
            seq.print()
        elif(not hasmatch and args.invert):
            seq.print()

def simplifyheader(args):
    if(hasattr(args.fields, '__upper__')):
        args.fields = (args.fields, )
    for seq in FSeqGenerator().next():
        values = [seq.getvalue(field) for field in args.fields]
        pairs = ['|'.join((args.fields[i], values[i])) for i in range(len(values))]
        header = '|'.join(pairs)
        FSeq(header, seq.seq).print()

def split(args):
    k = int(args.nfiles)
    p = args.prefix
    seqs = []
    for seq in FSeqGenerator().next():
        seqs.append(seq)
    for i in range(0, k):
        begin = i * (len(seqs) // k + 1)
        end = min(len(seqs), (i+1) * (len(seqs) // k + 1))
        outfile = p + str(i) + '.fasta'
        with open(outfile, 'w+') as fo:
            for seq in (seqs[x] for x in range(begin, end)):
                fo.write(seq.get_pretty_string() + '\n')

def subseq(args):
    ''' Index starting from 1 (not 0) '''
    for seq in FSeqGenerator().next():
        a = args.bounds[0]
        b = args.bounds[1]
        # If a > b, take reverse complement if that option is enabled,
        # if not die
        if(a > b):
            if(args.revcomp):
                newseq = FSeq.getrevcomp(seq.seq[b-1:a])
            else:
                print('Lower bound cannot be greater than upper bound ' + \
                      '(do you want reverse complement? See options)', file=sys.stderr)
                sys.exit()
        # If b >= a, this is a normal forward sequence
        else:
            newseq = seq.seq[a-1:b]
        FSeq(seq.header, newseq).print()

def fsubseq(args):
    '''
    Extracts many subsequences. Positional argument 'file'
    '''
    bounds = defaultdict(list)
    for row in csvrowGenerator(args.file):
        if(args.pattern_index >= 0):
            try:
                pat = row[args.pattern_index]
                del row[args.pattern_index]
            except IndexError as e:
                print(e)
                print("Given pattern index doesn't exist, fuck you!")
                sys.exit()
        else:
            pat = '.*'
        try:
            a,b = [int(x) for x in row]
        except TypeError as e:
            print(e)
            print("Bounds must be pair of integers")
            sys.exit()
        bounds[pat].append((a,b))

    for seq in FSeqGenerator().next():
        for pat in bounds.keys():
            if(not re.search(pat, seq.header)):
                continue
            for a,b in bounds[pat]:
                # If a > b, take reverse complement if that option is enabled,
                # if not die
                if(a > b):
                    if(args.revcomp):
                        newseq = FSeq.getrevcomp(seq.seq[(b-1):a])
                    else:
                        print('Lower bound cannot be greater than upper bound ' + \
                            '(do you want reverse complement? See options)', file=sys.stderr)
                        sys.exit()
                # If b >= a, this is a normal forward sequence
                else:
                    newseq = seq.seq[(a-1):b]
                FSeq(seq.header, newseq).print()

def unmask(args):
    ''' Converts to upcase '''
    for seq in FSeqGenerator().next():
        if(args.to_x):
            unmasked_seq = FSeq(seq.header, re.sub('[a-z]', 'X', seq.seq))
        else:
            unmasked_seq = FSeq(seq.header, seq.seq.upper())
        unmasked_seq.print()


# ==============
# FULL FUNCTIONS
# ==============

def sample(args):
    ''' Randomly sample n entries from input file '''
    seqs = [s for s in FSeqGenerator().next()]
    sample_indices = random.sample(range(len(seqs)), min(len(seqs), args.n))
    [seqs[i].print() for i in sample_indices]
def sort(args):
    seqs = [s for s in FSeqGenerator().next()]
    seqs.sort(key=lambda x: list(x.getvalue(y) for y in args.fields))
    [s.print() for s in seqs]


# =======
# EXECUTE
# =======

if __name__ == '__main__':
    args = parse()
    args.func(args)
