#! /usr/bin/python3

import csv
import random
import math
import re
import sys
import string
from collections import defaultdict

__version__ = "1.0"

# ================
# Argument Parsing
# ================

def parse(argv=None):
    import argparse
    parser = argparse.ArgumentParser(
        prog='smof',
        usage='<fastafile> | smof <subcommand> <options>',
        description='Tools for studying and manipulating fasta files')

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {}'.format(__version__))

    subparsers = parser.add_subparsers(
        metavar='[ for help on each: smof <subcommand> -h ]',
        title='subcommands')

    genusage = '<fastafile> | smof {} <options>'

    # CHSUM
    chsum_parser = subparsers.add_parser(
        'chksum',
        usage=genusage.format('chksum'),
        help="Calculate an md5 checksum for the input sequences"
    )
    chsum_parser.add_argument(
        '-i', '--ignore-case',
        help='Convert all to uppercase, before hashing',
        action='store_true',
        default=False
    )
    method = chsum_parser.add_mutually_exclusive_group(required=False)
    method.add_argument(
        '-w', '--whole-file',
        help='Calculate single md5sum for all headers and sequences (default action)',
        action='store_true',
        default=False
    )
    method.add_argument(
        '-q', '--each-sequence',
        help='Calculate md5sum for each sequence, write as TAB delimited list',
        action='store_true',
        default=False
    )
    method.add_argument(
        '-s', '--all-sequences',
        help='Calculate single md5sum for all sequences',
        action='store_true',
        default=False
    )
    method.add_argument(
        '-d', '--all-headers',
        help='Calculate single md5sum for all headers',
        action='store_true',
        default=False
    )
    chsum_parser.set_defaults(func=chsum)

    # COMPLEXITY
    complexity_parser = subparsers.add_parser(
        'complexity',
        usage=genusage.format('complexity'),
        help='Calculates the complexity of the given sequences'
    )
    complexity_parser.add_argument(
        '-k', '--alphabet-size',
        help='Number of letters in the alphabet (4 for DNA, 20 for proteins)',
        default=4
    )
    complexity_parser.add_argument(
        '-w', '--window-length',
        help='Window length (if provided, output will average of window complexities)',
        default=100
    )
    complexity_parser.add_argument(
        '-m', '--word-length',
        help='Length of each word',
        default=1
    )
    complexity_parser.add_argument(
        '-j', '--jump',
        help='Distance between adjacent windows',
        default=1
    )
    complexity_parser.add_argument(
        '-o', '--offset',
        help='Index of start point',
        default=0
    )
    complexity_parser.add_argument(
        '-d', '--drop',
        help="Drop sequence if contains this character (e.g. 'X' or 'N')",
        default=None
    )
    complexity_parser.set_defaults(func=complexity)

    # FSTAT
    fstat_parser = subparsers.add_parser(
        'fstat',
        usage=genusage.format('fstat'),
        help="Provides total counts for file sequence characters")
    fstat_parser.set_defaults(func=fstat)

    # UNMASK
    unmask_parser = subparsers.add_parser(
        'unmask',
        usage=genusage.format('unmask'),
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
        usage=genusage.format('hstat'),
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
        usage=genusage.format('idsearch'),
        help='Find sequences by field/value pair')
    idsearch_parser.add_argument(
        'field',
        help="Header field (e.g. 'gi' or 'locus')")
    idsearch_parser.add_argument(
        'value',
        help="Header field value (e.g. '15237703' or 'AT5G64430')")
    idsearch_parser.set_defaults(func=idsearch)

    # tounk
    tounk_parser = subparsers.add_parser(
        'tounk',
        usage=genusage.format('tounk'),
        help='Convert irregular characters to unknown character'
    )
    tounk_parser.add_argument(
        '-t', '--type',
        help='Sequence type [n, p]'
    )
    tounk_parser.add_argument(
        '-l', '--lc',
        help='Convert lower-case to unknown',
        action='store_true',
        default=False
    )
    tounk_parser.add_argument(
        '--nir',
        help='Nucleotide irregulars [default=(not ACGT)]',
        default=''.join(set(string.ascii_letters) - set('ACGTNacgtn'))
    )
    tounk_parser.add_argument(
        '--pir',
        help='Protein irregulars [default=BJOUZbjouz]',
        default='BJOUXZbjouxz'
    )
    tounk_parser.add_argument(
        '--nunk',
        help='Nucleotide unknown character (default=N)',
        default='N'
    )
    tounk_parser.add_argument(
        '--punk',
        help='Protein unknown character (default=X)',
        default='X'
    )
    tounk_parser.set_defaults(func=tounk)

    # PRETTYPRINT
    prettyprint_parser = subparsers.add_parser(
        'prettyprint',
        usage=genusage.format('prettyprint'),
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
        usage=genusage.format('qstat'),
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

    # RETRIEVE
    retrieve_parser = subparsers.add_parser(
        'retrieve',
        usage=genusage.format('retrieve'),
        help="Retrieve sequences with matches pattern"
    )
    retrieve_parser.add_argument(
        '-p', '--pattern',
        help="Perl regular expressions",
        nargs='+'
    )
    retrieve_parser.add_argument(
        '-P', '--patternfile',
        help="Perl regular expressions"
    )
    retrieve_parser.add_argument(
        '-g', '--groups',
        help="The acceptable values for parenthesized groups",
        nargs='*'
    )
    retrieve_parser.add_argument(
        '-G', '--groupfile',
        help="Read values from file"
    )
    retrieve_parser.add_argument(
        '-v', '--invert',
        help="Write sequences that don't match"
    )
    retrieve_parser.add_argument(
        '-q', '--match-sequence',
        help='Match sequence rather than header',
        action='store_true',
        default=False
    )
    retrieve_parser.add_argument(
        '-c', '--color',
        help='Color match',
        action='store_true',
        default=False
    )
    retrieve_parser.set_defaults(func=retrieve)

    # SAMPLE
    sample_parser = subparsers.add_parser(
        'sample',
        usage=genusage.format('sample'),
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
        usage=genusage.format('sort'),
        help="Sort sequences by given fields")
    sort_parser.add_argument(
        'fields',
        help="Header fields by which to sort sequences",
        nargs='+')
    sort_parser.set_defaults(func=sort)

    # SEARCH
    search_parser = subparsers.add_parser(
        'search',
        usage=genusage.format('search'),
        help='Search for pattern')
    search_parser.add_argument(
        'pattern',
        help='Perl regular expression search pattern')
    search_parser.add_argument(
        '-i', '--invert',
        help="Drop all not matching the pattern",
        action='store_true',
        default=False
    )
    search_parser.add_argument(
        '-q', '--seq',
        help='Search for pattern in the sequence',
        action='store_true',
        default=False
    )
    search_parser.add_argument(
        '-c', '--color',
        help='Highlight the matched sequence',
        action='store_true',
        default=False
    )
    search_parser.set_defaults(func=search)

    # SPLIT
    split_parser = subparsers.add_parser(
        'split',
        usage=genusage.format('split'),
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
        usage=genusage.format('subseq'),
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
        usage=genusage.format('fsubseq'),
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
        usage=genusage.format('fasta2csv'),
        help="Converts a fasta file to 2-column csv")
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
        usage=genusage.format('perm'),
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
        'rmfields',
        usage=genusage.format('rmfields'),
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
        usage=genusage.format('reverse'),
        help="Reverse each sequence")
    reverse_parser.set_defaults(func=reverse)

    # TRANSLATE
    translate_parser = subparsers.add_parser(
        'translate',
        usage=genusage.format('translate'),
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

    args = parser.parse_args(argv)

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
        self.colseq = ColorString()
        self.colheader = ColorString()

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
        if self.colheader.seq:
            self.colheader.print(column_width)
        else:
            print('>' + self.header)
        for i in range(0, len(self.seq), column_width):
            if self.colseq.seq:
                self.colseq.print(column_width)
                break
            else:
                print(self.seq[i:i + column_width])

    def get_pretty_string(self, column_width=80):
        out = ['>' + self.header]
        for i in range(0, len(self.seq), column_width):
            out.append(self.seq[i:i + column_width])
        outstr = '\n'.join(out)
        return(outstr)

    def seq_upper(self):
        self.seq = self.seq.upper()

    def header_upper(self):
        self.header = self.header.upper()

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

class Colors:
    HIGHLIGHT = chr(27) + '[32m'
    BACKGROUND = chr(27) + '[0m'

class ColorString:
    def __init__(self,
                 bgcolor=Colors.BACKGROUND,
                 default=Colors.HIGHLIGHT):
        self.bgcolor = bgcolor
        self.default = default
        self.seq = []

    def setseq(self, seq):
        if not self.seq:
            self.seq = [[self.bgcolor, s] for s in seq]

    def colorpos(self, pos, col=None):
        col = self.default if not col else col
        for i in pos:
            self.seq[i][0] = col

    def print(self, colwidth=None):
        lastcol = self.seq[0][0]
        for i in range(1, len(self.seq)):
            if(self.seq[i][0] == lastcol):
                self.seq[i][0] = ''
            else:
                lastcol = self.seq[i][0]
        for i in range(len(self.seq)):
            print(''.join(self.seq[i]), end='')
            if(colwidth and i % colwidth == 0 and i != 0):
                print()
        print()

    def colormatch(self, pattern, col=None):
        col = self.default if not col else col
        s = ''.join([x[1] for x in self.seq])
        for m in re.compile(pattern).finditer(s):
            a = m.start()
            b = m.start() + len(m.group())
            self.colorpos(list(range(a,b)), col)


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

def chsum(args, gen):
    '''
    Prints various hashes of the input fasta file
    '''
    import hashlib
    h = hashlib.md5()
    # Hash the sequences only (in input order)
    if(args.all_sequences):
        fun = lambda s: h.update(s.seq.encode('ascii'))
    # Hash the headers only (in input order)
    elif(args.all_headers):
        fun = lambda s: h.update(s.header.encode('ascii'))
    # Write <header>\t<sequence hash> for each sequence
    elif(args.each_sequence):
        fun = lambda s: print('\t'.join((s.header, hashlib.md5(s.seq.encode('ascii')).hexdigest())))
    # DEFAULT: Hash each header-sequence pair
    else:
        fun = lambda s: h.update('\n'.join((s.header, s.seq)).encode('ascii'))

    for seq in gen.next():
        if args.ignore_case:
            if args.all_headers or args.whole_file :
                seq.header_upper()
            if not args.all_headers:
                seq.seq_upper()
        fun(seq)

    # Print output hash for cumulative options
    if not args.each_sequence:
        print(h.hexdigest())

def complexity(args, gen):
    try:
        w = int(args.window_length)
        m = int(args.word_length)
        k = pow(int(args.alphabet_size), m)
        p = int(args.jump)
        offset = int(args.offset)
    except ValueError:
        print('All values must be integers')
        raise SystemExit

    # Calculates the variable component of a single window's score
    def varscore_generator(seq, words):
        for i in range(offset, len(seq) - w + 1, p):
            counts = defaultdict(int)
            for j in range(i, i+w):
                try:
                    counts[words[j]] += 1
                except IndexError:
                    break
            yield sum([math.log(math.factorial(x), k) for x in counts.values()])

    w_fact = math.log(math.factorial(w), k)

    seq_id = 0
    for seq in gen.next():
        seq_id += 1
        mean = 'NA'
        var = 'NA'
        if(len(seq.seq) < w + offset): pass
        elif(args.drop and args.drop in seq.seq): pass
        else:
            words = tuple(seq.seq[i:i+m] for i in range(offset, len(seq.seq) - m + 1))
            varscores = tuple(score for score in varscore_generator(seq.seq, words))
            winscores = tuple((1 / w) * (w_fact - v) for v in varscores)
            mean = sum(winscores) / len(winscores)
            try:
                var = sum([pow(mean - x, 2) for x in winscores]) / (len(varscores) - 1)
            except ZeroDivisionError:
                var = 'NA'
            try:
                col1 = seq.getvalue('gb', quiet=True)
            except:
                col1 = seq_id
            print("{},{},{}".format(col1, mean, var))

def fasta2csv(args, gen):
    w = csv.writer(sys.stdout, delimiter=args.delimiter)
    out = []
    if(args.header):
        if(args.fields):
            w.writerow(args.fields + ['seq'])
        else:
            w.writerow(['header', 'seq'])
    for seq in gen.next():
        if(args.fields):
            elements = [seq.getvalue(field) for field in args.fields]
        else:
            elements = [seq.header]
        w.writerow(tuple(elements + [seq.seq]))

def fstat(args, gen):
    stats = Stat()
    nseqs = 0
    nchars = 0
    for seq in gen.next():
        nseqs += 1
        nchars += len(seq.seq)
        stats.update_counts(seq)
    for k, v in sorted(stats.getdict().items(), key=lambda x: x[1], reverse=True):
        print("{}: {} {}".format(k, v, round(v/nchars, 4)))
    print("nseqs: {}".format(nseqs))
    print("nchars: {}".format(nchars))
    print("mean length: {}".format(round(nchars/nseqs, 4)))

def hstat(args, gen):
    ''' Writes chosen header and seq length data to csv '''
    fieldnames = list(args.fields)
    if(args.length):
        fieldnames.append('length')
    w = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    w.writeheader()
    for seq in gen.next():
        row = {}
        for field in args.fields:
            row[field] = seq.getvalue(field)
        if(args.length):
            row['length'] = len(seq.seq)
        w.writerow(row)

def idsearch(args, gen):
    ''' Print entries whose headers contain a field with a given value '''
    # TODO make mass search
    for seq in gen.next():
        if(seq.field_is(args.field, args.value)):
            seq.print()

def perm(args, gen):
    w = args.word_size
    start = args.start_offset
    end = args.end_offset
    for seq in gen.next():
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

def prettyprint(args, gen):
    ''' Print each sequence with even columns of length cwidth '''
    for seq in gen.next():
        seq.print(args.cwidth)

def qstat(args, gen):
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
    for seq in gen.next():
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

def reverse(args, gen):
    ''' Reverse each sequence '''
    for seq in gen.next():
        rseq = FSeq(seq.header, seq.seq[::-1])
        rseq.print()

def retrieve(args, gen):
    '''
    Extracts sequences matching a certain pattern
    '''
    def getset(values, filename):
        s = set()
        if filename:
            with open(filename, 'r') as f:
                s.update([line.rstrip('\n') for line in f.readlines()])
        if values:
            s.update(values)
        return(s)

    groups = getset(args.groups, args.groupfile)
    pat = set((re.compile(p) for p in getset(args.pattern, args.patternfile)))

    if(len(pat) > 1 and groups):
        print('Cannot process multiple patterns with groups', file=sys.stderr)
        raise SystemExit
    if not pat:
        print('Please provide a pattern (-p <regex>)', file=sys.stderr)

    for seq in gen.next():
        if(args.match_sequence):
            text = seq.seq
        else:
            text = seq.header
        for p in pat:
            m = re.search(p, text)
            if not m and not args.invert:
                continue
            elif groups:
                try:
                    # If at least one group defined
                    match = m.group(1)
                except:
                    # If no groups defined
                    match = m.group(0)
                if (match in groups and args.invert) or \
                   (match not in groups and not args.invert):
                    continue
            if not args.color or args.invert:
                seq.print()
                break
            # KLUDGE, TODO: Make seq and header classes
            else:
                if(args.match_sequence):
                    seq.colseq.setseq(text)
                    seq.colseq.colormatch(p)
                else:
                    seq.colheader.setseq(text)
                    seq.colheader.colormatch(p)
        if(seq.colseq.seq or seq.colheader.seq):
            seq.print()

def search(args, gen):
    '''
    Print entries whose headers contain a given pattern. Similar to `retrieve` but lighterweight.
    '''
    prog = re.compile(args.pattern)
    for seq in gen.next():
        text = seq.seq if args.seq else seq.header
        m = prog.search(text)
        if (not m and not args.invert) or (m and args.invert):
            continue
        if(args.color and not args.invert):
            if(args.seq):
                seq.colseq.setseq(text)
                seq.colseq.colormatch(prog)
            else:
                seq.colheader.setseq(text)
                seq.colheader.colormatch(prog)
        seq.print()

def simplifyheader(args, gen):
    if(hasattr(args.fields, '__upper__')):
        args.fields = (args.fields, )
    for seq in gen.next():
        values = [seq.getvalue(field) for field in args.fields]
        pairs = ['|'.join((args.fields[i], values[i])) for i in range(len(values))]
        header = '|'.join(pairs)
        FSeq(header, seq.seq).print()

def split(args, gen):
    k = int(args.nfiles)
    p = args.prefix
    seqs = []
    for seq in gen.next():
        seqs.append(seq)
    for i in range(0, k):
        begin = i * (len(seqs) // k + 1)
        end = min(len(seqs), (i+1) * (len(seqs) // k + 1))
        outfile = p + str(i) + '.fasta'
        with open(outfile, 'w+') as fo:
            for seq in (seqs[x] for x in range(begin, end)):
                fo.write(seq.get_pretty_string() + '\n')

def subseq(args, gen):
    ''' Index starting from 1 (not 0) '''
    for seq in gen.next():
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

def tounk(args, gen):
    if(args.type.lower()[0] in ['a', 'p']):
        irr = args.pir
        unk = args.punk
    elif(args.type.lower()[0] in ['n', 'd']):
        irr = args.nir
        unk = args.nunk
    if(args.lc):
        irr = ''.join(set(irr) | set(string.ascii_lowercase))
    trans = str.maketrans(irr, unk * len(irr))
    for seq in gen.next():
        seq.seq = seq.seq.translate(trans)
        seq.print()

def fsubseq(args, gen):
    '''
    Extracts many subsequences. Positional argument file
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

    for seq in gen.next():
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

def unmask(args, gen):
    ''' Converts to upcase '''
    for seq in gen.next():
        if(args.to_x):
            unmasked_seq = FSeq(seq.header, re.sub('[a-z]', 'X', seq.seq))
        else:
            unmasked_seq = FSeq(seq.header, seq.seq.upper())
        unmasked_seq.print()


# ==============
# FULL FUNCTIONS
# ==============

def sample(args, gen):
    ''' Randomly sample n entries from input file '''
    seqs = [s for s in gen.next()]
    sample_indices = random.sample(range(len(seqs)), min(len(seqs), args.n))
    [seqs[i].print() for i in sample_indices]

def sort(args, gen):
    seqs = [s for s in gen.next()]
    seqs.sort(key=lambda x: list(x.getvalue(y) for y in args.fields))
    [s.print() for s in seqs]


# =======
# EXECUTE
# =======

if __name__ == '__main__':
    gen = FSeqGenerator()
    args = parse()
    args.func(args, gen)
