#!/usr/bin/env python3

import argparse
import math
import re
import sys
import string
import copy
from collections import Counter
from collections import defaultdict
from hashlib import md5

__version__ = "1.9.0"

# ================
# Argument Parsing
# ================

class Parser:
    def __init__(self):
        self.parser = self._build_base_parser()
        self.subparsers = self._build_subparsers()
        self.usage = '<fastafile> | smof {} <options>'

    def _build_base_parser(self):
        parser = argparse.ArgumentParser(
            prog='smof',
            usage='<fastafile> | smof <subcommand> <options>',
            description='Tools for studying and manipulating fasta files')
        parser.add_argument(
            '-v', '--version',
            action='version',
            version='%(prog)s {}'.format(__version__))
        return(parser)

    def _build_subparsers(self):
        subparsers = self.parser.add_subparsers(
            metavar='[ for help on each: smof <subcommand> -h ]',
            title='subcommands')
        return(subparsers)

def parse(argv=None):

    parser = Parser()

    subcommands = [Chksum, Clean, Complexity, Fasta2csv, Grep,
                   Head, Perm, Rename, Reverse, Sample, Sniff,
                   Sort, Split, Stat, Subseq, Tail, Uniq, Wc, Winnow]
    for cmd in subcommands:
        cmd(parser)

    argv = argv if argv else sys.argv[1:]

    if not argv:
        parser.parser.print_help()
        sys.exit(0)

    if argv[0] in ['idsearch', 'retrieve', 'search', 'rmfields']:
        err("{} is deprecated, use 'smof grep'".format(argv[0]))

    args = parser.parser.parse_args(argv)

    return(args)


# =================
# CLASS DEFINITIONS
# =================

class Alphabet:
    # Including U, for selenocysteine, and * for STOP
    PROT     = set('ACDEFGHIKLMNPQRSTUVWYX*')
    PROT_UNK = set('X')
    PROT_AMB = set('BZJ')
    PROT_EXC = set('EFILQPXJZ*') # unique to protein sequences
    DNA      = set('ACGTN')
    DNA_UNK  = set('N')
    DNA_AMB  = set('RYSWKMDBHV')
    RNA      = set('ACGUN')
    RNA_UNK  = set('N')
    RNA_AMB  = set('RYSWKMDBHV')
    GAP      = set('.-_')
    STOP     = {'TAG', 'TAA', 'TGA', 'UAG', 'UAA', 'UGA'}
    START    = {'ATG', 'AUG'}

class Colors:
    OFF          = chr(27) + '[0;0m'
    RED          = chr(27) + '[0;31m'
    GREEN        = chr(27) + '[0;32m'
    YELLOW       = chr(27) + '[0;33m'
    MAGENTA      = chr(27) + '[0;35m'
    CYAN         = chr(27) + '[0;36m'
    WHITE        = chr(27) + '[0;37m'
    BLUE         = chr(27) + '[0;34m'
    BOLD_RED     = chr(27) + '[1;31m'
    BOLD_GREEN   = chr(27) + '[1;32m'
    BOLD_YELLOW  = chr(27) + '[1;33m'
    BOLD_MAGENTA = chr(27) + '[1;35m'
    BOLD_CYAN    = chr(27) + '[1;36m'
    BOLD_WHITE   = chr(27) + '[1;37m'
    BOLD_BLUE    = chr(27) + '[1;34m'

    patstr = chr(27) + '\[[0-9;]+m'
    pat = re.compile(patstr)

    COLORS = {
        'red'          : RED,
        'green'        : GREEN,
        'yellow'       : YELLOW,
        'magenta'      : MAGENTA,
        'cyan'         : CYAN,
        'white'        : WHITE,
        'blue'         : BLUE,
        'bold_red'     : BOLD_RED,
        'bold_green'   : BOLD_GREEN,
        'bold_yellow'  : BOLD_YELLOW,
        'bold_magenta' : BOLD_MAGENTA,
        'bold_cyan'    : BOLD_CYAN,
        'bold_white'   : BOLD_WHITE,
        'bold_blue'    : BOLD_BLUE
    }

class ColorAA:
    def __init__(self):
        self.group = [
            ['LVGAIP',    'aliphatic', Colors.BOLD_BLUE],
            ['FYW',       'aromatic',  Colors.BOLD_RED],
            ['SEKDRTNQH', 'polar',     Colors.BOLD_GREEN],
            ['MCU',       'thiol',     Colors.BOLD_YELLOW]
        ]
        # add lower cases
        self.group = [[l + l.lower(),g,c] for l,g,c in self.group]

    def color(a):
        for chars, group, color in self.group:
            if a in chars:
                return(Colors.OFF + a + color)
        return(a)

class ColorString:
    def __init__(self,
                 seq=None,
                 bgcolor=Colors.OFF,
                 default=Colors.BOLD_RED):
        self.bgcolor = bgcolor
        self.default = default
        self.cind = []
        if seq:
            self.append(seq)

    def append(self, thing, bg=None):
        bg = bg if bg else Colors.OFF
        if isinstance(thing, FSeq):
            thing = thing.colseq
        if isinstance(thing, ColorString):
            newcind = thing.cind
        elif isinstance(thing, str):
            newcind = []
            pos = []
            color_on = False
            for s in re.split('((?:{})+)'.format(Colors.patstr), thing):
                if s and s[0] == chr(27):
                    # if there are multiple colors, take only the last
                    col = chr(27) + s.split(chr(27))[-1]
                    # start a new color
                    if not color_on and col != self.bgcolor:
                        newcind.append([pos[-1], None, col])
                    # end old color (at previous index), start new color
                    elif color_on and col != self.bgcolor:
                        newcind[-1][1] = pos[-1]
                        newcind.append([pos[-1], None, col])
                    # end old color
                    elif color_on and col == self.bgcolor:
                        newcind[-1][1] = pos[-1]
                    color_on = bool(col != self.bgcolor)
                else:
                    last = pos[-1] if pos else 0
                    pos.append(len(s) + last)
            if newcind and not newcind[-1][1]:
                newcind[-1][1] = len(thing)
        else:
            err('ColorString can only append strings, Fseq, or ColorString objects')
        self.cind += [[b + len(self.cind), e + len(self.cind), c] for b,e,c in newcind]

    def subseq(self, a, b):
        self.cind = [x for x in self.cind if not (x[0] > b-1 or x[1]-1 < a)]
        self.cind = [[max(0, s - a), e - a, c] for s,e,c in self.cind]

    def reverse(self, length):
        self.cind = [[length - e, length - s, c] for s,e,c in self.cind]

    def colorpos(self, a, b, col=None):
        col = self.default if not col else col
        for i in reversed(range(len(self.cind))):
            start,end,color = self.cind[i]
            # prior starts in interval and extends beyond it
            if (a <= start < b) and end > b:
                # reset start site
                self.cind[i][0] = b
            # prior begins before interval and ends within it
            elif start < a and (a < end < b):
                # reset end site
                self.cind[i][1] = a
            # prior is contained within interval
            elif a < start and b > end:
                # coloror over the prior
                del self.cind[i]
            # prior contains interval
            elif a > start and b < end:
                # replace prior with colorored regions flanking interval
                del self.cind[i]
                self.cind.append([start, a, color])
                self.cind.append([b, end, color])
        self.cind.append([a, b, col])

    def print(self, seq, colwidth=None, out=sys.stdout):
        starts = {x[0]:x[2] for x in self.cind}
        ends = set([x[1] for x in self.cind])
        colored = False
        for i in range(len(seq)):
            try:
                out.write(starts[i])
                colored = True
            except:
                if i in ends:
                    colored -= 1
                    out.write(self.bgcolor)
                    colored = False
            if(colwidth and i % colwidth == 0 and i != 0):
               out.write("\n")
            out.write(seq[i])
        out.write(self.bgcolor if colored else '')
        out.write("\n")


    def colormatch(self, pattern, col=None):
        col = self.default if not col else col
        s = ''.join([x[1] for x in self.seq])
        for m in re.compile(pattern).finditer(s):
            a = m.start()
            b = m.start() + len(m.group())
            self.colorpos(a, b, col)

    def copy(self):
        new_obj = ColorString(bgcolor=self.bgcolor, default=self.default)
        new_obj.cind = self.cind
        return(new_obj)

class FileDescription:
    def __init__(self):
        self.seqs = set()
        self.headers = set()
        self.ntype = {
            'prot':0,
            'dna':0,
            'rna':0,
            'illegal':0,
            'ambiguous':0}
        self.ncase = {
            'uppercase':0,
            'lowercase':0,
            'mixedcase':0}
        self.pfeat = {
            'selenocysteine':0,
            'initial-Met':0,
            'internal-stop':0,
            'terminal-stop':0}
        from itertools import product
        self.nfeat = {''.join(c):0 for c in product('01', repeat=4)}
        self.ufeat = {
            'gapped':0,
            'unknown':0,
            'ambiguous':0}

    def add_seq(self, seq):
        '''
        Calculates properties for one sequence
        @type seq: FSeq object
        '''
        # Add md5 hashes of sequences and headers to respective sets
        self.seqs.update([md5(bytes(seq.seq, 'ascii')).digest()])
        self.headers.update([md5(bytes(seq.header, 'ascii')).digest()])

        counts = Counter(seq.seq)

        # Ungapped the sequence is required for downstream analysis
        is_gapped = self._handle_gaps(seq, counts)
        if is_gapped:
            seq.ungap()
            counts = Counter(seq.seq)

        scase = self._handle_case(counts)
        # Sum upper and lowercase
        if scase not in 'upper':
            counts = counter_caser(counts)

        s = seq.seq.upper()

        # ('prot'|'dna'|'rna'|'amb'|'bad')
        stype = self._handle_type(counts)

        if stype == 'prot':
            tstop = bool('*' == s[-1])
            self.pfeat['terminal-stop']  += tstop
            self.pfeat['internal-stop']  += (counts['*'] - tstop) > 0
            self.pfeat['selenocysteine'] += 'U' in counts
            self.pfeat['initial-Met']    += 'M' == s[0]
            self.ufeat['unknown']        += 'X' in counts
            self.ufeat['ambiguous']      += bool(Alphabet.PROT_AMB & set(counts))
        elif stype in ('dna', 'rna'):
            self.ufeat['unknown'] += bool('N' in counts)
            self.ufeat['ambiguous'] += bool(Alphabet.DNA_AMB & set(counts))

            start  = self._has_start(s)
            stop   = self._has_stop(s)
            triple = self._is_triple(s)
            sense  = self._is_sense(s)

            profile = ''.join([str(int(x)) for x in (start, stop, triple, sense)])
            self.nfeat[profile] += 1

    @classmethod
    def _has_start(cls, s):
        '''
        Tests if the first codon is the START codon, assumes uppercase
        '''
        return(s[0:3] == 'ATG')

    @classmethod
    def _has_stop(cls, s):
        '''
        Tests if the last three letters are a STOP codon, assumes uppercase
        '''
        return(s[-3:] in Alphabet.STOP)

    @classmethod
    def _is_triple(cls, s):
        '''
        Tests is the sequence is multiple of three
        '''
        return(len(s) % 3 == 0)

    @classmethod
    def _is_sense(cls, s):
        '''
        Tests if there are no internal STOP codons in the first frame
        '''
        codons = set(s[i:i+3] for i in range(0, len(s) - 3, 3))
        return(bool(not codons & Alphabet.STOP))

    def _handle_gaps(self, seq, counts):
        # Handle gaps
        if Alphabet.GAP & set(counts):
            self.ufeat['gapped'] += 1
            return(True)
        return(False)

    def _handle_case(self, counts):
        has_lower = set(string.ascii_lowercase) & set(counts)
        has_upper = set(string.ascii_uppercase) & set(counts)
        if has_lower and has_upper:
            case = 'mixedcase'
        elif has_lower:
            case = 'lowercase'
        else:
            case = 'uppercase'
        self.ncase[case] += 1
        return(case)

    def _handle_type(self, counts):
        stype = guess_type(counts)
        self.ntype[stype] += 1
        return(stype)

class FileStat:
    def __init__(self):
        self.counts = Counter()
        self.nseqs = 0
        self.lengths = []

    def add_seq(self, stat):
        if stat.counts:
            self.counts += stat.counts
        self.nseqs += 1
        self.lengths.append(stat.length)

class FSeq:
    # The translator for taking reverse complements
    revtrans = str.maketrans('acgtACGT', 'tgcaTGCA')
    # Translator of ungapping
    ungapper = str.maketrans('', '', ''.join(Alphabet.GAP))
    def __init__(self, header, seq, handle_color=False, purge_color=False):
        self.seq = seq
        self.header = header
        self.colseq = None
        self.colheader = None
        self.handle_color = handle_color
        if purge_color or handle_color:
            self._process_color(handle_color)

    def _process_color(self, handle_color=True):
        if not self.colseq:
            self.colseq = ColorString()
        if not self.colheader:
            self.colheader = ColorString()
        if bool(re.search(Colors.pat, self.seq)):
            if handle_color:
                self.colseq.append(self.seq)
            self.seq = self._clear_color(self.seq)
        if bool(re.search(Colors.pat, self.header)):
            if handle_color:
                self.colheader.append(self.header)
            self.header = self._clear_color(self.header)

    def __hash__(self):
        return(hash((self.header, self.seq)))

    def __eq__(self, other):
        return((self.header, self.seq) == (other.header, other.seq))

    def _clear_color(self, text):
        return(re.sub(Colors.pat, '', text))

    def subseq(self, a, b):
        header_suffix = '|SUBSEQ({}..{})'.format(a + 1, b)
        header = ParseHeader.firstword(self.header) + header_suffix
        sub = FSeq(header, self.seq[a:b])
        if self.colheader:
            # TODO implement header color preservation
            pass
        if self.colseq:
            sub.colseq.seq = self.colseq.seq[a:b]
        return(sub)

    def color_seq(self, *args, **kwargs):
        if not self.colseq:
            self._process_color()
        self.colseq.colorpos(*args, **kwargs)

    def color_header(self, *args, **kwargs):
        if not self.colheader:
            self._process_color()
        self.colheader.colorpos(*args, **kwargs)

    def ungap(self, suffix='|UNGAPPED'):
        self.seq = self.seq.translate(FSeq.ungapper)
        if suffix:
            self.header = self.header + suffix

    def print(self, column_width=80, color=True, out=sys.stdout):
        out.write('>')
        if self.colheader and color:
            self.colheader.print(self.header, column_width, out=out)
        else:
            out.write('%s\n' % self.header)
        for i in range(0, len(self.seq), column_width):
            if self.colseq and color:
                self.colseq.print(self.seq, column_width, out=out)
                break
            else:
                out.write('%s\n' % self.seq[i:i + column_width])

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

    def subseq(self, a, b, suffix=None):
        suffix = suffix if suffix else '|SUBSEQ({}..{})'.format(a+1,b)
        header = ParseHeader.firstword(self.header) + suffix
        newseq = FSeq(header, self.seq[a:b])
        if self.colseq:
            newseq.colseq = self.colseq.copy()
            newseq.colseq.subseq(a, b)
        return(newseq)

    def reverse(self, suffix='|REVERSE'):
        self.seq = self.seq[::-1]
        if self.handle_color:
            self.colseq.reverse(len(self.seq))
        self.header += suffix

    @classmethod
    def getrevcomp(cls, seq):
        trans = lambda s: s.translate(FSeq.revtrans)
        if isinstance(seq, str):
            return(trans(seq[::-1]))
        elif isinstance(seq, FSeq):
            newheader = ParseHeader.firstword(seq.header) + '|REVCOM'
            newseq = FSeq(newheader, trans(seq.seq[::-1]))
            if seq.colseq:
                newseq.colseq = seq.colseq
                newseq.colseq.reverse(len(seq.seq))
            if seq.colheader:
                # TODO implement this
                pass
            return(newseq)

class FSeqGenerator:
    def __init__(self, fh=sys.stdin):
        '''
        fh can be any iterable object
        '''
        self.fh = fh

    def next(self, *args, **kwargs):
        seq_list = []
        header = ''
        for line in self.fh:
            line = line.strip()
            if not line or line[0] == '#':
                continue
            if ">" == line[0]:
                if seq_list:
                    yield FSeq(header, ''.join(seq_list), *args, **kwargs)
                elif header:
                    # If no sequence but there is a header ...
                    err("Illegally empty sequence")
                seq_list = []
                header = line[1:]
            elif header:
                seq_list.append(line)
            else:
                err("First fasta line must begin with '>'")
        if header:
            if seq_list:
                yield FSeq(header, ''.join(seq_list), *args, **kwargs)
            else:
                err("Illegally empty sequence")
        else:
            err("smof could not retrieve any sequence from this file, exiting")

class Maps:
    DNA_AMB = {
        'R':'AG',
        'Y':'CT',
        'S':'GC',
        'W':'AT',
        'K':'GT',
        'M':'AC',
        'B':'CGT',
        'D':'AGT',
        'H':'ACT',
        'V':'ACG',
        'N':'ACGT'
    }

class ParseHeader:
    def firstword(h):
        return(re.sub('^(\S+).*', '\\1', h))

    def ncbi_format(h, fields):
        raise NotImplemented

    def regex_group(h, regex):
        raise NotImplemented

class SeqStat:
    def __init__(self, seq, count=True):
        self.counts = Counter(seq.seq) if count else None
        self.header = seq.header
        self.length = len(seq.seq)

    def aslist(self,
               charset=None,
               length=False,
               masked=False,
               header_fun=None,
               ignorecase=False):

        if not charset:
            charset = set(self.counts)
        line = []
        if header_fun:
            try:
                line.append(header_fun(self.header))
            except TypeError:
                err("Cannot process header: '{}'".format(self.header))
        if length:
            line.append(self.length)

        if masked:
            line.append(sum_lower(self.counts))

        if ignorecase:
            charset = set(''.join(charset).upper())
            self.counts = counter_caser(self.counts)

        line += [self.counts[c] for c in sorted(charset)]

        return(line)

    @classmethod
    def getheader(cls, charset, length=False, masked=False, ignorecase=False):
        header = ['seqid']

        if ignorecase:
            charset = set(''.join(charset).upper())

        if length:
            header += ['length']

        if masked:
            header += ['masked']

        header += list(sorted(charset))

        return(header)

class StatFun:
    @classmethod
    def N50(cls, x, issorted=False):
        x = sorted(x) if not issorted else x
        N = sum(x)
        total = 0
        for i in range(len(x)-1, -1, -1):
            total += x[i]
            if total > N / 2:
                return(x[i])

    @classmethod
    def mean(cls, x):
        return(sum(x) / len(x))

    @classmethod
    def median(cls, x, issorted=False):
        return(cls.quantile(x, 0.5, issorted=issorted))

    @classmethod
    def sd(cls, x):
        if(len(x) < 2):
            return('NA')
        else:
            mean = sum(x) / len(x)
            stdev = (sum((y - mean) ** 2 for y in x) / (len(x) - 1)) ** 0.5
            return(stdev)
        return(sd)

    @classmethod
    def quantile(cls, x, q, issorted=False):
        '''
        Calculates quantile as the weighted average between indices
        '''
        # Die if out of bounds
        if not (0 <= q <= 1):
            err('quantile must be between 0 and 1')

        # Ensure the vector is sorted
        x = sorted(x) if not issorted else x

        # Return max or min for q = 1 or 0
        if q == 1:
            return(x[-1])
        elif q == 0:
            return(x[0])

        v = (len(x) - 1) * q
        r = v % 1
        i = math.floor(v)
        quantile = x[i] * (1-r) + x[i+1] * r
        return(quantile)

    @classmethod
    def summary(cls, x):
        x = sorted(x)
        out = {
            'min':x[0],
            'max':x[-1],
            '1st_qu':cls.quantile(x, 0.25, issorted=True),
            'median':cls.quantile(x, 0.50, issorted=True),
            '3rd_qu':cls.quantile(x, 0.75, issorted=True),
            'mean':cls.mean(x),
            'sd':cls.sd(x),
            'N50':cls.N50(x, issorted=True)
        }
        return(out)


# =================
# UTILITY FUNCTIONS
# =================

def csvrowGenerator(filename):
    import csv
    try:
        f = open(filename, 'r')
        dialect = csv.Sniffer().sniff(f.read(1024))
        f.seek(0)
        reader = csv.reader(f, dialect)
    except IOError as e:
        print(e, file=sys.stderr)
        sys.exit(0)
    for row in reader:
        yield row

def counter_caser(counter, lower=False):
    '''
    Sums cases in Collections.Counter object
    '''
    if lower:
        out = counter + Counter({k.lower():v for k,v in counter.items() if k.isupper()})
        out = out - Counter({k:v for k,v in counter.items() if k.isupper()})
    else:
        out = counter + Counter({k.upper():v for k,v in counter.items() if k.islower()})
        out = out - Counter({k:v for k,v in counter.items() if k.islower()})
    return(out)

def sum_lower(counter):
    lc = [v for k,v in counter.items() if k in string.ascii_lowercase]
    return(sum(lc))

def guess_type(counts):
    '''
    Predict sequence type from character counts (dna|rna|prot|ambiguous|illegal)
    '''
    if(isinstance(counts, str)):
        counts = Counter(counts)
    elif(isinstance(counts, FSeq)):
        counts = Counter(counts.seq)

    # Convert all to upper case
    counts = counter_caser(counts)
    # Remove gaps from Counter
    counts = Counter({k:n for k,n in counts.items() if k not in Alphabet.GAP})

    # If all chars are in ACGT
    if set(counts) <= Alphabet.DNA:
        stype = 'ambiguous' if sum(counts.values()) < 3 else 'dna'
    # If all chars in ACGU
    elif set(counts) <= Alphabet.RNA:
        stype = 'rna'
    # If has any chars unique to proteins (EFILQPXJZ*)
    elif set(counts) & Alphabet.PROT_EXC:
        if set(counts) <= Alphabet.PROT | Alphabet.PROT_AMB:
            stype = 'prot'
        else:
            stype = 'illegal'
    # If all the residues could be aa, DNA, or RNA
    elif set(counts) <= (Alphabet.PROT | Alphabet.PROT_AMB):
        # If more than 80% look like nucleic acids, set 'amb_nucl'
        if (sum([counts[x] for x in 'ACGTUN' if x in counts]) / sum(counts.values())) > 0.8:
            if 'U' in counts:
                stype = 'illegal' if 'T' in counts else 'rna'
            else:
                stype = 'dna'
        # Otherwise set as ambibuous
        else:
            stype = 'ambiguous'
    # If none of these match, something is horribly wrong with your
    # sequence
    else:
        stype = 'illegal'
    return(stype)

def headtailtrunk(seq, first=None, last=None):
    '''
    This function is used by the Head and Tail classes to portray partial
    sections of sequences.
    '''
    outseq = FSeq(seq.header, seq.seq)
    if first and last:
        if first + last < len(seq.seq):
            outseq.header = ParseHeader.firstword(seq.header) + \
                    '|TRUNCATED:first-{}_last-{}'.format(first, last)
            outseq.seq = '{}{}{}'.format(
                seq.seq[0:first],
                '...',
                seq.seq[-last:]
            )
    elif first:
        outseq.header = ParseHeader.firstword(seq.header) + \
                '|TRUNCATED:first-{}'.format(first)
        outseq.seq = seq.seq[0:first]
    elif last:
        outseq.header = ParseHeader.firstword(seq.header) + \
                '|TRUNCATED:last-{}'.format(last)
        outseq.seq = seq.seq[-last:]
    elif first == 0 and last == 0:
        err("Illegal empty sequence, dying ...")
    return(outseq)

def ascii_histchar(dif, chars=' .~*O'):
    if dif <= 0:
        return(chars[0])
    elif dif < 0.25:
        return(chars[1])
    elif dif < 0.5:
        return(chars[2])
    elif dif < 0.75:
        return(chars[3])
    else:
        return(chars[4])

def counting_number(i):
    i = int(i)
    if i < 1:
         raise argparse.ArgumentTypeError("%s is an invalid counting number" % i)
    return i

def positive_int(i):
    i = int(i)
    if i < 0:
         raise argparse.ArgumentTypeError("%s is an invalid positive number" % i)
    return i

def err(msg):
    print(msg, file=sys.stderr)
    sys.exit(1)

def ambiguous2perl(pattern):
    perlpat = []
    in_bracket = False
    escaped = False
    for c in pattern:
        amb = c in Maps.DNA_AMB
        if c == '\\':
            escaped = True
            continue
        elif escaped:
            c = c if amb else '\\' + c
            escaped = False
        elif amb:
            v = Maps.DNA_AMB[c]
            c = v if in_bracket else '[%s]' % v
        elif c == '[':
            in_bracket = True
        elif c == ']':
            in_bracket = False
        perlpat.append(c)
    return(''.join(perlpat))

# ====================
# ONE-BY-ONE FUNCTIONS
# ====================

class Subcommand:
    def __init__(self, parser_obj, force_color=False):
        self.force_color = force_color
        self.func = self.write
        self.usage = parser_obj.usage
        self.subparsers = parser_obj.subparsers
        self._parse()

    def _parse(self):
        raise NotImplemented

    def generator(self, args, gen):
        raise NotImplemented

    def write(self, args, gen, out=sys.stdout):
        for output in self.generator(args, gen):
            if(isinstance(output, FSeq)):
                color = sys.stdout.isatty() or self.force_color
                output.print(color=color, out=out)
            else:
                out.write('%s\n' % output)

class Chksum(Subcommand):
    def _parse(self):
        cmd_name = 'chksum'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="calculate an md5 checksum for the input sequences"
        )
        parser.add_argument(
            '-i', '--ignore-case',
            help='convert all to uppercase, before hashing',
            action='store_true',
            default=False
        )
        method = parser.add_mutually_exclusive_group(required=False)
        method.add_argument(
            '-w', '--whole-file',
            help='calculate single md5sum for all headers and sequences (default action)',
            action='store_true',
            default=False
        )
        method.add_argument(
            '-q', '--each-sequence',
            help='calculate md5sum for each sequence, write as TAB delimited list',
            action='store_true',
            default=False
        )
        method.add_argument(
            '-s', '--all-sequences',
            help='calculate single md5sum for all sequences',
            action='store_true',
            default=False
        )
        method.add_argument(
            '-d', '--all-headers',
            help='calculate single md5sum for all headers',
            action='store_true',
            default=False
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        md5hash = md5()
        # Hash the sequences only (in input order)
        if(args.all_sequences):
            fun = lambda s,h: md5hash.update(s)
        # Hash the headers only (in input order)
        elif(args.all_headers):
            fun = lambda s,h: md5hash.update(h)
        # DEFAULT: Hash each header-sequence pair
        else:
            fun = lambda s,h: md5hash.update(h + b'\n' + s)

        for seq in gen.next():
            if args.ignore_case:
                if args.all_headers or args.whole_file :
                    seq.header_upper()
                if not args.all_headers:
                    seq.seq_upper()
            s = seq.seq.encode('ascii')
            h = seq.header.encode('ascii')
            # Write <header>\t<sequence hash> for each sequence
            if(args.each_sequence):
                yield h.decode() + '\t' + md5(s).hexdigest()
            else:
                fun(s,h)

        # Print output hash for cumulative options
        if not args.each_sequence:
            yield md5hash.hexdigest()

class Clean(Subcommand):
    def _parse(self):
        cmd_name = 'clean'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="masks things (optionally), pretty prints")
        parser.add_argument(
            '-t', '--type',
            metavar='n|p',
            help='sequence type'
        )
        parser.add_argument(
            '-u', '--toupper',
            help="convert to uppercase",
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-l', '--tolower',
            help="convert to lowercase",
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-x', '--toseq',
            help="removes all nonletter characters (gaps, stops, etc.)",
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-r', '--mask-irregular',
            help="converts irregular letters to unknown",
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-m', '--mask-lowercase',
            help='convert lower-case to unknown',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '--nir',
            help='nucleotide irregulars [default=(not ACGT)]',
            metavar='STR',
            default=''.join(set(string.ascii_letters) - set('ACGTNacgtn'))
        )
        parser.add_argument(
            '--pir',
            metavar='STR',
            help='protein irregulars [default=BJOUZbjouz]',
            default='BJOUXZbjouxz'
        )
        parser.add_argument(
            '--nunk',
            metavar='CHAR',
            help='nucleotide unknown character (default=N)',
            default='N'
        )
        parser.add_argument(
            '--punk',
            metavar='CHAR',
            help='protein unknown character (default=X)',
            default='X'
        )
        parser.set_defaults(func=self.func)

    def _process_args(self, args):
        if (args.mask_lowercase or args.mask_irregular) and not args.type:
            err('Please provide sequence type (--type)')

        if args.tolower and args.toupper:
            err('Err, you want me to convert to lower AND upper?')

    def generator(self, args, gen):

        # Catches illegal combinations of arguments
        self._process_args(args)

        # Make translator, this is only possible if type is given
        trans = None
        if args.type:
            if args.type.lower()[0] in ['a', 'p']:
                irr = args.pir
                unk = args.punk
            elif args.type.lower()[0] in ['n', 'd']:
                irr = args.nir
                unk = args.nunk
            else:
                err('Type not recognized')

            a = ''
            # Get irregular characters
            if args.mask_irregular:
                a += irr

            # convert lowercase to unknown
            if args.mask_lowercase:
                a += string.ascii_lowercase

            a = ''.join(set(a))

            b = unk * len(a)
            if irr:
                trans = str.maketrans(a, b)

        for seq in gen.next(purge_color=True):
            # WARNING: order is important here, don't swap thoughtlesly
            # Remove all nonletters or wanted, otherwise just remove space
            if args.toseq:
                seq.seq = re.sub('[^A-Za-z]', '', seq.seq)
            else:
                seq.seq = re.sub('[^\S]', '', seq.seq)

            # Irregular or lowercase to unknown
            if trans:
                seq.seq = seq.seq.translate(trans)

            # Change everything to desired case
            if args.toupper:
                seq.seq = seq.seq.upper()
            elif args.tolower:
                seq.seq = seq.seq.lower()

            yield seq

class Complexity(Subcommand):
    def _parse(self):
        cmd_name = 'complexity'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help='calculates linguistic complexity'
        )
        parser.add_argument(
            '-k', '--alphabet-size',
            help='number of letters in the alphabet (4 for DNA, 20 for proteins)',
            type=counting_number,
            metavar='INT',
            default=4
        )
        parser.add_argument(
            '-w', '--window-length',
            help='window length (if provided, output will average of window complexities)',
            type=counting_number,
            metavar='INT',
            default=100
        )
        parser.add_argument(
            '-m', '--word-length',
            help='length of each word',
            type=counting_number,
            metavar='INT',
            default=1
        )
        parser.add_argument(
            '-j', '--jump',
            help='distance between adjacent windows',
            type=counting_number,
            metavar='INT',
            default=1
        )
        parser.add_argument(
            '-o', '--offset',
            help='index of start pocounting_number',
            type=positive_int,
            metavar='INT',
            default=0
        )
        parser.add_argument(
            '-d', '--drop',
            help="drop sequence if contains this character (e.g. 'X' or 'N')",
            metavar='CHAR',
            default=None
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        try:
            w = int(args.window_length)
            m = int(args.word_length)
            k = pow(int(args.alphabet_size), m)
            p = int(args.jump)
            offset = int(args.offset)
        except ValueError:
            err('All values must be integers')

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

        for seq in gen.next():
            mean = 'NA'
            var = 'NA'
            seqid = ParseHeader.firstword(seq.header)
            if len(seq.seq) >= w + offset:
                words = tuple(seq.seq[i:i+m] for i in range(offset, len(seq.seq) - m + 1))
                varscores = tuple(score for score in varscore_generator(seq.seq, words))
                winscores = tuple((1 / w) * (w_fact - v) for v in varscores)
                mean = sum(winscores) / len(winscores)
                try:
                    var = sum([pow(mean - x, 2) for x in winscores]) / (len(varscores) - 1)
                except ZeroDivisionError:
                    var = 'NA'
            elif args.drop and args.drop in seq.seq:
                continue

            mean = mean if mean == 'NA' else '{:.5f}'.format(mean)
            var = var if var == 'NA' else '{:.5e}'.format(var)

            yield '\t'.join((seqid, mean, var))

class Fasta2csv(Subcommand):
    def _parse(self):
        cmd_name = 'fasta2csv'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="converts a fasta file to 2-column csv")
        parser.add_argument(
            '-d', '--delimiter',
            help="set delimiter (TAB by default)",
            metavar='CHAR',
            default='\t')
        parser.add_argument(
            '-r', '--header',
            help='write header (default=False)',
            action='store_true',
            default=False)
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        if(args.header):
            yield 'seqid{}seq'.format(args.delimiter)

        for seq in gen.next():
            seqid = ParseHeader.firstword(seq.header)
            yield '{}{}{}'.format(seqid, args.delimiter, seq.seq)

class Perm(Subcommand):
    def _parse(self):
        cmd_name = 'perm'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="randomly order sequence"
        )
        parser.add_argument(
            '-w', '--word-size',
            help='size of each word (default=1)',
            type=counting_number,
            metavar='INT',
            default=1
        )
        parser.add_argument(
            '-s', '--start-offset',
            help='number of letters to ignore at beginning (default=0)',
            type=positive_int,
            metavar='INT',
            default=0
        )
        parser.add_argument(
            '-e', '--end-offset',
            help='number of letters to ignore at end (default=0)',
            type=positive_int,
            metavar='INT',
            default=0
        )
        parser.add_argument(
            '--seed',
            help='set random seed (for reproducibility/debugging)',
            type=counting_number
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        import random
        if args.seed:
            random.seed(args.seed)
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
            header='{}|PERMUTATION:start={};end={};word_size={}'.format(
                ParseHeader.firstword(seq.header), start, end, w)
            yield FSeq(header, out)

class Reverse(Subcommand):
    def _parse(self):
        cmd_name = 'reverse'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="reverse each sequence (NOT reverse complement)")
        parser.add_argument(
            '-S', '--preserve-color',
            help='Preserve incoming color',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-Y', '--force-color',
            help='print in color even to non-tty (DANGEROUS)',
            action='store_true',
            default=False
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        ''' Reverse each sequence '''
        self.force_color = args.force_color
        for seq in gen.next(handle_color=args.preserve_color):
            seq.reverse()
            yield seq

class Sniff(Subcommand):
    def _parse(self):
        cmd_name = 'sniff'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="extract info about the sequence"
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        seqsum = FileDescription()
        for seq in gen.next():
            seqsum.add_seq(seq)
        yield seqsum

    def write(self, args, gen, out=sys.stdout):
        seqsum = next(self.generator(args, gen))
        nseqs = sum(seqsum.ncase.values())

        has_degen_headers = nseqs != len(seqsum.headers)
        has_degen_seqs = nseqs != len(seqsum.seqs)
        has_bad = seqsum.ntype['illegal'] != 0

        if has_degen_seqs:
            out.write("{} uniq sequences ({} total)\n".format(len(seqsum.seqs), nseqs))
        else:
            out.write("Total sequences: {}\n".format(nseqs))

        if has_degen_headers:
            out.write("WARNING: headers are not unique ({}/{})\n".format(len(seqsum.headers), nseqs))
        if has_bad:
            out.write("WARNING: illegal characters found\n")

        def write_dict(d, name, N):
            uniq = [[k,v] for  k,v in d.items() if v != 0]
            if len(uniq) == 1:
                out.write("All {}\n".format(uniq[0][0]))
            else:
                out.write("{}:\n".format(name))
                for k,v in sorted(uniq, key=lambda x: -x[1]):
                    out.write("  {:<20} {:<10} {:>7.4%}\n".format(k + ':', v, v/N))

        write_dict(seqsum.ntype, 'Sequence types', nseqs)
        write_dict(seqsum.ncase, 'Sequences cases', nseqs)

        nnucl = sum([v for k,v in seqsum.ntype.items() if k in {'dna', 'rna'}])
        nprot = sum([v for k,v in seqsum.ntype.items() if k == 'prot'])

        def write_feat(d, text, N, drop=False):
            if not N:
                return
            out.write('%s\n' % text)
            for k,v in sorted(list(d.items()), key=lambda x: -x[1]):
                if (drop and v != 0) or not drop:
                    out.write("  {:<20} {:<10} {:>7.4%}\n".format(k + ':', v, v/N))

        write_feat(seqsum.nfeat, "Nucleotide Features", nnucl, drop=True)
        write_feat(seqsum.pfeat, "Protein Features:", nprot)
        write_feat(seqsum.ufeat, "Universal Features:", nseqs)

class Stat(Subcommand):
    def _parse(self):
        cmd_name = 'stat'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="calculate sequence statistics")
        parser.add_argument(
            '-d', '--delimiter',
            help='output delimiter'
        )
        parser.add_argument(
            '-q', '--byseq',
            help='write a line for each sequence',
            default=False,
            action='store_true'
        )
        parser.add_argument(
            '-I', '--case-sensitive',
            help='match case',
            default=False,
            action='store_true'
        )
        parser.add_argument(
            '-m', '--count-lower',
            help='count the number of lowercase characters',
            default=False,
            action='store_true'
        )
        parser.add_argument(
            '-c', '--counts',
            help='write counts of all characters',
            default=False,
            action='store_true'
        )
        parser.add_argument(
            '-t', '--type',
            help='guess sequence type',
            default=False,
            action='store_true'
        )
        parser.add_argument(
            '-l', '--length',
            help='write sequence length',
            default=False,
            action='store_true'
        )
        parser.add_argument(
            '-p', '--proportion',
            help='write proportion of each character',
            default=False,
            action='store_true'
        )
        parser.add_argument(
            '-C', '--aa-profile',
            help='display protein profile',
            default=False,
            action='store_true'
        )
        parser.add_argument(
            '-g', '--hist',
            help='write ascii histogram of sequence lengths',
            default=False,
            action='store_true'
        )
        parser.add_argument(
            '-G', '--log-hist',
            help='write ascii histogram of sequence log2 lengths',
            default=False,
            action='store_true'
        )
        parser.set_defaults(func=self.func)

    def _process_args(self, args):
        # If no output options are specified, do length stats
        if not any((args.counts, args.type, args.length, args.proportion, args.count_lower)):
            args.length = True
        return(args)

    def _get_length_lines(self, args, g):
        lines = []
        total = sum(g.lengths)
        N = len(g.lengths)
        if N > 1:
            s = StatFun.summary(g.lengths)

            # Yield total number of sequences
            lines.append("{:10s} {}".format('nseq:', len(g.lengths)))

            # lines.append totla number of letters
            lines.append("{:10s} {}".format('nchars:', sum(g.lengths)))

            # lines.append five number summary of sequence lengths
            fivesum = [round(s[x]) for x in ('min','1st_qu','median','3rd_qu','max')]
            fivesum_str = "{:10s} {} {} {} {} {}"
            lines.append(fivesum_str.format('5sum:', *fivesum))

            # lines.append mean and standard deviation
            meansd_str="{:10s} {:d} ({:d})"
            lines.append(meansd_str.format('mean(sd):', round(s['mean']), round(s['sd'])))

            # lines.append N50
            lines.append("{:10s} {}".format('N50:', s['N50']))
        else:
            lstr = ', '.join([str(x) for x in sorted(g.lengths)])
            lines.append("nchars: {}".format(lstr))
        return(lines)

    def _get_hist_lines(self, args, g, title=None, height=10, width=60, log=False):
        lines = []
        try:
            import numpy
        except ImportError:
            err('Please install numpy (needed for histograms)')

        if title:
            lines.append('')
            lines.append(title)

        if log:
            lengths = [math.log(x, 2) for x in g.lengths]
        else:
            lengths = g.lengths

        y = numpy.histogram(lengths, bins=width)[0]
        y = [height * x / max(y) for x in y]

        for row in reversed(range(height)):
            out = ''.join([ascii_histchar(h - row) for h in y])
            lines.append('|{}|'.format(out))
        return(lines)

    def _get_aaprofile_lines(self, args, g, title=None, height=10):
        lines = []
        if title:
            lines.append('')
            lines.append(title)

        colorAA = ColorAA()
        aacols = []
        for chars, group, color in colorAA.group:
            for c in chars:
                if not args.case_sensitive and c.islower():
                    continue
                cheight = height * g.counts[c] / max(g.counts.values())
                aacols.append([c, cheight, color])
        # Draw histogram
        for row in reversed(range(height)):
            out = ''.join([c + ascii_histchar(y - row, chars=" .:'|") for l,y,c in aacols])
            out = '{}{}'.format(out, Colors.OFF)
            lines.append(out)
        names = ''.join([l for l,y,c in aacols])
        lines.append(names + Colors.OFF)
        return(lines)

    def _get_count_lines(self, args, g):
        lines = []
        lower = sum_lower(g.counts) if args.count_lower else None
        if not args.case_sensitive:
            g.counts = counter_caser(g.counts)

        if args.type:
            lines.append(guess_type(g.counts))

        N = sum(g.lengths)
        slen = str(len(str(max(g.counts.values()))) + 2)
        count_iter = sorted(g.counts.items(), key=lambda x: -x[1])
        if args.counts ^ args.proportion:
            for k,v in count_iter:
                val = v/N if args.proportion else v
                if args.counts:
                    exp = "{}{:>%sd}" % slen
                else:
                    exp = "{}{:>11.5%}"
                lines.append(exp.format(k,val))
        elif args.counts and args.proportion:
            for k,v in count_iter:
                outstr = "{}{:>" + slen + "d}{:>11.5%}"
                lines.append(outstr.format(k,v,v/N))

        if args.count_lower:
            lines.append("{:10s} {} ({:.1%})".format('lower:', lower, lower/N))
        return(lines)

    def _byfile(self, args, gen):
        g = FileStat()
        # Do I need to count the characters? (much faster if I don't)
        need_count = any((args.counts, args.proportion, args.count_lower, args.type, args.aa_profile))
        for seq in gen.next():
            g.add_seq(SeqStat(seq, count=need_count))

        if need_count:
            lines = self._get_count_lines(args, g)
            yield '\n'.join(lines)

        if args.length:
            lines = self._get_length_lines(args, g)
            yield '\n'.join(lines)

        if args.hist:
            if args.log_hist:
                lines = self._get_hist_lines(args, g, title='Flat histogram')
            else:
                lines = self._get_hist_lines(args, g)
            yield '\n'.join(lines)

        if args.log_hist:
            if args.hist:
                lines = self._get_hist_lines(args, g, title='Log2 histogram', log=True)
            else:
                lines = self._get_hist_lines(args, g, log=True)
            yield '\n'.join(lines)

        if args.aa_profile:
            if args.hist or args.log_hist:
                lines = self._get_aaprofile_lines(args, g, title='AA profile')
            else:
                lines = self._get_aaprofile_lines(args, g)
            yield '\n'.join(lines)

    def _byseq(self, args, gen):
        seqlist = []
        charset = set()
        if args.length and not (args.counts or args.proportion):
            delimiter = args.delimiter if args.delimiter else '\t'
            for seq in gen.next():
                seqid = ParseHeader.firstword(seq.header)
                yield("{}{}{}".format(seqid, delimiter, len(seq.seq)))
        else:
            for seq in gen.next():
                seqstat = SeqStat(seq)
                seqlist.append(seqstat)
                charset.update(seqstat.counts)

            delimiter = args.delimiter if args.delimiter else ','
            joiner = lambda s,d: '{}'.format(d).join([str(x) for x in s])

            ignorecase = not args.case_sensitive
            kwargs = {'masked':args.count_lower,
                    'length':args.length,
                    'ignorecase':ignorecase}

            yield joiner(SeqStat.getheader(charset, **kwargs), delimiter)

            for q in seqlist:
                line = q.aslist(charset=charset,
                                header_fun=ParseHeader.firstword,
                                **kwargs)
                yield(joiner(line, delimiter))

    def generator(self, args, gen):
        args = self._process_args(args)
        if args.byseq:
            g = self._byseq(args, gen)
        else:
            g = self._byfile(args, gen)
        for item in g:
            yield(item)

class Subseq(Subcommand):
    def _parse(self):
        cmd_name = 'subseq'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="extract subsequence from each entry (revcomp if a<b)"
        )
        parser.add_argument(
            '-c', '--color',
            help='color subsequence (do not extract)',
            choices=Colors.COLORS.keys(),
            metavar='STR',
            default=None
        )
        parser.add_argument(
            '-f', '--gff',
            metavar='FILE',
            help='get bounds from this gff3 file',
            type=argparse.FileType('r')
        )
        parser.add_argument(
            '-k', '--keep',
            help='With --gff, keep sequences with no matches',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-b', '--bounds',
            metavar='N',
            help="from and to values (indexed from 1)",
            nargs=2,
            type=counting_number
        )
        parser.add_argument(
            '-Y', '--force-color',
            help='print in color even to non-tty (DANGEROUS)',
            action='store_true',
            default=False
        )
        parser.set_defaults(func=self.func)

    def _subseq(self, seq, a, b, color=None):
        start,end = sorted([a,b])
        end = min(end, len(seq.seq))

        # Check boundaries
        if start > len(seq.seq):
            err('Start position must be less than seq length')

        if color:
            c = Colors.COLORS[color]
            seq.colseq.colorpos(start-1, end, c)
            return(seq)
        else:
            outseq = seq.subseq(start-1, end)
            if (a > b) and guess_type(seq.seq) == 'dna':
                outseq = FSeq.getrevcomp(outseq)
            return(outseq)

    def _gff_generator(self, args, gen):
        subseqs = defaultdict(list)
        for line in args.gff:
            row = line.split('\t')
            try:
                subseqs[row[0]].append({'start':int(row[3]),
                                        'end':int(row[4]),
                                        'strand':row[6]})
            except IndexError:
                err('Improper gff3 file')
            except ValueError:
                err('gff bounds must be integers')

        for seq in gen.next(handle_color=True):
            seqid = ParseHeader.firstword(seq.header)
            try:
                if seqid not in subseqs.keys():
                    raise KeyError
            except KeyError:
                if args.keep:
                    yield seq
                continue

            if args.color:
                for s in subseqs[seqid]:
                    seq = self._subseq(seq, s['start'], s['end'], args.color)
                yield seq
            else:
                for s in subseqs[seqid]:
                    yield self._subseq(seq, s['start'], s['end'])

    def _bound_generator(self, args, gen):
        for seq in gen.next(handle_color=True):
            yield self._subseq(seq, *args.bounds, color=args.color)

    def generator(self, args, gen):
        self.force_color = args.force_color
        if args.gff:
            sgen = self._gff_generator
        else:
            sgen = self._bound_generator

        for item in sgen(args, gen):
            yield item

class Winnow(Subcommand):
    def _parse(self):
        cmd_name = 'winnow'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="remove sequences which meet given conditions"
        )
        parser.add_argument(
            '-c', '--contain',
            metavar='STR',
            help="remove if contains any of these characters"
        )
        parser.add_argument(
            '-C', '--not-contain',
            metavar='STR',
            help="remove if not contains any of these characters"
        )
        parser.add_argument(
            '-s', '--shorter-than',
            help="remove if sequence is shorter than i",
            type=counting_number,
            metavar='INT'
        )
        parser.add_argument(
            '-S', '--longer-than',
            help="remove if sequence is longer than i",
            type=counting_number,
            metavar='INT'
        )
        parser.add_argument(
            '-v', '--invert',
            help="invert selection",
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-p', '--composition',
            metavar='EXPR',
            help="composition (e.g. -p 'GC < 0.5')"
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        tests = []
        if args.contain:
            tests.append(lambda s, v=args.contain: bool(set(v) & set(s)))
        if args.not_contain:
            tests.append(lambda s, v=args.not_contain: bool(set(v) & set(s)))
        if args.longer_than:
            tests.append(lambda s, v=args.longer_than: len(s) > v)
        if args.shorter_than:
            tests.append(lambda s, v=args.shorter_than: len(s) < v)
        if args.composition:
            try:
                ch,sign,per = args.composition.split()
            except ValueError:
                err('Composition argument have 3 values')
            legal_signs = ('<', '<=', '>=', '>', '==')
            if not sign in legal_signs:
                err("Middle term must be a comparison symbol ('<', '<=', '>=', '>', '==')")
            try:
                per = float(per)
            except ValueError:
                err("Third value must be a float")
            if not 0 <= per <= 1:
                err("Third value must be between 0 and 1")
            ch = set(str(ch))

            def evaluate(s):
                c = Counter(s)
                p = sum([c[x] for x in ch]) / len(s)
                return(eval('p {} {}'.format(sign, per)))

            tests.append(evaluate)


        for seq in gen.next():
            reject = [x(seq.seq) for x in tests]
            if (any(reject) and args.invert) or (not any(reject) and not args.invert):
                yield seq


# ==============
# FULL FUNCTIONS
# ==============

class Sample(Subcommand):
    def _parse(self):
        cmd_name = 'sample'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="randomly select entries from fasta file")
        parser.add_argument(
            'n',
            help="sample size",
            type=counting_number,
            nargs='?',
            default=1)
        parser.add_argument(
            '--seed',
            help='set random seed (for reproducibility/debugging)',
            type=counting_number
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        ''' Randomly sample n entries from input file '''
        import random
        if args.seed:
            random.seed(args.seed)
        seqs = [s for s in gen.next()]
        sample_indices = random.sample(range(len(seqs)), min(len(seqs), args.n))
        for i in sample_indices:
            yield seqs[i]

class Sort(Subcommand):
    def _parse(self):
        cmd_name = 'sort'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="sort sequences")
        parser.add_argument(
            '-x', '--regex',
            help="sort by single regex capture")
        parser.add_argument(
            '-r', '--reverse',
            help="reverse sort",
            action='store_true',
            default=False)
        parser.add_argument(
            '-n', '--numeric',
            help="numeric sort",
            action='store_true',
            default=False)
        parser.add_argument(
            '-l', '--length',
            help='sort by sequence length',
            action='store_true',
            default=False)
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        seqs = [s for s in gen.next()]

        if args.numeric and not args.regex:
            err('--numeric does nothing unless with --regex')

        # Set type of order determining variable
        if args.numeric:
            def typer(x):
                try:
                    return(float(x))
                except ValueError:
                    err("'{}' cannot be numerically sorted".format(x))
        else:
            def typer(x):
                return(x)

        # Set search term
        if args.regex:
            r = re.compile(args.regex)
            def sortterm(x):
                try:
                    capture = re.search(r, x.header).groups()[0]
                    return(typer(capture))
                except AttributeError:
                    err("No match for regex '{}'".format(args.regex))
                except IndexError:
                    err("Nothing was captured in regex '{}'".format(args.regex))
        elif args.length:
            def sortterm(x):
                return(len(x.seq))
        else:
            def sortterm(x):
                return(x.header)

        seqs.sort(key=lambda x: sortterm(x), reverse=args.reverse)

        for s in seqs:
            yield s

class Split(Subcommand):
    def _parse(self):
        cmd_name = 'split'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help='split a multifasta file into k smaller files'
        )
        parser.add_argument(
            '-n', '--nfiles',
            help='number of output files'
        )
        parser.add_argument(
            '-p', '--prefix',
            help='prefix for output files',
            default='xxx'
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        for s in gen.next():
            yield s

    def write(self, args, gen, out=None):
        k = int(args.nfiles)
        p = args.prefix
        seqs = list(self.generator(args, gen))
        for i in range(0, k):
            begin = i * (len(seqs) // k + 1)
            end = min(len(seqs), (i+1) * (len(seqs) // k + 1))
            outfile = p + str(i) + '.fasta'
            with open(outfile, 'w+') as fo:
                for seq in (seqs[x] for x in range(begin, end)):
                    fo.write(seq.get_pretty_string() + '\n')


# ==============
# UNIX EMULATORS
# ==============

class Head(Subcommand):
    def _parse(self):
        cmd_name = 'head'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="Write first N sequences (default=1)"
        )
        parser.add_argument(
            '-n', '--nseqs',
            help='print N sequences',
            metavar='N',
            type=counting_number,
            default=1
        )
        parser.add_argument(
            '-f', '--first',
            help='print first K letters of each sequence',
            metavar='K',
            type=counting_number
        )
        parser.add_argument(
            '-l', '--last',
            help='print last K letters of each sequence',
            metavar='K',
            type=counting_number
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        i = 1
        for seq in gen.next():
            yield headtailtrunk(seq, args.first, args.last)
            if i == args.nseqs:
                break
            i += 1

class Grep(Subcommand):
    def _parse(self):
        cmd_name = 'grep'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="roughly emulates the UNIX grep command"
        )
        parser.add_argument(
            'patterns',
            metavar='PATTERNS',
            help='patterns to match',
            nargs='*'
        )
        parser.add_argument(
            '-q', '--match-sequence',
            help='match sequence rather than header',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-f', '--file',
            metavar='FILE',
            type=argparse.FileType('r'),
            help='obtain patterns from FILE, one per line'
        )
        parser.add_argument(
            '-w', '--wrap',
            metavar='REG',
            help='a regular expression to capture PATTERNS'
        )
        parser.add_argument(
            '-P', '--perl-regexp',
            help='treat PATTERNS as perl regular expressions',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-B', '--ambiguous-nucl',
            help='parse extended nucleotide alphabet',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-I', '--case-sensitive',
            help='DO NOT ignore case distinctions (ignore by default)',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-v', '--invert-match',
            help='print non-matching entries',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-c', '--count',
            help='print number of entries with matches',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-m', '--count-matches',
            help='print number of non-overlapping matches',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-b', '--both-strands',
            help='search both strands',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-r', '--reverse-only',
            help='only search the reverse strand',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-y', '--no-color',
            help='do not print in color',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-Y', '--force-color',
            help='print in color even to non-tty (DANGEROUS)',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-S', '--preserve-color',
            help='Preserve incoming color',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '--gff',
            help='output matches in gff format',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '--gff-type',
            help='name of searched feature',
            metavar='STR',
            default='regex_match'
        )
        parser.set_defaults(func=self.func)

    def _process_arguments(self, args):
        # Stop if there are any incompatible options
        if args.count_matches and args.invert_match:
            err('--count-matches argument is incompatible with --invert-matches')

        if not args.match_sequence and (args.both_strands or args.ambiguous_nucl or args.reverse_only):
            err('If you want to search sequence, set -q flag')

        if args.wrap and args.perl_regexp:
            err("PATTERNS found in --wrap captures must be literal (-P and -w incompatible)")

        if args.ambiguous_nucl:
            args.perl_regexp = True

        if args.force_color and args.no_color:
            err("WTF? --force-color AND --no-color?")

        # Some things just don't make sense in header searches ...
        if args.gff or args.ambiguous_nucl:
            args.match_sequence = True

        # Decide when to color
        if args.force_color or (sys.stdout.isatty() and not args.no_color):
            args.color = True
        else:
            args.color = False

        # Others don't make sense with color
        if args.gff or args.count_matches and args.color:
            args.color = False

        if args.gff:
            args.count = False
            args.count_matches = False

        return(args)

    def _create_matcher(self, args, pat, wrapper):
        # Check existence for matches to wrapper captures
        def swrpmatcher(text, strand='.'):
            for m in re.finditer(wrapper, text):
                if m.group(1) in pat:
                    return(True)
            return(False)

        # Check existence of matches
        def spatmatcher(text, strand='.'):
            for p in pat:
                if re.search(p, text):
                    return(True)
            return(False)

        # Find locations of matches to wrappers
        def gwrpmatcher(text, strand='.'):
            pos = []
            for m in re.finditer(wrapper, text):
                if m.group(1) in pat:
                    match = {'pos':(m.start(), m.end()), 'strand':strand}
                    pos.append(match)
            return(pos)

        # Find locations of matches
        def gpatmatcher(text, strand='.'):
            pos = []
            for p in pat:
                for m in re.finditer(p, text):
                    match = {'pos':(m.start(), m.end()), 'strand':strand}
                    pos.append(match)
            return(pos)

        if args.gff or args.count_matches or args.color:
            matcher = gwrpmatcher if wrapper else gpatmatcher
        else:
            matcher = swrpmatcher if wrapper else spatmatcher

        if args.reverse_only or args.both_strands:
            if matcher.__name__ in ('swrmatcher', 'spatmatcher'):
                if args.reverse_only:
                    def rmatcher(text):
                        match = matcher(FSeq.getrevcomp(text))
                        return(match)
                else:
                    def rmatcher(text):
                        match = matcher(text) + matcher(FSeq.getrevcomp(text))
                        return(match)
            else:
                def rev(matcher, text):
                    rmatch = []
                    for d in matcher(FSeq.getrevcomp(text), strand='-'):
                        d['pos'] = (len(text) - d['pos'][1], len(text) - d['pos'][0])
                        rmatch.append(d)
                    return(rmatch)
                if args.reverse_only:
                    def rmatcher(text):
                        return(rev(matcher, text))
                else:
                    def rmatcher(text):
                        fmatch = matcher(text)
                        rmatch = rev(matcher, text)
                        return(fmatch + rmatch)
            return(rmatcher)
        else:
            return(matcher)

    def _get_pattern(self, args):
        pat = set()
        if args.file:
            pat.update([l.rstrip('\n') for l in args.file])
        if args.patterns:
            pat.update(args.patterns)

        if args.ambiguous_nucl:
            apat = set()
            for p in pat:
                perlpat = ambiguous2perl(p)
                apat.update([perlpat])
            pat = apat

        if not pat:
            err('Please provide a pattern')

        # TODO searching for perfect matches would be faster without using
        # regex (just <str>.find(pat))
        if not args.perl_regexp and not args.wrap:
            pat = [re.escape(p) for p in pat]

        return(pat)

    def _makegen(self, args):
        if args.match_sequence:
            gettext = lambda x: x.seq
        else:
            gettext = lambda x: x.header

        if args.gff:
            def sgen(gen, matcher):
                source = "smof-{}".format(__version__)
                gfftype = args.gff_type
                row = [None,    #1 seqid
                       source,  #2 source
                       gfftype, #3 type
                       None,    #4 start
                       None,    #5 end
                       '.',     #6 score
                       None,    #7 strand
                       '.',     #8 phase
                       '.'      #9 attributes
                      ]
                for seq in gen.next():
                    row[0] = re.sub('^>(\S+)', '\1', seq.header)
                    for m in matcher(seq.seq):
                        row[3] = m['pos'][0] + 1
                        row[4] = m['pos'][1]
                        row[6] = m['strand']
                        yield('\t'.join([str(s) for s in row]))

        elif args.count or args.count_matches:
            def sgen(gen, matcher):
                count, matches = 0, 0
                for seq in gen.next():
                    text = gettext(seq)
                    m = matcher(text)
                    if (m and not args.invert_match) or (not m and args.invert_match):
                        count += 1
                        try:
                            matches += len(m)
                        except:
                            pass
                if args.count and args.count_matches:
                    yield "{}\t{}".format(count, matches)
                elif args.count:
                    yield count
                elif args.count_matches:
                    yield matches
        else:
            def sgen(gen, matcher):
                for seq in gen.next(handle_color=args.preserve_color):
                    text = gettext(seq)
                    m = matcher(text)
                    if (m and not args.invert_match) or (not m and args.invert_match):
                        if args.color:
                            for pos in [d['pos'] for d in m]:
                                if args.match_sequence:
                                    seq.color_seq(*pos, col=Colors.BOLD_RED)
                                else:
                                    seq.color_header(*pos, col=Colors.BOLD_RED)
                        yield(seq)
        return(sgen)

    def generator(self, args, gen):
        args = self._process_arguments(args)

        pat = self._get_pattern(args)

        flags = re.IGNORECASE if not args.case_sensitive else 0

        if args.wrap:
            wrapper = re.compile(args.wrap, flags=flags)
        else:
            pat = set((re.compile(p, flags=flags) for p in pat))
            wrapper = None

        matcher = self._create_matcher(args, pat, wrapper)

        sgen = self._makegen(args)

        self.force_color = args.force_color
        for item in sgen(gen, matcher):
            yield item

class Rename(Subcommand):
    def _parse(self):
        cmd_name = 'rename'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="Larry Wall's rename pythonized for headers"
        )
        parser.add_argument(
            'pattern',
            metavar='PATTERN',
            help="regex pattern"
        )
        parser.add_argument(
            'replacement',
            metavar='REPLACEMENT',
            nargs='?',
            default='',
            help="regex replacement (default='', i.e. delete PATTERN)"
        )
        parser.add_argument(
            'headers',
            metavar='HEADERS',
            nargs='?',
            default=None,
            help="regex to select headers (if not supplied, select all)"
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        matcher = re.compile(args.headers) if args.headers else None
        for seq in gen.next():
            if not matcher or re.search(matcher, seq.header):
                seq.header = re.sub(args.pattern, args.replacement, seq.header)
            yield seq

class Uniq(Subcommand):
    def _parse(self):
        cmd_name = 'uniq'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="prints entries if header/sequence pair is unique"
        )
        parser.add_argument(
            '-c', '--count',
            help='writes (count|header) in tab-delimited format',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-d', '--repeated',
            help='print only repeated entries',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-u', '--uniq',
            help='print only unique entries',
            action='store_true',
            default=False
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        from collections import OrderedDict
        seqs = OrderedDict()
        for seq in gen.next():
            try:
                seqs[seq] += 1
            except KeyError:
                seqs[seq] = 1

        if args.repeated:
            sgen = ((k,v) for k,v in seqs.items() if v > 1)
        elif args.uniq:
            sgen = ((k,v) for k,v in seqs.items() if v == 1)
        else:
            sgen = seqs.items()

        if args.count:
            for k,v in sgen:
                yield("{}\t{}".format(v, k.header))
        else:
            for k,v in sgen:
                yield(k)

class Wc(Subcommand):
    def _parse(self):
        cmd_name = 'wc'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="roughly emulates the UNIX wc command"
        )
        parser.add_argument(
            '-m', '--chars',
            help='writes the summed length of all sequences',
            action='store_true',
            default=False
        )
        parser.add_argument(
            '-l', '--lines',
            help='writes the total number of sequences',
            action='store_true',
            default=False
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        nchars, nseqs = 0, 0
        for seq in gen.next():
            nchars += len(seq.seq)
            nseqs += 1
        yield nseqs
        yield nchars

    def write(self, args, gen, out=sys.stdout):
        nseqs, nchars = list(self.generator(args, gen))
        if args.chars and not args.lines:
            out.write('%s\n' % nchars)
        elif args.lines and not args.chars:
            out.write('%s\n' % nseqs)
        else:
            out.write("{}\t{}\n".format(nseqs, nchars))

class Tail(Subcommand):
    def _parse(self):
        cmd_name = 'tail'
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="Write last N sequences (default=1)"
        )
        parser.add_argument(
            '-n', '--nseqs',
            help='print N sequences',
            metavar='N',
            type=counting_number,
            default=1
        )
        parser.add_argument(
            '-f', '--first',
            help='print first K letters of each sequence',
            metavar='K',
            type=counting_number
        )
        parser.add_argument(
            '-l', '--last',
            help='print last K letters of each sequence',
            metavar='K',
            type=counting_number
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        from collections import deque
        try:
            lastseqs = deque(maxlen=args.nseqs)
        except ValueError:
            err('--nseqs argument must be positive')
        for seq in gen.next():
            lastseqs.append(seq)

        for s in lastseqs:
            yield headtailtrunk(s, args.first, args.last)


# =======
# EXECUTE
# =======

if __name__ == '__main__':
    gen = FSeqGenerator()
    args = parse()
    args.func(args, gen, out=sys.stdout)
