import argparse
import math
import re
import sys
import string
import os
import signal
import textwrap
from itertools import chain
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
from hashlib import md5

from smof.version import __version__

# ================
# Argument Parsing
# ================


class Parser:
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            prog="smof",
            usage="<fastafile> | smof <subcommand> <options>",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="Tools for studying and manipulating fasta files",
            epilog=textwrap.dedent(
                """\
                Project site: https://github.com/incertae-sedis/smof
                Report bugs/requests via https://github.com/incertae-sedis/smof/issues
                Author: Zebulun Arendsee (zbwrnz@gmail.com)
            """
            ),
        )
        self.parser.add_argument(
            "-v",
            "--version",
            action="version",
            version="%(prog)s {}".format(__version__),
        )
        self.subparsers = self.parser.add_subparsers(
            metavar="[ for help on each: smof <subcommand> -h ]", title="subcommands"
        )
        self.usage = "<fastafile> | smof {} <options>"


def parse(argv=None):

    parser = Parser()

    # A list of valid subcommands, each of which is a class (defined below)
    subcommands = [
        Cut,
        Clean,
        Consensus,
        Filter,
        Grep,
        Md5sum,
        Head,
        Permute,
        Reverse,
        Sample,
        Sniff,
        Sort,
        Split,
        Stat,
        Subseq,
        Tail,
        Translate,
        Uniq,
        Wc,
    ]
    # add a subparser to parser for each subcommand
    for cmd in subcommands:
        cmd(parser)

    # the argv variable will be defined only for testing purposes
    # when it is not defined, take arguments from the terminal
    argv = argv if argv else sys.argv[1:]

    # if no argument is provided, print help
    if not argv:
        parser.parser.print_help()
        sys.exit(0)

    # handle attempted access to deprecated subcommands
    if argv[0] in ["rename", "fasta2csv"]:
        err("{} is deprecated".format(argv[0]))
    if argv[0] in ["idsearch", "retrieve", "search", "rmfields"]:
        err("{} is deprecated, use 'smof grep'".format(argv[0]))
    if argv[0] == "winnow":
        err("`winnow` is deprecated, use `smof filter`")
    if argv[0] == "chksum":
        err("`winnow` is deprecated, use `smof md5sum`")
    if argv[0] == "perm":
        err("`perm` is deprecated, use `smof permute`")

    # parse arguments, the default function that will ultimately be run is
    # selected according to tthe user-chosen subcommand
    args = parser.parser.parse_args(argv)

    return args


# =================
# CLASS DEFINITIONS
# =================


class Alphabet:
    # Including U, for selenocysteine, and * for STOP
    PROT = set("ACDEFGHIKLMNPQRSTUVWYX*")
    PROT_UNK = set("X")
    PROT_AMB = set("BZJ")
    PROT_EXC = set("EFILQPXJZ*")  # unique to protein sequences
    DNA = set("ACGTN")
    DNA_UNK = set("N")
    DNA_AMB = set("RYSWKMDBHV")
    RNA = set("ACGUN")
    RNA_UNK = set("N")
    RNA_AMB = set("RYSWKMDBHV")
    GAP = set(".-_")
    STOP = {"TAG", "TAA", "TGA", "UAG", "UAA", "UGA"}
    START = {"ATG", "AUG"}


class Colors:
    OFF = chr(27) + r"[0;0m"
    RED = chr(27) + r"[0;31m"
    GREEN = chr(27) + r"[0;32m"
    YELLOW = chr(27) + r"[0;33m"
    MAGENTA = chr(27) + r"[0;35m"
    CYAN = chr(27) + r"[0;36m"
    WHITE = chr(27) + r"[0;37m"
    BLUE = chr(27) + r"[0;34m"
    BOLD_RED = chr(27) + r"[1;31m"
    BOLD_GREEN = chr(27) + r"[1;32m"
    BOLD_YELLOW = chr(27) + r"[1;33m"
    BOLD_MAGENTA = chr(27) + r"[1;35m"
    BOLD_CYAN = chr(27) + r"[1;36m"
    BOLD_WHITE = chr(27) + r"[1;37m"
    BOLD_BLUE = chr(27) + r"[1;34m"

    patstr = chr(27) + r"\[[0-9;]+m"
    pat = re.compile(patstr)

    COLORS = {
        "red": RED,
        "green": GREEN,
        "yellow": YELLOW,
        "magenta": MAGENTA,
        "cyan": CYAN,
        "white": WHITE,
        "blue": BLUE,
        "bold_red": BOLD_RED,
        "bold_green": BOLD_GREEN,
        "bold_yellow": BOLD_YELLOW,
        "bold_magenta": BOLD_MAGENTA,
        "bold_cyan": BOLD_CYAN,
        "bold_white": BOLD_WHITE,
        "bold_blue": BOLD_BLUE,
    }


class ColorAA:
    def __init__(self):
        self.group = [
            ["LVGAIP", "aliphatic", Colors.BOLD_BLUE],
            ["FYW", "aromatic", Colors.BOLD_RED],
            ["SEKDRTNQH", "polar", Colors.BOLD_GREEN],
            ["MCU", "thiol", Colors.BOLD_YELLOW],
        ]
        # add lower cases
        self.group = [[l + l.lower(), g, c] for l, g, c in self.group]

    def color(self, a):
        for chars, group, color in self.group:
            if a in chars:
                return Colors.OFF + a + color
        return a


class ColorString:
    def __init__(self, seq=None, bgcolor=Colors.OFF, default=Colors.BOLD_RED):
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
            for s in re.split("((?:{})+)".format(Colors.patstr), thing):
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
            err("ColorString can only append strings, FSeq, or ColorString objects")
        self.cind += [
            [b + len(self.cind), e + len(self.cind), c] for b, e, c in newcind
        ]

    def subseq(self, a, b):
        self.cind = [x for x in self.cind if not (x[0] > b - 1 or x[1] - 1 < a)]
        self.cind = [[max(0, s - a), e - a, c] for s, e, c in self.cind]

    def reverse(self, length):
        self.cind = [[length - e, length - s, c] for s, e, c in self.cind]

    def colorpos(self, a, b, col=None):
        col = self.default if not col else col
        for i in reversed(range(len(self.cind))):
            start, end, color = self.cind[i]
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
        starts = {x[0]: x[2] for x in self.cind}
        ends = {x[1] for x in self.cind}
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
            if colwidth and i % colwidth == 0 and i != 0:
                out.write("\n")
            out.write(seq[i])
        out.write(self.bgcolor if colored else "")
        out.write("\n")

    def copy(self):
        new_obj = ColorString(bgcolor=self.bgcolor, default=self.default)
        new_obj.cind = self.cind
        return new_obj


class FileDescription:
    def __init__(self):
        self.seqs = set()
        self.headers = set()
        self.ntype = {"prot": 0, "dna": 0, "rna": 0, "illegal": 0, "ambiguous": 0}
        self.ncase = {"uppercase": 0, "lowercase": 0, "mixedcase": 0}
        self.pfeat = {
            "selenocysteine": 0,
            "initial-Met": 0,
            "internal-stop": 0,
            "terminal-stop": 0,
        }
        from itertools import product

        self.nfeat = {"".join(c): 0 for c in product("01", repeat=4)}
        self.ufeat = {"gapped": 0, "unknown": 0, "ambiguous": 0}

    def add_seq(self, seq):
        """
        Calculates properties for one sequence
        @type seq: FSeq object
        """
        # Add md5 hashes of sequences and headers to respective sets
        self.seqs.update([md5(bytes(seq.seq, "ascii")).digest()])
        self.headers.update([md5(bytes(seq.header, "ascii")).digest()])

        counts = Counter(seq.seq)

        # Ungapped the sequence is required for downstream analysis
        is_gapped = self._handle_gaps(seq, counts)
        if is_gapped:
            seq.ungap()
            counts = Counter(seq.seq)

        scase = self._handle_case(counts)
        # Sum upper and lowercase
        if scase not in "upper":
            counts = counter_caser(counts)

        s = seq.seq.upper()

        # ('prot'|'dna'|'rna'|'amb'|'bad')
        stype = self._handle_type(counts)

        if stype == "prot":
            tstop = "*" == s[-1]
            self.pfeat["terminal-stop"] += tstop
            self.pfeat["internal-stop"] += (counts["*"] - tstop) > 0
            self.pfeat["selenocysteine"] += "U" in counts
            self.pfeat["initial-Met"] += "M" == s[0]
            self.ufeat["unknown"] += "X" in counts
            self.ufeat["ambiguous"] += bool(Alphabet.PROT_AMB & set(counts))
        elif stype in ("dna", "rna"):
            self.ufeat["unknown"] += bool("N" in counts)
            self.ufeat["ambiguous"] += bool(Alphabet.DNA_AMB & set(counts))

            start = self._has_start(s)
            stop = self._has_stop(s)
            triple = self._is_triple(s)
            sense = self._is_sense(s)

            profile = "".join([str(int(x)) for x in (start, stop, triple, sense)])
            self.nfeat[profile] += 1

    def get_nseqs(self):
        return sum(self.ntype.values())

    def count_degenerate_headers(self):
        return self.get_nseqs() - len(self.headers)

    def count_degenerate_seqs(self):
        return self.get_nseqs() - len(self.seqs)

    @classmethod
    def _has_start(cls, s):
        """
        Tests if the first codon is the START codon, assumes uppercase
        """
        return s[0:3] == "ATG"

    @classmethod
    def _has_stop(cls, s):
        """
        Tests if the last three letters are a STOP codon, assumes uppercase
        """
        return s[-3:] in Alphabet.STOP

    @classmethod
    def _is_triple(cls, s):
        """
        Tests is the sequence is multiple of three
        """
        return len(s) % 3 == 0

    @classmethod
    def _is_sense(cls, s):
        """
        Tests if there are no internal STOP codons in the first frame
        """
        codons = set(s[i : i + 3] for i in range(0, len(s) - 3, 3))
        return bool(not codons & Alphabet.STOP)

    def _handle_gaps(self, seq, counts):
        # Handle gaps
        if Alphabet.GAP & set(counts):
            self.ufeat["gapped"] += 1
            return True
        return False

    def _handle_case(self, counts):
        has_lower = set(string.ascii_lowercase) & set(counts)
        has_upper = set(string.ascii_uppercase) & set(counts)
        if has_lower and has_upper:
            case = "mixedcase"
        elif has_lower:
            case = "lowercase"
        else:
            case = "uppercase"
        self.ncase[case] += 1
        return case

    def _handle_type(self, counts):
        stype = guess_type(counts)
        self.ntype[stype] += 1
        return stype


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
    # Extended alphabet:
    # W = [AT]  <--> S = [GC]
    # M = [AC]  <--> K = [GT]
    # R = [AG]  <--> Y = [CT]
    # B = [GTC] <--> V = [ACG]
    # D = [AGT] <--> H = [ACT]
    # N <--> N
    revtrans = str.maketrans(
        "acgtuwsmkrybvdhnACGTUWSMKRYBVDHN", "tgcaaswkmyrvbhdnTGCAASWKMYRVBHDN"
    )
    # Translator of ungapping
    ungapper = str.maketrans("", "", "".join(Alphabet.GAP))
    # Codon table for translation
    codon_table = {
        "TTT": "F",
        "TCT": "S",
        "TAT": "Y",
        "TGT": "C",
        "TTC": "F",
        "TCC": "S",
        "TAC": "Y",
        "TGC": "C",
        "TTA": "L",
        "TCA": "S",
        "TAA": "*",
        "TGA": "*",
        "TTG": "L",
        "TCG": "S",
        "TAG": "*",
        "TGG": "W",
        "CTT": "L",
        "CCT": "P",
        "CAT": "H",
        "CGT": "R",
        "CTC": "L",
        "CCC": "P",
        "CAC": "H",
        "CGC": "R",
        "CTA": "L",
        "CCA": "P",
        "CAA": "Q",
        "CGA": "R",
        "CTG": "L",
        "CCG": "P",
        "CAG": "Q",
        "CGG": "R",
        "ATT": "I",
        "ACT": "T",
        "AAT": "N",
        "AGT": "S",
        "ATC": "I",
        "ACC": "T",
        "AAC": "N",
        "AGC": "S",
        "ATA": "I",
        "ACA": "T",
        "AAA": "K",
        "AGA": "R",
        "ATG": "M",
        "ACG": "T",
        "AAG": "K",
        "AGG": "R",
        "GTT": "V",
        "GCT": "A",
        "GAT": "D",
        "GGT": "G",
        "GTC": "V",
        "GCC": "A",
        "GAC": "D",
        "GGC": "G",
        "GTA": "V",
        "GCA": "A",
        "GAA": "E",
        "GGA": "G",
        "GTG": "V",
        "GCG": "A",
        "GAG": "E",
        "GGG": "G",
    }

    def __init__(
        self, header, seq, filename=None, handle_color=False, purge_color=False
    ):
        self.seq = seq
        self.header = header
        self.colseq = None
        self.colheader = None
        self.handle_color = handle_color
        self.filename = filename
        self.moltype = None
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
        return hash((self.header, self.seq))

    def __eq__(self, other):
        return (self.header, self.seq) == (other.header, other.seq)

    @staticmethod
    def _clear_color(text):
        return re.sub(Colors.pat, "", text)

    def color_seq(self, *args, **kwargs):
        if not self.colseq:
            self._process_color()
        self.colseq.colorpos(*args, **kwargs)

    def color_header(self, *args, **kwargs):
        if not self.colheader:
            self._process_color()
        self.colheader.colorpos(*args, **kwargs)

    def ungap(self):
        self.seq = self.seq.translate(FSeq.ungapper)
        self.header = ParseHeader.add_suffix(self.header, "ungapped")

    def print(self, col_width=80, color=True, out=sys.stdout):
        out.write(">")
        if self.colheader and color:
            self.colheader.print(self.header, colwidth=None, out=out)
        else:
            out.write("%s\n" % self.header)
        for i in range(0, len(self.seq), col_width):
            if self.colseq and color:
                self.colseq.print(self.seq, col_width, out=out)
                break
            else:
                out.write("%s\n" % self.seq[i : i + col_width])

    def get_pretty_string(self, col_width=80):
        out = [">" + self.header]
        for i in range(0, len(self.seq), col_width):
            out.append(self.seq[i : i + col_width])
        outstr = "\n".join(out)
        return outstr

    def seq_upper(self):
        self.seq = self.seq.upper()

    def set_moltype(self, moltype=None):
        if moltype:
            self.moltype = moltype
        else:
            self.moltype = guess_type(self)

    def get_moltype(self):
        if not self.moltype:
            self.set_moltype()
        return self.moltype

    def header_upper(self):
        self.header = self.header.upper()

    def subseq(self, a, b):
        header = ParseHeader.subseq(self.header, a + 1, b)
        newseq = FSeq(header, self.seq[a:b])
        if self.colseq:
            newseq.colseq = self.colseq.copy()
            newseq.colseq.subseq(a, b)
        return newseq

    def add_filename(self):
        if self.filename:
            self.header = ParseHeader.add_tag(
                h=self.header, tag="filename", value=self.filename
            )
            self.colheader = None

    def reverse(self):
        self.seq = self.seq[::-1]
        if self.handle_color:
            self.colseq.reverse(len(self.seq))
        self.header = ParseHeader.add_suffix(self.header, "reverse")

    @classmethod
    def getrevcomp(cls, seq):
        trans = lambda s: s[::-1].translate(FSeq.revtrans)
        if isinstance(seq, str):
            return trans(seq)
        elif isinstance(seq, FSeq):
            newheader = ParseHeader.add_suffix(seq.header, "revcom")
            newseq = FSeq(newheader, trans(seq.seq))
            if seq.colseq:
                newseq.colseq = seq.colseq
                newseq.colseq.reverse(len(seq.seq))
            if seq.colheader:
                # TODO implement this
                pass
            return newseq


class FSeqGenerator:
    def __init__(self, args):
        self.args = args

    def next(self, *args, **kwargs):
        # If no input is given,
        # and if smof is not reading user input from stdin,
        # assume piped input is from STDIN
        try:
            if not self.args.fh:
                fh = [sys.stdin]
            else:
                fh = self.args.fh
        # If args does not have a .fh argument, then try treating args itself
        # as the input
        except AttributeError:
            fh = [self.args]

        for fastafile in fh:
            seq_list = []
            header = None
            # If there are multiple input files, store the filename
            # If there is only one, e.g. STDIN, don't store a name
            filename = None if len(fh) == 1 else fastafile
            try:
                f = open(fastafile, "r")
            except TypeError:
                f = fastafile
            except FileNotFoundError:
                err("File '%s' not found" % fastafile)

            for line in f:
                line = line.strip()
                if not line or line[0] == "#":
                    continue
                if line[0] == ">":
                    if seq_list:
                        yield FSeq(
                            header,
                            "".join(seq_list),
                            filename=filename,
                            *args,
                            **kwargs
                        )
                    elif header:
                        # NOTE: yields an empty sequence! This is usually
                        # a BAD THING, but it can happen in the wild
                        yield FSeq(header, "", filename=filename, *args, **kwargs)
                    seq_list = []
                    header = line[1:]
                # '' is valid for a header
                elif header is not None:
                    seq_list.append(line)
                else:
                    err("First fasta line must begin with '>'")
            # process the last sequence
            if header is not None:
                if seq_list:
                    yield FSeq(
                        header, "".join(seq_list), filename=filename, *args, **kwargs
                    )
                else:
                    # NOTE: yields empty sequence!
                    yield FSeq(header, "", filename=filename, *args, **kwargs)

            try:
                f.close()
            except AttributeError:
                pass


class Maps:
    DNA_AMB = {
        "R": "AG",
        "Y": "CT",
        "S": "GC",
        "W": "AT",
        "K": "GT",
        "M": "AC",
        "B": "CGT",
        "D": "AGT",
        "H": "ACT",
        "V": "ACG",
        "N": "ACGT",
    }


class ParseHeader:
    @staticmethod
    def firstword(h, delimiter=" \t"):
        return re.sub("^([^%s]+).*" % delimiter, "\\1", h)

    @staticmethod
    def description(h):
        return re.sub(r"^\S+\s*", "", h)

    @staticmethod
    def add_suffix(h, suffix):
        return re.sub(r"^(\S+)(.*)", "\\1|%s\\2" % suffix, h)

    @staticmethod
    def add_tag(h, tag, value):
        return re.sub(r"^(\S+)(.*)", "\\1 %s=%s\\2" % (tag, value), h)

    @staticmethod
    def subseq(h, a, b):
        header = "%s|subseq(%d..%d) %s" % (
            ParseHeader.firstword(h),
            a,
            b,
            ParseHeader.description(h),
        )
        return header.strip()

    @staticmethod
    def permute(h, start, end, wordsize):
        header = "%s|permutation:start=%d;end=%d;word_size=%d %s" % (
            ParseHeader.firstword(h),
            start,
            end,
            wordsize,
            ParseHeader.description(h),
        )
        return header.strip()

    @staticmethod
    def ncbi_format(h, fields):
        raise NotImplementedError

    @staticmethod
    def regex_group(h, regex):
        raise NotImplementedError


class SeqStat:
    def __init__(self, seq, count=True):
        self.counts = Counter(seq.seq) if count else None
        self.header = seq.header
        self.length = len(seq.seq)

    def aslist(
        self,
        charset=None,
        length=False,
        masked=False,
        header_fun=None,
        ignorecase=False,
    ):

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
            charset = set("".join(charset).upper())
            self.counts = counter_caser(self.counts)

        line += [self.counts[c] for c in sorted(charset)]

        return line

    @classmethod
    def getheader(cls, charset, length=False, masked=False, ignorecase=False):
        header = ["seqid"]

        if ignorecase:
            charset = set("".join(charset).upper())

        if length:
            header += ["length"]

        if masked:
            header += ["masked"]

        header += list(sorted(charset))

        return header


class StatFun:
    @classmethod
    def N50(cls, xs, issorted=False):
        xs = sorted(xs) if not issorted else xs
        N = sum(xs)
        total = 0
        for i in range(len(xs) - 1, -1, -1):
            total += xs[i]
            if total > N / 2:
                return xs[i]

    @classmethod
    def mean(cls, xs):
        if not xs:
            mu = float("nan")
        else:
            mu = sum(xs) / len(xs)
        return mu

    @classmethod
    def median(cls, xs, issorted=False):
        return cls.quantile(xs, 0.5, issorted=issorted)

    @classmethod
    def sd(cls, xs):
        if len(xs) < 2:
            stdev = float("nan")
        else:
            mean = sum(xs) / len(xs)
            stdev = (sum((y - mean) ** 2 for y in xs) / (len(xs) - 1)) ** 0.5
        return stdev

    @classmethod
    def quantile(cls, xs, q, issorted=False):
        """
        Calculates quantile as the weighted average between indices
        """
        # Die if out of bounds
        if not 0 <= q <= 1:
            err("quantile must be between 0 and 1")

        # Ensure the vector is sorted
        xs = sorted(xs) if not issorted else xs

        # Return max or min for q = 1 or 0
        if q == 1:
            return xs[-1]
        elif q == 0:
            return xs[0]

        v = (len(xs) - 1) * q
        r = v % 1
        i = math.floor(v)
        quantile = xs[i] * (1 - r) + xs[i + 1] * r
        return quantile

    @classmethod
    def summary(cls, xs):
        xs = sorted(xs)
        out = {
            "min": xs[0],
            "max": xs[-1],
            "1st_qu": cls.quantile(xs, 0.25, issorted=True),
            "median": cls.quantile(xs, 0.50, issorted=True),
            "3rd_qu": cls.quantile(xs, 0.75, issorted=True),
            "mean": cls.mean(xs),
            "sd": cls.sd(xs),
            "N50": cls.N50(xs, issorted=True),
        }
        return out


# =================
# UTILITY FUNCTIONS
# =================


def counter_caser(counter, lower=False):
    """
    Sums cases in Collections.Counter object
    """
    if lower:
        out = counter + Counter(
            {k.lower(): v for k, v in counter.items() if k.isupper()}
        )
        out = out - Counter({k: v for k, v in counter.items() if k.isupper()})
    else:
        out = counter + Counter(
            {k.upper(): v for k, v in counter.items() if k.islower()}
        )
        out = out - Counter({k: v for k, v in counter.items() if k.islower()})
    return out


def sum_lower(counter):
    lc = [v for k, v in counter.items() if k in string.ascii_lowercase]
    return sum(lc)


def guess_type(counts):
    """
    Predict sequence type from character counts (dna|rna|prot|ambiguous|illegal)
    """
    if isinstance(counts, str):
        counts = Counter(counts)
    elif isinstance(counts, FSeq):
        counts = Counter(counts.seq)

    # Convert all to upper case
    counts = counter_caser(counts)
    # Remove gaps from Counter
    counts = Counter({k: n for k, n in counts.items() if k not in Alphabet.GAP})

    # If all chars are in ACGT
    if set(counts) <= Alphabet.DNA:
        stype = "ambiguous" if sum(counts.values()) < 3 else "dna"
    # If all chars in ACGU
    elif set(counts) <= Alphabet.RNA:
        stype = "rna"
    # If has any chars unique to proteins (EFILQPXJZ*)
    elif set(counts) & Alphabet.PROT_EXC:
        if set(counts) <= Alphabet.PROT | Alphabet.PROT_AMB:
            stype = "prot"
        else:
            stype = "illegal"
    # If all the residues could be aa, DNA, or RNA
    elif set(counts) <= (Alphabet.PROT | Alphabet.PROT_AMB):
        # If more than 80% look like nucleic acids, set 'amb_nucl'
        if (
            sum([counts[x] for x in "ACGTUN" if x in counts]) / sum(counts.values())
        ) > 0.8:
            if "U" in counts:
                stype = "illegal" if "T" in counts else "rna"
            else:
                stype = "dna"
        # Otherwise set as ambibuous
        else:
            stype = "ambiguous"
    # If none of these match, something is horribly wrong with your
    # sequence
    else:
        stype = "illegal"
    return stype


def headtailtrunk(seq, first=None, last=None):
    """
    This function is used by the Head and Tail classes to portray partial
    sections of sequences.
    """
    outseq = FSeq(seq.header, seq.seq)
    if first and last:
        if first + last < len(seq.seq):
            outseq.header = ParseHeader.firstword(
                seq.header
            ) + "|TRUNCATED:first-{}_last-{}".format(first, last)
            outseq.seq = "{}{}{}".format(seq.seq[0:first], "...", seq.seq[-last:])
    elif first:
        outseq.header = ParseHeader.firstword(
            seq.header
        ) + "|TRUNCATED:first-{}".format(first)
        outseq.seq = seq.seq[0:first]
    elif last:
        outseq.header = ParseHeader.firstword(seq.header) + "|TRUNCATED:last-{}".format(
            last
        )
        outseq.seq = seq.seq[-last:]
    elif first == 0 and last == 0:
        err("Illegal empty sequence, dying ...")
    return outseq


def ascii_histchar(dif, chars=" .~*O"):
    if dif <= 0:
        return chars[0]
    elif dif < 0.25:
        return chars[1]
    elif dif < 0.5:
        return chars[2]
    elif dif < 0.75:
        return chars[3]
    else:
        return chars[4]


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
    sys.exit(msg)


def ambiguous2perl(pattern):
    perlpat = []
    in_bracket = False
    escaped = False
    for c in pattern:
        amb = c in Maps.DNA_AMB
        if c == "\\":
            escaped = True
            continue
        elif escaped:
            c = c if amb else "\\" + c
            escaped = False
        elif amb:
            v = Maps.DNA_AMB[c]
            c = v if in_bracket else "[%s]" % v
        elif c == "[":
            in_bracket = True
        elif c == "]":
            in_bracket = False
        perlpat.append(c)
    return "".join(perlpat)


def find_max_orf(dna, from_start=False):
    dna = dna.translate(FSeq.ungapper).upper()
    max_start = None
    max_length = 0
    for offset in [2, 3, 4]:
        start = None
        length = 0
        for i in range(offset, len(dna), 3):
            codon = dna[i - 2 : i + 1]
            # if we encounter a STOP codon
            if codon in Alphabet.STOP:
                # if this ORF is the longest, record it
                if start != None and ((length, max_start) > (max_length, start)):
                    max_start, max_length = start, length
                # reset the ORF
                start, length = None, 0
                continue
            # if all ORFs are required to start with START codons, and if this
            # is not a START codon, then do not start a new ORF
            if from_start and start == None and not (codon in Alphabet.START):
                continue
            # if we are not currently in an ORF, initialize one
            if not start != None:
                start, length = i - 2, 0
            # increment the length of the current ORF
            length += 3
        # if the final ORF is the longest, record it
        # break ties by preferring the lower start position
        if start != None and ((length, max_start) > (max_length, start)):
          max_start, max_length = start, length
    return (max_start, max_length)


def translate_dna(dna):
    # remove gaps
    dna = dna.translate(FSeq.ungapper).upper()
    aa = []
    for i in range(2, len(dna), 3):
        codon = dna[i - 2 : i + 1]
        if codon in FSeq.codon_table:
            aa.append(FSeq.codon_table[codon])
        else:
            aa.append("X")
    return ''.join(aa)


def get_orf(dna, all_frames=False, from_start=False, translate=True):
    if not all_frames:
        aa = translate_dna(dna)
        if from_start and not aa[0] == "M":
          aa = ""
        return aa
    else:
        cds_start, cds_length = find_max_orf(dna, from_start=from_start)
        if cds_start is None:
          return ""
        cds = dna[cds_start : cds_start + cds_length]
        if translate:
          return translate_dna(cds)
        else:
          return cds


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
        raise NotImplementedError

    def generator(self, args, gen):
        raise NotImplementedError

    def write(self, args, gen, out=sys.stdout):
        for output in self.generator(args, gen):
            if isinstance(output, FSeq):
                color = sys.stdout.isatty() or self.force_color
                output.print(color=color, out=out)
            else:
                out.write("%s\n" % output)


class Clean(Subcommand):
    def _parse(self):
        cmd_name = "clean"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="cleans fasta files",
            description="""Remove all space within the sequences and write them
            in even columns (default width of 80 characters). Case and all
            characters (except whitespace) are preserved by default.""",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument("-t", "--type", metavar="n|p", help="sequence type")
        parser.add_argument(
            "-u",
            "--toupper",
            help="convert to uppercase",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-l",
            "--tolower",
            help="convert to lowercase",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-x",
            "--toseq",
            help="removes all non-letter characters (gaps, stops, etc.)",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-s",
            "--reduce-header",
            help="Remove all text from header that follows the first word (delimited by [ |])",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-r",
            "--mask-irregular",
            help="converts irregular letters to unknown",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-m",
            "--mask-lowercase",
            help="convert lower-case to unknown",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-w",
            "--col_width",
            help="width of the sequence output (0 indicates no wrapping)",
            metavar="W",
            type=positive_int,
            default=80,
        )
        parser.add_argument(
            "-d",
            "--standardize",
            help="Convert 'X' in DNA to 'N' and '[._]' to '-' (for gaps)",
            action="store_true",
            default=False,
        )
        parser.set_defaults(func=self.func)

    @staticmethod
    def _process_args(args):
        if (args.mask_lowercase or args.mask_irregular) and not args.type:
            err("Please provide sequence type (--type)")

        if args.tolower and args.toupper:
            err("Err, you want me to convert to lower AND upper?")

    def generator(self, args, gen):

        # Catches illegal combinations of arguments
        self._process_args(args)

        # Make translator, this is only possible if type is given
        trans = None
        if args.type:
            if args.type.lower()[0] in ["a", "p"]:
                # irregulars proteins include selenocysteine (U) and protein
                # ambiguous characters
                irr = "".join(Alphabet.PROT_AMB) + "U"
                unk = "X"
                standard_trans = str.maketrans("._", "--")
            elif args.type.lower()[0] in ["n", "d"]:
                # irregular nucleotides are ambiguous characters
                irr = "".join(Alphabet.DNA_AMB)
                unk = "N"
                standard_trans = str.maketrans("Xx._", "Nn--")
            else:
                err("Type not recognized")

            a = ""
            # Get irregular characters
            if args.mask_irregular:
                a += irr.upper() + irr.lower()

            # convert lowercase to unknown
            if args.mask_lowercase:
                a += string.ascii_lowercase

            a = "".join(set(a))

            b = unk * len(a)
            if irr:
                trans = str.maketrans(a, b)

        for seq in gen.next(purge_color=True):
            if args.reduce_header:
                seq.header = ParseHeader.firstword(seq.header, delimiter=" \t|")

            if args.standardize:
                try:
                    seq.seq = seq.seq.translate(standard_trans)
                except UnboundLocalError:
                    err("Please provide a type argument (-t)")

            # WARNING: order is important here, don't swap thoughtlesly
            # Remove all nonletters or wanted, otherwise just remove space
            if args.toseq:
                seq.seq = re.sub(r"[^A-Za-z]", "", seq.seq)
            else:
                seq.seq = re.sub(r"[^\S]", "", seq.seq)

            # Irregular or lowercase to unknown
            if trans:
                seq.seq = seq.seq.translate(trans)

            # Change everything to desired case
            if args.toupper:
                seq.seq = seq.seq.upper()
            elif args.tolower:
                seq.seq = seq.seq.lower()

            yield seq

    def write(self, args, gen, out=sys.stdout):
        if args.col_width == 0:
            args.col_width = int(1e12)  # Approximation of infinity, i.e. no wrap
        for seq in self.generator(args, gen):
            seq.print(col_width=args.col_width, color=False, out=out)


class Filter(Subcommand):
    def _parse(self):
        cmd_name = "filter"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="extracts sequences meeting the given conditions",
            description="""Prints every entry by default. You may add one or
            more criteria to filter the results (e.g. `smof filter -s 200 -l
            100 -c 'GC > .5'` will print only sequences between 100 and 200
            resides in length and greater than 50% GC content).""",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "-s",
            "--shorter-than",
            help="keep only if length is less than or equal to LEN",
            type=positive_int,  # 0 is valid, it would keep only 0-length sequences
            metavar="LEN",
        )
        parser.add_argument(
            "-l",
            "--longer-than",
            help="keep only if length is greater than or equal to LEN",
            type=positive_int,  # 0 is valid (but pointless) -- it will keep everything
            metavar="LEN",
        )
        parser.add_argument(
            "-c",
            "--composition",
            metavar="EXPR",
            help="keep only if composition meets the condition (e.g. 'GC > .5')",
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        tests = []
        if args.shorter_than is not None:
            # args.shorter_than CAN be 0, so explicitly match to None
            tests.append(lambda s, v=args.shorter_than: len(s) <= v)
        if args.longer_than is not None:
            # args.longer_than CAN be 0, so explicitly match to None
            tests.append(lambda s, v=args.longer_than: len(s) >= v)
        if args.composition:
            try:
                ch, sign, per = args.composition.split()
            except:
                err(
                    'The argument for --composition must be three space separated values, e.g. "GC > .5"'
                )
            legal_signs = ("<", "<=", ">=", ">", "==", "=", "!=")
            if not sign in legal_signs:
                err(
                    "Middle term must be a comparison symbol ('<', '<=', '>=', '>', '==', '=', '!=')"
                )
            if sign == "=":
                sign = "=="
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
                return eval("p {} {}".format(sign, per))

            tests.append(evaluate)

        for seq in gen.next():
            accept = all([x(seq.seq) for x in tests])
            if accept:
                yield seq


class Md5sum(Subcommand):
    def _parse(self):
        cmd_name = "md5sum"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="calculate an md5 checksum for the input sequences",
            description="""By default, `smof md5sum` concantenates all headers
            and sequences and calculates the md5sum for the resulting string.
            This is identical to `tr -d '\\n>' < a.fa | md5sum`.""",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "-i",
            "--ignore-case",
            help="convert all to uppercase, before hashing",
            action="store_true",
            default=False,
        )
        method = parser.add_mutually_exclusive_group(required=False)
        method.add_argument(
            "-q",
            "--each-sequence",
            help="calculate md5sum for each sequence (TAB delimited)",
            action="store_true",
            default=False,
        )
        method.add_argument(
            "-s",
            "--all-sequences",
            help="calculate one md5sum for all concatenated sequences",
            action="store_true",
            default=False,
        )
        method.add_argument(
            "-d",
            "--all-headers",
            help="calculate one md5sum for all concatenated headers",
            action="store_true",
            default=False,
        )
        method.add_argument(
            "-r",
            "--replace-header",
            help="replace the header of each entry with the checksum of the sequence",
            action="store_true",
            default=False,
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        md5hash = md5()
        # Hash the sequences only (in input order)
        if args.all_sequences:
            fun = lambda s, h: md5hash.update(s)
        # Hash the headers only (in input order)
        elif args.all_headers:
            fun = lambda s, h: md5hash.update(h)
        # DEFAULT: Hash headers and sequences (concatenated)
        # Equivalent to:
        # $ tr -d '\n>' < myfile.fa | md5sum
        else:
            fun = lambda s, h: md5hash.update(h + s)

        for seq in gen.next():
            if args.ignore_case:
                seq.header_upper()
                seq.seq_upper()
            s = seq.seq.encode("ascii")
            h = seq.header.encode("ascii")
            # Write <header>\t<sequence hash> for each sequence
            if args.replace_header:
                yield FSeq(md5(s).hexdigest(), seq.seq)
            elif args.each_sequence:
                yield "{}\t{}".format(
                    ParseHeader.firstword(seq.header), md5(s).hexdigest()
                )
            else:
                fun(s, h)

        # Print output hash for cumulative options
        if not (args.each_sequence or args.replace_header):
            yield md5hash.hexdigest()


class Permute(Subcommand):
    def _parse(self):
        cmd_name = "permute"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="randomly order sequence",
            description="""Randomly order letters in each sequence. The
            --word-size option allows random ordering of words of the given
            size. The --start-offset and --end-offset options are useful if,
            for example, you want to rearrange the letters within a coding
            sequence but want to preserve the start and stop codons.""",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "-w",
            "--word-size",
            help="size of each word (default=1)",
            type=counting_number,
            metavar="INT",
            default=1,
        )
        parser.add_argument(
            "-s",
            "--start-offset",
            help="number of letters to ignore at beginning (default=0)",
            type=positive_int,
            metavar="INT",
            default=0,
        )
        parser.add_argument(
            "-e",
            "--end-offset",
            help="number of letters to ignore at end (default=0)",
            type=positive_int,
            metavar="INT",
            default=0,
        )
        parser.add_argument(
            "--seed",
            help="set random seed (for reproducibility/debugging)",
            type=counting_number,
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
            suffix = s[L - end : L]
            rseq = s[start : L - end]
            M = len(rseq)
            words = list(rseq[i : i + w] for i in range(0, M - w + 1, w))
            words.append(rseq[(M - M % w) : M])
            random.shuffle(words)
            out = "".join(prefix + "".join(words) + suffix)
            header = ParseHeader.permute(seq.header, start, end, w)
            yield FSeq(header, out)


class Reverse(Subcommand):
    def _parse(self):
        cmd_name = "reverse"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="reverse each sequence (or reverse complement)",
            description="""Reverse the letters in each sequence. The complement
            is NOT taken unless the -c flag is set. The extended nucleotide
            alphabet is supported.""",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "-c",
            "--complement",
            help="take the reverse complement of the sequence",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-V",
            "--no-validate",
            help="do not check whether the sequence is DNA before reverse complement",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-S",
            "--preserve-color",
            help="Preserve incoming color",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-Y",
            "--force-color",
            help="print in color even to non-tty (DANGEROUS)",
            action="store_true",
            default=False,
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        """ Reverse each sequence """
        self.force_color = args.force_color
        if args.complement:

            def f(s):
                return FSeq.getrevcomp(s)

        else:

            def f(s):
                s.reverse()
                return s

        if args.complement and not args.no_validate:

            def func(s):
                if s.get_moltype() == "dna":
                    return f(s)
                else:
                    msg = "Cannot take reverse complement of the sequence '%s' since it does not appear to DNA"
                    err(msg % ParseHeader.firstword(s.header))

        else:
            func = f

        for seq in gen.next(handle_color=args.preserve_color):
            yield func(seq)


class Sniff(Subcommand):
    def _parse(self):
        cmd_name = "sniff"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="extract info about the sequence",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="Identifies the sequence type and aids in diagnostics.",
            epilog=textwrap.dedent(
                """\
                The output can be divided into 6 sections

                1. Overview and warnings

                  smof sniff counts the number of unique sequences and the number
                  of total sequences. It warns if there are any sequences with
                  illegal characters or if there are any duplicate headers. Example:

                  > 23 uniq sequences (24 total)
                  > WARNING: headers are not unique (23/24)
                  > WARNING: illegal characters found

                  Illegal characters include any character that is neither
                  standard, ambiguous, a gap [_-.], or a stop [*].

                2. Sequence types

                  For each entry, it predicts whether it is protein, DNA, RNA, or
                  illegal. Example:

                  > Sequence types:
                  >   prot:                20         83.3333%
                  >   dna:                 2          8.3333%
                  >   illegal:             1          4.1667%
                  >   rna:                 1          4.1667%

                  The 2nd column is the count, 3rd percentage

                3. Sequence cases

                  Reports the case of the sequences, example:

                  > Sequences cases:
                  >   uppercase:           21         87.5000%
                  >   lowercase:           2          8.3333%
                  >   mixedcase:           1          4.1667%

                4. Nucleotide features

                  Reports a summary nucleotide features

                  The nucleotide features entry is comprised of four flags
                  which will all equal 1 for a proper nucleotide coding sequence
                  (0 otherwise). A sequence will be counted as 1111 if it:

                    1) starts with a start codon (ATG)
                    2) ends with a stop codon (TAG, TAA, or TGA)
                    3) has a length that is a multiple of three
                    4) has no internal stop codon. If a sequence lacks a
                       start codon, but otherwise looks like a coding sequence,
                       it will have the value 0111.

                  For example:

                  > Nucleotide Features
                  >   0000:                2          66.6667%
                  >   1100:                1          33.3333%


                5. Protein features

                  1) terminal-stop - does the sequence end with '*'?
                  2) initial-Met - does the sequence start with 'M'?
                  3) internal-stop - does '*' appear within the sequence?
                  4) selenocysteine - does the sequence include 'U'?

                  Example:

                  > Protein Features:
                  >   terminal-stop:       20         100.0000%
                  >   initial-Met:         19         95.0000%
                  >   internal-stop:       0          0.0000%
                  >   selenocysteine:      0          0.0000%

                6. Universal features

                  Example:

                  > Universal Features:
                  >   ambiguous:           1          4.1667%
                  >   unknown:             0          0.0000%
                  >   gapped:              0          0.0000%

                Ambiguous characters are RYSWKMDBHV for nucleotides and BJZ
                for proteins.

                Unknown characters are X for proteins and N for nucleotides

                Gaps are '-_.'

            """
            ),
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        seqsum = FileDescription()
        for seq in gen.next():
            seqsum.add_seq(seq)
        yield seqsum

    def write(self, args, gen, out=sys.stdout):
        """
        This function basically just formats and prints the information in a
        FileDescription object
        """
        # The generator yields only this one item: a FileDescription object
        seqsum = next(self.generator(args, gen))

        # Total number of sequences
        nseqs = seqsum.get_nseqs()

        # Print number of uniq and total sequences
        if seqsum.count_degenerate_seqs():
            uniq = nseqs - seqsum.count_degenerate_seqs()
            out.write("{} uniq sequences ({} total)\n".format(uniq, nseqs))
        else:
            out.write("Total sequences: {}\n".format(nseqs))

        # Warn if there are any duplicate headers
        if seqsum.count_degenerate_headers():
            uniq = nseqs - seqsum.count_degenerate_headers()
            out.write("WARNING: headers are not unique ({}/{})\n".format(uniq, nseqs))

        # Warn if there are any illegal characters
        if seqsum.ntype["illegal"]:
            out.write("WARNING: illegal characters found\n")

        def write_dict(d, name, N):
            # Print keys if value is greater than 0
            uniq = [[k, v] for k, v in d.items() if v > 0]
            # E.g. If all of the sequences are proteins, print 'All prot'
            if len(uniq) == 1:
                out.write("All {}\n".format(uniq[0][0]))
            # Otherwise print the count and proportion of each represented type
            else:
                out.write("{}:\n".format(name))
                for k, v in sorted(uniq, key=lambda x: -x[1]):
                    out.write("  {:<20} {:<10} {:>7.4%}\n".format(k + ":", v, v / N))

        def write_feat(d, text, N, drop=False):
            # If no sequences are of this type (e.g. 'prot'), do nothing
            if N == 0:
                return
            out.write("%s\n" % text)
            # Sort the dictionary by value
            for k, v in sorted(list(d.items()), key=lambda x: -x[1]):
                # If the key is represented, print its count and proportion
                if (drop and v != 0) or not drop:
                    out.write("  {:<20} {:<10} {:>7.4%}\n".format(k + ":", v, v / N))

        write_dict(seqsum.ntype, "Sequence types", nseqs)
        write_dict(seqsum.ncase, "Sequences cases", nseqs)

        nnucl = seqsum.ntype["dna"] + seqsum.ntype["rna"]
        nprot = seqsum.ntype["prot"]
        write_feat(seqsum.nfeat, "Nucleotide Features", nnucl, drop=True)
        write_feat(seqsum.pfeat, "Protein Features:", nprot)
        write_feat(seqsum.ufeat, "Universal Features:", nseqs)


class Stat(Subcommand):
    def _parse(self):
        cmd_name = "stat"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="calculate sequence statistics",
            description="""The default action is to count the lengths of all
            sequences and output summary statistics including: 1) the number of
            sequences, 2) the number of characters, 3) the five-number summary
            of sequence lengths (minimum, 25th quantile, median, 75th quantile,
            and maximum), 4) the mean and standard deviation of lengths, and 5)
            the N50 (if you don't know what that is, you don't need to
            know).""",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument("-d", "--delimiter", help="output delimiter", default="\t")
        parser.add_argument(
            "-q",
            "--byseq",
            help="write a line for each sequence",
            default=False,
            action="store_true",
        )
        parser.add_argument(
            "-I",
            "--case-sensitive",
            help="match case",
            default=False,
            action="store_true",
        )
        parser.add_argument(
            "-m",
            "--count-lower",
            help="count the number of lowercase characters",
            default=False,
            action="store_true",
        )
        parser.add_argument(
            "-c",
            "--counts",
            help="write counts of all characters",
            default=False,
            action="store_true",
        )
        parser.add_argument(
            "-t",
            "--type",
            help="guess sequence type",
            default=False,
            action="store_true",
        )
        parser.add_argument(
            "-l",
            "--length",
            help="write sequence length",
            default=False,
            action="store_true",
        )
        parser.add_argument(
            "-p",
            "--proportion",
            help="write proportion of each character",
            default=False,
            action="store_true",
        )
        parser.add_argument(
            "-C",
            "--aa-profile",
            help="display protein profile",
            default=False,
            action="store_true",
        )
        parser.add_argument(
            "-g",
            "--hist",
            help="write ascii histogram of sequence lengths",
            default=False,
            action="store_true",
        )
        parser.add_argument(
            "-G",
            "--log-hist",
            help="write ascii histogram of sequence log2 lengths",
            default=False,
            action="store_true",
        )
        parser.set_defaults(func=self.func)

    @staticmethod
    def _process_args(args):
        # If no output options are specified, do length stats
        if not any(
            (args.counts, args.type, args.length, args.proportion, args.count_lower)
        ):
            args.length = True
        return args

    @staticmethod
    def _get_length_lines(args, g):
        lines = []
        total = sum(g.lengths)
        N = len(g.lengths)
        if N > 1:
            s = StatFun.summary(g.lengths)

            # Yield total number of sequences
            lines.append("{:10s} {}".format("nseq:", len(g.lengths)))

            # lines.append totla number of letters
            lines.append("{:10s} {}".format("nchars:", sum(g.lengths)))

            # lines.append five number summary of sequence lengths
            fivesum = [
                round(s[x]) for x in ("min", "1st_qu", "median", "3rd_qu", "max")
            ]
            fivesum_str = "{:10s} {} {} {} {} {}"
            lines.append(fivesum_str.format("5sum:", *fivesum))

            # lines.append mean and standard deviation
            meansd_str = "{:10s} {:d} ({:d})"
            lines.append(
                meansd_str.format("mean(sd):", round(s["mean"]), round(s["sd"]))
            )

            # lines.append N50
            lines.append("{:10s} {}".format("N50:", s["N50"]))
        else:
            lstr = ", ".join([str(x) for x in sorted(g.lengths)])
            lines.append("nchars: {}".format(lstr))
        return lines

    @staticmethod
    def _get_hist_lines(args, g, title=None, height=10, width=60, log=False):
        lines = []
        try:
            import numpy
        except ImportError:
            err("Please install numpy (needed for histograms)")

        if title:
            lines.append("")
            lines.append(title)

        if log:
            lengths = [math.log(x, 2) for x in g.lengths]
        else:
            lengths = g.lengths

        y = numpy.histogram(lengths, bins=width)[0]
        y = [height * x / max(y) for x in y]

        for row in reversed(range(height)):
            out = "".join([ascii_histchar(h - row) for h in y])
            lines.append("|{}|".format(out))
        return lines

    @staticmethod
    def _get_aaprofile_lines(args, g, title=None, height=10):
        lines = []
        if title:
            lines.append("")
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
            out = "".join(
                [c + ascii_histchar(y - row, chars=" .:'|") for l, y, c in aacols]
            )
            out = "{}{}".format(out, Colors.OFF)
            lines.append(out)
        names = "".join([l for l, y, c in aacols])
        lines.append(names + Colors.OFF)
        return lines

    @staticmethod
    def _get_count_lines(args, g):
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
            for k, v in count_iter:
                val = v / N if args.proportion else v
                if args.counts:
                    exp = "{}{:>%sd}" % slen
                else:
                    exp = "{}{:>11.5%}"
                lines.append(exp.format(k, val))
        elif args.counts and args.proportion:
            for k, v in count_iter:
                outstr = "{}{:>" + slen + "d}{:>11.5%}"
                lines.append(outstr.format(k, v, v / N))

        if args.count_lower:
            lines.append("{:10s} {} ({:.1%})".format("lower:", lower, lower / N))
        return lines

    def _byfile(self, args, gen):
        g = FileStat()
        # Do I need to count the characters? (much faster if I don't)
        need_count = any(
            (args.counts, args.proportion, args.count_lower, args.type, args.aa_profile)
        )
        for seq in gen.next():
            g.add_seq(SeqStat(seq, count=need_count))

        if need_count:
            lines = self._get_count_lines(args, g)
            yield "\n".join(lines)

        if args.length:
            lines = self._get_length_lines(args, g)
            yield "\n".join(lines)

        if args.hist:
            if args.log_hist:
                lines = self._get_hist_lines(args, g, title="Flat histogram")
            else:
                lines = self._get_hist_lines(args, g)
            yield "\n".join(lines)

        if args.log_hist:
            if args.hist:
                lines = self._get_hist_lines(args, g, title="Log2 histogram", log=True)
            else:
                lines = self._get_hist_lines(args, g, log=True)
            yield "\n".join(lines)

        if args.aa_profile:
            if args.hist or args.log_hist:
                lines = self._get_aaprofile_lines(args, g, title="AA profile")
            else:
                lines = self._get_aaprofile_lines(args, g)
            yield "\n".join(lines)

    @staticmethod
    def _byseq(args, gen):
        seqlist = []
        charset = set()
        if args.length and not (args.counts or args.proportion):
            for seq in gen.next():
                seqid = ParseHeader.firstword(seq.header)
                yield "{}{}{}".format(seqid, args.delimiter, len(seq.seq))
        else:
            for seq in gen.next():
                seqstat = SeqStat(seq)
                seqlist.append(seqstat)
                charset.update(seqstat.counts)

            ignorecase = not args.case_sensitive
            kwargs = {
                "masked": args.count_lower,
                "length": args.length,
                "ignorecase": ignorecase,
            }

            out = SeqStat.getheader(charset, **kwargs)
            yield args.delimiter.join(out)

            count_offset = args.length + args.count_lower + 1
            for q in seqlist:
                line = q.aslist(
                    charset=charset, header_fun=ParseHeader.firstword, **kwargs
                )
                if args.counts:
                    out = [str(x) for x in line]
                elif args.proportion:
                    total = sum(line[count_offset:])
                    props = [c / total for c in line[count_offset:]]
                    out = [str(x) for x in chain(line[0:count_offset], props)]
                yield args.delimiter.join(out)

    def generator(self, args, gen):
        args = self._process_args(args)
        if args.byseq:
            g = self._byseq(args, gen)
        else:
            g = self._byfile(args, gen)
        for item in g:
            yield item


class Split(Subcommand):
    def _parse(self):
        cmd_name = "split"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="split a fasta file into smaller files",
            description="""Breaks a multiple sequence fasta file into several
            smaller files.""",
        )
        parser.add_argument(
            "-n",
            "--number",
            help="Number of output files or sequences per file",
            type=counting_number,
            default=2,
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "-q",
            "--seqs",
            help="split by maximum sequences per file",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-p",
            "--prefix",
            help='prefix for output files (default="xxx")',
            default="xxx",
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        for s in gen.next():
            yield s

    def write(self, args, gen, out=None):
        p = args.prefix
        N = args.number
        used = set()
        for i, seq in enumerate(self.generator(args, gen)):
            fnum = i // N if args.seqs else i % N
            outfile = "%s%s.fasta" % (p, str(fnum))
            if not outfile in used and os.path.isfile(outfile):
                err('Split refuses to overwrite "%s"' % outfile)
            used.add(outfile)
            with open(outfile, "a") as fo:
                fo.write(seq.get_pretty_string() + "\n")


class Subseq(Subcommand):
    def _parse(self):
        cmd_name = "subseq"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="extract subsequence from each entry (revcomp if a<b)",
            description="""The current default action is unfortunately
            excruciating death. The simplest usage is `smof subseq -b START
            STOP`, where START and STOP are two integers. If START is greater
            than STOP, and if the sequence appears to be nucleic, `subseq` will
            write the reverse complement. Subseq can also read start and stop
            positions from a GFF file, where column 1 in the GFF is checked
            against the sequence id (the first word in the fasta header). In
            addition to sequence subsetting, `subseq` can color the matched
            regions.""",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "-b",
            "--bounds",
            metavar="N",
            help="from and to values (indexed from 1)",
            nargs=2,
            type=counting_number,
        )
        parser.add_argument(
            "-f",
            "--gff",
            metavar="FILE",
            help="get bounds from this gff3 file",
            type=argparse.FileType("r"),
        )
        parser.add_argument(
            "-k",
            "--keep",
            help="With --gff, keep sequences with no matches",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-c",
            "--color",
            help="color subsequence (do not extract)",
            choices=Colors.COLORS.keys(),
            metavar="STR",
            default=None,
        )
        parser.add_argument(
            "-Y",
            "--force-color",
            help="print in color even to non-tty (DANGEROUS)",
            action="store_true",
            default=False,
        )
        parser.set_defaults(func=self.func)

    @staticmethod
    def _subseq(seq, a, b, color=None):
        start, end = sorted([a, b])
        end = min(end, len(seq.seq))

        # Check boundaries
        if start > len(seq.seq):
            err("Start position must be less than seq length")

        if color:
            c = Colors.COLORS[color]
            seq.colseq.colorpos(start - 1, end, c)
            outseq = seq
        else:
            outseq = seq.subseq(start - 1, end)
            if (a > b) and seq.get_moltype() == "dna":
                outseq = FSeq.getrevcomp(outseq)
        return outseq

    def _gff_generator(self, args, gen):
        subseqs = defaultdict(list)
        for line in args.gff:
            row = line.split("\t")
            try:
                a, b = int(row[3]), int(row[4])
                if row[6] == "-":
                    subseqs[row[0]].append({"start": max(a, b), "end": min(a, b)})
                else:
                    subseqs[row[0]].append({"start": a, "end": b})
            except IndexError:
                err("Improper gff3 file")
            except ValueError:
                err("gff bounds must be integers")

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
                    seq = self._subseq(seq, s["start"], s["end"], args.color)
                yield seq
            else:
                for s in subseqs[seqid]:
                    yield self._subseq(seq, s["start"], s["end"])

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


class Translate(Subcommand):
    def _parse(self):
        cmd_name = "translate"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="translate a DNA sequence into a protein sequence",
            description="""Only the standard gene code table is supported. Any
            codons with ambiguous characters will be translated as X. Trailing
            characters will be ignored. All gaps [_.-] will be removed. When -f
            is True, then the longest product will be found. Ties are resolved
            by comparing position (earlier positions are preferred) and then
            frame (first frame is preferred). By default, translation starts at
            the first nucleotide.""",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "-s",
            "--from-start",
            help="Require each product begin with a start codon",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-f",
            "--all-frames",
            help="Translate in all frames, keep longest",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-c",
            "--cds",
            help="Write the DNA coding sequence",
            action="store_true",
            default=False,
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        """ Reverse each sequence """
        for seq in gen.next():
            orf = get_orf(
                seq.seq, all_frames=args.all_frames, from_start=args.from_start, translate=not args.cds
            )
            yield FSeq(header=seq.header, seq=orf)


# ==============
# FULL FUNCTIONS
# ==============


class Sample(Subcommand):
    def _parse(self):
        cmd_name = "sample"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="randomly select entries from fasta file",
            description="""Randomly sample entries. `sample` reads the entire
            file into memory, so should not be used for extremely large
            files.""",
        )
        parser.add_argument(
            "-n",
            "--number",
            help="sample size (default=%(default)s)",
            type=counting_number,
            default=1,
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "--seed",
            help="set random seed (for reproducibility/debugging)",
            type=counting_number,
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        """ Randomly sample n entries from input file """
        import random

        if args.seed:
            random.seed(args.seed)
        seqs = [s for s in gen.next()]
        sample_indices = random.sample(range(len(seqs)), min(len(seqs), args.number))
        for i in sample_indices:
            yield seqs[i]


class Sort(Subcommand):
    def _parse(self):
        cmd_name = "sort"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="sort sequences",
            description="""Sorts the entries in a fasta file. By default, it
            sorts by the header strings. `sort` reads the entire file into
            memory, so should not be used for extremely large files.""",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "-x", "--regex", metavar="REG", help="sort by single regex capture"
        )
        parser.add_argument(
            "-r", "--reverse", help="reverse sort", action="store_true", default=False
        )
        parser.add_argument(
            "-n",
            "--numeric-sort",
            help="numeric sort",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-l",
            "--length-sort",
            help="sort by sequence length",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-R",
            "--random-sort",
            help="randomly sort sequences",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-k", "--key", default=None, help="Key to sort on (column number or tag)"
        )
        parser.add_argument(
            "-t",
            "--field-separator",
            help="The field separator (default: '|')",
            default="|",
        )
        parser.add_argument(
            "-p",
            "--pair-separator",
            help="The separator between a tag and value  (default: '=')",
            default="=",
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        seqs = [s for s in gen.next()]

        if args.numeric_sort and not args.regex:
            err("--numeric does nothing unless with --regex")

        # Set type of order determining variable
        if args.numeric_sort:

            def typer(x):
                try:
                    return float(x)
                except ValueError:
                    err("'{}' cannot be numerically sorted".format(x))

        else:

            def typer(x):
                return x

        # Set search term
        if args.regex:
            r = re.compile(args.regex)

            def sortterm(x):
                try:
                    capture = re.search(r, x.header).groups()[0]
                    return typer(capture)
                except AttributeError:
                    err("No match for regex '{}'".format(args.regex))
                except IndexError:
                    err("Nothing was captured in regex '{}'".format(args.regex))

        elif args.key:
            try:
                key = int(args.key) - 1

                def sortterm(x):
                    try:
                        return typer(x.header.split(args.field_separator)[key])
                    except IndexError:
                        err("Cannot sort by column '{}', two few columns".format(key))

            except:
                key = args.key

                def sortterm(x):
                    if args.pair_separator == args.field_separator:
                        xs = x.header.split(args.field_separator)
                        d = {xs[i]: xs[i + 1] for i in range(0, len(xs), 2)}
                    else:
                        xs = [
                            y.split(args.pair_separator)
                            for y in x.header.split(args.field_separator)
                        ]
                        try:
                            d = {k: v for k, v in xs}
                        except ValueError as e:
                            err(str(e))
                    try:
                        return typer(d[key])
                    except KeyError:
                        err("Could not find key '{}'".format(key))
                    except Exception as e:
                        err(str(e))

        elif args.random_sort:
            import random

            def sortterm(x):
                return random.uniform(0, 1)

        elif args.length_sort:

            def sortterm(x):
                return len(x.seq)

        else:

            def sortterm(x):
                return x.header

        seqs.sort(key=sortterm, reverse=args.reverse)

        for s in seqs:
            yield s


class Consensus(Subcommand):
    def _parse(self):
        cmd_name = "consensus"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="finds the consensus sequence for aligned sequence",
            description="""Given input in aligned FASTA file format, where all
            sequences are of equal length (possibly with gaps), `consensus`
            will find the most common character in each column. Ties are
            resolved alphabetically. Optionally, it will instead provide the
            counts or proportions of each character at each position.""",
        )
        parser.add_argument(
            "-t",
            "--table",
            help="Print count table instead of consensus",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        raise NotImplementedError

    def write(self, args, gen, out=sys.stdout):
        seqs = [s for s in gen.next()]
        imax = max([len(s.seq) for s in seqs])
        try:
            transpose = [[s.seq[i] for s in seqs] for i in range(0, imax)]
        except IndexError:
            err("All sequences must be of equivalent length")

        counts = Counter(("").join([s.seq for s in seqs]))
        characters = list(counts.keys())

        if args.table:
            out.write("\t".join(characters))
            out.write("\n")
            for column in transpose:
                c = Counter(column)
                out.write("\t".join([str(c[x]) for x in characters]))
                out.write("\n")

        else:
            consensus = [Counter(c).most_common()[0][0] for c in transpose]
            header = "Consensus"
            FSeq(header, "".join(consensus)).print()


# ==============
# UNIX EMULATORS
# ==============


class Cut(Subcommand):
    def _parse(self):
        cmd_name = "cut"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="emulates UNIX cut command, where fields are entries",
            description="""Prints sequences by index""",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "-f", "--fields", help="Indices to print, comma delimited, with ranges"
        )
        parser.add_argument(
            "-v",
            "--complement",
            help="Invert selection",
            action="store_true",
            default=False,
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        fields = args.fields.split(",")
        indices = set()
        for f in fields:
            try:
                indices.add(int(f) - 1)
            except ValueError:
                try:
                    s, t = f.split("-")
                    try:
                        indices |= set(range(int(s) - 1, int(t)))
                    except TypeError:
                        err(
                            "'{}-{}' does not specify a valid range".format(
                                str(s), str(t)
                            )
                        )
                except ValueError:
                    err("Cannot parse '{}'".format(args.fields))

        i = 0
        if args.complement:
            for seq in gen.next():
                if not i in indices:
                    yield seq
                i += 1
        else:
            m = max(indices)
            for seq in gen.next():
                if i > m:
                    break
                if i in indices:
                    yield seq
                i += 1


class Head(Subcommand):
    def _parse(self):
        cmd_name = "head"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="writes the first sequences in a file",
            description="""`smof head` is modeled after GNU tail and follows
            the same basic conventions except it is entry-based rather than
            line-based. By default, `smof head` outputs ONE sequence (rather
            than the 10 line default for `head`)""",
        )
        parser.add_argument(
            "nseqs", help="-K print first K entries", metavar="K", nargs="?"
        )
        parser.add_argument(
            "-n",
            "--entries",
            metavar="K",
            help="print first K entries; or use -n -K to print all but the last K entries",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "-f",
            "--first",
            help="print first K letters of each sequence",
            metavar="K",
            type=counting_number,
        )
        parser.add_argument(
            "-l",
            "--last",
            help="print last K letters of each sequence",
            metavar="K",
            type=counting_number,
        )

        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        i = 1
        allbut = False

        if args.entries:
            if args.nseqs:
                err("Please don't use nseqs with --entries")
            try:
                if args.entries[0] == "-":
                    allbut = True
                    nseqs = int(args.entries[1:])
                else:
                    nseqs = int(args.entries)
            except AttributeError:
                err("-n (--entries) must be a number")
        elif args.nseqs:
            # This resolve cases where there is a positional filename and no
            # nseqs given.
            # If the first positional argument is a readable filename, treat
            # it as input. Otherwise, try to interpret it as a number
            if os.access(args.nseqs, os.R_OK):
                args.fh = [args.nseqs] + args.fh
                nseqs = 1
            else:
                try:
                    nseqs = int(re.match(r"-(\d+)", args.nseqs).group(1))
                except AttributeError:
                    err("N must be formatted as '-12'")
        else:
            nseqs = 1

        if allbut:
            seqs = list()
            for seq in gen.next():
                if i > nseqs:
                    yield headtailtrunk(seqs.pop(0), args.first, args.last)
                seqs.append(seq)
                i += 1
        else:
            for seq in gen.next():
                yield headtailtrunk(seq, args.first, args.last)
                if i == nseqs:
                    break
                i += 1


class Grep(Subcommand):
    def _parse(self):
        cmd_name = "grep"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="roughly emulates the UNIX grep command",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                """\
                Smof grep is based on GNU grep but operates on fasta entries.
                It allows you to extract entries where either the header or the
                sequence match the search term. For sequence searches, it can
                produce GFF formatted output, which specifies the location of
                each match.

                The --wrap option limits search space to expressions captured
                by a Perl regular expression. This, coupled with the --file
                option, allows thousands of sequences to be rapidly extracted
                based on terms from a file.

                Smof grep can also take a fasta file as a search term input
                (--fastain) and return sequences containing exact matches to
                the sequences in the search fasta file. See the documentation
                for examples.
                """
            ),
        )
        parser.add_argument(
            "pattern", metavar="PATTERN", help="pattern to match", nargs="?"
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "-q",
            "--match-sequence",
            help="match sequence rather than header",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-f",
            "--file",
            metavar="FILE",
            type=argparse.FileType("r"),
            help="obtain patterns from FILE, one per line",
        )
        parser.add_argument(
            "-L",
            "--files-without-match",
            help="print names files with no matches",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-l",
            "--files-with-matches",
            help="print names input files with matches",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-w",
            "--wrap",
            metavar="REG",
            help="a regular expression to capture patterns",
        )
        parser.add_argument(
            "-P",
            "--perl-regexp",
            help="treat patterns as perl regular expressions",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-G",
            "--ambiguous-nucl",
            help="parse extended nucleotide alphabet",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-I",
            "--case-sensitive",
            help="DO NOT ignore case distinctions (ignore by default)",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-v",
            "--invert-match",
            help="print non-matching entries",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-o",
            "--only-matching",
            help="show only the part that matches",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-B",
            "--before-context",
            help="Include N characters before match",
            metavar="N",
            type=positive_int,
            default=0,
        )
        parser.add_argument(
            "-A",
            "--after-context",
            help="Include N characters after match",
            metavar="N",
            type=positive_int,
            default=0,
        )
        parser.add_argument(
            "-C",
            "--context",
            help="Include N characters before and after match",
            metavar="N",
            type=positive_int,
            default=0,
        )
        parser.add_argument(
            "-c",
            "--count",
            help="print number of entries with matches",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-m",
            "--count-matches",
            help="print number of non-overlapping matches",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-x",
            "--line-regexp",
            help="force PATTERN to match the whole text (regex allowed)",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-X",
            "--exact",
            help="force PATTERN to literally equal the text (fast)",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-g",
            "--gapped",
            help="match across gaps when searching aligned sequences",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-b",
            "--both-strands",
            help="search both strands",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-r",
            "--reverse-only",
            help="only search the reverse strand",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-y",
            "--no-color",
            help="do not print in color",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-Y",
            "--force-color",
            help="print in color even to non-tty (DANGEROUS)",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-S",
            "--preserve-color",
            help="Preserve incoming color",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "--color",
            help="Choose a highlight color",
            choices=Colors.COLORS.keys(),
            metavar="STR",
            default="bold_red",
        )
        parser.add_argument(
            "--gff",
            help="output matches in gff format",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "--gff-type",
            help="name of searched feature",
            metavar="STR",
            default="regex_match",
        )
        parser.add_argument(
            "--fastain",
            help="Search for exact sequence matches against FASTA",
            metavar="FASTA",
        )
        parser.set_defaults(func=self.func)

    @staticmethod
    def _process_arguments(args):
        # If the pattern is readable, it is probably meant to be an input, not
        # a pattern
        if args.pattern and os.access(args.pattern, os.R_OK):
            args.fh = [args.pattern] + args.fh
            args.pattern = None

        # Stop if there are any incompatible options
        if args.count_matches and args.invert_match:
            err("--count-matches argument is incompatible with --invert-matches")

        if args.line_regexp and (args.wrap):
            err("--line_regexp is incompatible with --wrap")

        if args.gff and (args.exact or args.line_regexp):
            err("--gff is incompatible with --exact and --line_regexp")

        if args.gff and (args.files_without_match or args.files_with_matches):
            err("--gff is incompatible with -l and -L options")

        if args.fastain:
            args.match_sequence = True

        if not args.match_sequence and (
            args.both_strands or args.ambiguous_nucl or args.reverse_only
        ):
            err("If you want to search sequence, set -q flag")

        if args.wrap and args.perl_regexp:
            err(
                "Patterns found in --wrap captures must be literal (-P and -w incompatible)"
            )

        if args.ambiguous_nucl:
            args.perl_regexp = True

        if args.force_color and args.no_color:
            err("WTF? --force-color AND --no-color?")

        if args.only_matching and (
            args.exact or args.gff or args.count or args.count_matches
        ):
            args.only_matching = False

        if args.only_matching and args.invert_match:
            err("--only-matching is incompatible with --inver-match")

        if (args.perl_regexp or args.ambiguous_nucl) and args.exact:
            err("--exact works only with literal strings (incompatible with -P or -G")

        # Some things just don't make sense in header searches ...
        if args.gff or args.ambiguous_nucl:
            args.match_sequence = True

        # Decide when to color
        if args.force_color or (sys.stdout.isatty() and not args.no_color):
            args.color_output = True
        else:
            args.color_output = False

        args.color = Colors.COLORS[args.color]

        # Others don't make sense with color
        if any(
            (
                args.gff,
                args.count_matches,
                args.only_matching,
                args.line_regexp,
                args.exact,
                args.files_without_match,
                args.files_with_matches,
            )
        ):
            args.color_output = False

        # gff overides certain other options
        if args.gff:
            args.count = False
            args.count_matches = False
        # -A and -B take priority over -C
        args.before_context = (
            args.before_context if args.before_context else args.context
        )
        args.after_context = args.after_context if args.after_context else args.context

        return args

    @staticmethod
    def _create_matcher(args, pat, wrapper):

        # Select a search space preparation function
        if args.match_sequence:
            if args.gapped:
                gettext = lambda x: "".join(re.split(r"-+", x.seq))
            else:
                gettext = lambda x: x.seq
        else:
            gettext = lambda x: x.header

        def context(func):
            has_context = args.before_context or args.after_context
            if has_context:

                def inner(seq, **kw):
                    text = gettext(seq)
                    matches = []
                    for m in func(text, **kw):
                        m["pos"][0] = max(0, m["pos"][0] - args.before_context)
                        m["pos"][1] = min(len(text), m["pos"][1] + args.after_context)
                        matches.append(m)
                    return matches

            else:

                def inner(seq, **kw):
                    text = gettext(seq)
                    return func(text, **kw)

            return inner

        def seqtotext(func):
            def inner(seq, **kw):
                text = gettext(seq)
                return func(text, **kw)

            return inner

        # Check existence for matches to wrapper captures
        @seqtotext
        def swrpmatcher(text, **kw):
            for m in re.finditer(wrapper, text):
                if m.group(1) in pat:
                    return True
            return False

        # Check existence of matches
        @seqtotext
        def spatmatcher(text, **kw):
            for p in pat:
                if re.search(p, text):
                    return True
            return False

        # Check if pattern matches entire text
        @seqtotext
        def linematcher(text, **kw):
            for p in pat:
                m = re.match(p, text)
                if m and m.end() == len(text):
                    return True
            return False

        @seqtotext
        def exactmatcher(text, **kw):
            return text in pat

        @context
        def gwrpmatcher(text, strand="."):
            pos = []
            for m in re.finditer(wrapper, text):
                if m.group(1) in pat:
                    match = {"pos": [m.start(1), m.end(1)], "strand": strand}
                    pos.append(match)
            return pos

        @context
        def gpatmatcher(text, strand="."):
            pos = []
            for p in pat:
                for m in re.finditer(p, text):
                    match = {"pos": [m.start(), m.end()], "strand": strand}
                    pos.append(match)
            return pos

        # the matchers are of two types:
        # 1. boolean - is the pattern present in the given sequence?
        # 2. position - where are the patterns located?
        # this flag records which type of pattern is used
        by_position = False

        # Select a base regular expression function
        if args.exact:
            matcher = exactmatcher
        elif args.line_regexp:
            matcher = linematcher
        elif args.gff or args.count_matches or args.color_output or args.only_matching:
            matcher = gwrpmatcher if wrapper else gpatmatcher
            by_position = True
        else:
            matcher = swrpmatcher if wrapper else spatmatcher

        # Prepare gapped or ungapped search function
        def search_function(matcher, **kw):
            if not args.gapped or not by_position:

                def inner(seq, **kw):  # ungapped
                    return matcher(seq, **kw)

            else:

                def inner(seq, **kw):  # gapped
                    """
                    Maps from positions on an ungapped sequence to positions on a
                    gapped sequence. For example, it can convert the ATA match to
                    'GATACA' on (1,3), to the '-GA--TACA' match on (2,5).
                    """
                    matches = matcher(seq, **kw)
                    gaps = list(re.finditer(r"-+", seq.seq))
                    if not gaps:
                        return matches
                    for g in gaps:
                        g0 = g.span()[0]
                        g1 = g.span()[1]
                        glen = g1 - g0
                        for m in matches:
                            if m["pos"][0] >= g0:
                                m["pos"][0] += glen
                                m["pos"][1] += glen
                        for m in matches:
                            if (m["pos"][0] < g0) and (m["pos"][1] > g0):
                                m["pos"][1] += glen
                    return matches

            return inner

        # Process functions that include reverse complements
        def stranded_function(matcher, by_position):
            if by_position:

                def rev(matcher, seq):
                    rmatch = []
                    text_length = len(seq.seq) if args.gapped else len(gettext(seq))
                    for d in matcher(FSeq.getrevcomp(seq), strand="-"):
                        d["pos"] = text_length - d["pos"][1], text_length - d["pos"][0]
                        rmatch.append(d)
                    return rmatch

                if args.reverse_only:

                    def rmatcher(seq):
                        return rev(matcher, seq)

                else:

                    def rmatcher(seq):
                        fmatch = matcher(seq, strand="+")
                        rmatch = rev(matcher, seq)
                        return fmatch + rmatch

            else:
                f = lambda x: matcher(x, strand="+")
                r = lambda x: matcher(FSeq.getrevcomp(x), strand="-")
                if args.reverse_only:

                    def rmatcher(seq):
                        return r(seq)

                else:

                    def rmatcher(seq):
                        return r(seq) or f(seq)

            return rmatcher

        matcher = search_function(matcher)

        if args.reverse_only or args.both_strands:
            matcher = stranded_function(matcher, by_position)

        return matcher

    @staticmethod
    def _get_pattern(args):
        pat = set()
        if args.fastain:
            # read patterns from a fasta file
            pat.update((s.seq for s in FSeqGenerator(args.fastain).next()))
        if args.file:
            # read patterns from a file (stripping whitespace from the end)
            pat.update([l.rstrip("\n\t\r ") for l in args.file])
        if args.pattern:
            # read pattern from command line
            pat.update([args.pattern])

        if args.ambiguous_nucl:
            apat = set()
            for p in pat:
                perlpat = ambiguous2perl(p)
                apat.update([perlpat])
            pat = apat

        if not pat and not (args.fastain or args.file):
            err("Please provide a pattern")

        # TODO searching for perfect matches would be faster without using
        # regex (just <str>.find(pat))
        if not (args.perl_regexp or args.wrap or args.exact):
            pat = [re.escape(p) for p in pat]

        flags = re.IGNORECASE if not args.case_sensitive else 0

        if args.wrap:
            wrapper = re.compile(args.wrap, flags=flags)
        else:
            wrapper = None

        if not (args.wrap or args.exact):
            pat = set((re.compile(p, flags=flags) for p in pat))

        return (pat, wrapper)

    @staticmethod
    def _makegen(args):

        if args.gff:

            def sgen(gen, matcher):
                source = "smof-{}".format(__version__)
                gfftype = args.gff_type
                row = [
                    None,  # 1 seqid
                    source,  # 2 source
                    gfftype,  # 3 type
                    None,  # 4 start
                    None,  # 5 end
                    ".",  # 6 score
                    None,  # 7 strand
                    ".",  # 8 phase
                    ".",  # 9 attributes
                ]
                for seq in gen.next():
                    row[0] = ParseHeader.firstword(seq.header)
                    matches = list(matcher(seq))
                    for m in matches:
                        row[3] = m["pos"][0] + 1
                        row[4] = m["pos"][1]
                        row[6] = m["strand"]
                        yield "\t".join([str(s) for s in row])

        elif args.count or args.count_matches:

            def sgen(gen, matcher):
                seqcount = OrderedDict()
                seqmatch = OrderedDict()
                for seq in gen.next():
                    m = matcher(seq)
                    if seq.filename not in seqcount:
                        seqcount[seq.filename] = 0
                        seqmatch[seq.filename] = 0
                    if m:
                        try:
                            seqmatch[seq.filename] += len(m)
                            seqcount[seq.filename] += 1
                        except TypeError:
                            pass
                for filename, count in seqcount.items():
                    match = seqmatch[filename]
                    if args.count or args.count_matches:
                        if len(seqcount) > 1:
                            if args.count and args.count_matches:
                                yield "{}\t{}\t{}".format(count, match, filename)
                            elif args.count:
                                yield "{}\t{}".format(count, filename)
                            elif args.count_matches:
                                yield "{}\t{}".format(match, filename)
                        else:
                            if args.count and args.count_matches:
                                yield "{}\t{}".format(count, match)
                            elif args.count:
                                yield count
                            elif args.count_matches:
                                yield match

        elif args.files_without_match or args.files_with_matches:

            def sgen(gen, matcher):
                seqmat = OrderedDict()
                for seq in gen.next():
                    matched = matcher(seq)
                    if seq.filename not in seqmat:
                        seqmat[seq.filename] = False
                    if matched:
                        seqmat[seq.filename] = True
                for filename, matched in seqmat.items():
                    if (not matched and args.files_without_match) or (
                        matched and args.files_with_matches
                    ):
                        yield filename

        elif args.only_matching:

            def sgen(gen, matcher):
                for seq in gen.next():
                    matches = matcher(seq)
                    text = seq.seq if args.match_sequence else seq.header
                    for m in matches:
                        match = text[m["pos"][0] : m["pos"][1]]
                        if args.match_sequence:
                            header = "%s|subseq(%d..%d) %s" % (
                                ParseHeader.firstword(seq.header),
                                m["pos"][0],
                                m["pos"][1],
                                ParseHeader.description(seq.header),
                            )
                            yield FSeq(header, match)
                        else:
                            yield match

        else:

            def sgen(gen, matcher):
                for seq in gen.next(handle_color=args.preserve_color):
                    matches = matcher(seq)
                    if (matches and not args.invert_match) or (
                        not matches and args.invert_match
                    ):
                        if args.color_output:
                            for pos in [m["pos"] for m in matches]:
                                if args.match_sequence:
                                    seq.color_seq(*pos, col=args.color)
                                else:
                                    seq.color_header(*pos, col=args.color)
                        yield seq

        return sgen

    def generator(self, args, gen):
        args = self._process_arguments(args)

        pat, wrapper = self._get_pattern(args)

        matcher = self._create_matcher(args, pat, wrapper)

        sgen = self._makegen(args)

        self.force_color = args.force_color
        for item in sgen(gen, matcher):
            yield item


class Uniq(Subcommand):
    def _parse(self):
        cmd_name = "uniq"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="count, omit, or merge repeated entries",
            description="""Emulates the GNU uniq command. Two entries are
            considered equivalent only if their sequences AND headers are
            exactly equal. Newlines are ignored but all comparisons are
            case-sensitive. The pack/unpack option is designed to be compatible
            with the conventions used in the NCBI-BLAST non-redundant
            databases: entries with identical sequences are merged and their
            headers are joined with SOH (0x01) as a delimiter (by default).""",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        group = parser.add_mutually_exclusive_group()
        group.add_argument(
            "-c",
            "--count",
            help="writes (count|header) in tab-delimited format",
            action="store_true",
            default=False,
        )
        group.add_argument(
            "-d",
            "--repeated",
            help="print only repeated entries",
            action="store_true",
            default=False,
        )
        group.add_argument(
            "-u",
            "--uniq",
            help="print only unique entries",
            action="store_true",
            default=False,
        )
        group.add_argument(
            "-p",
            "--pack",
            help="combine redundant sequences by concatenating the headers",
            action="store_true",
            default=False,
        )
        group.add_argument(
            "-P",
            "--unpack",
            help="reverse the pack operation",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-z",
            "--pack-sep",
            help="set delimiting string for pack/unpack operations (SOH, 0x01, by default)",
            default="\x01",
        )
        parser.add_argument(
            "-f",
            "--final-header",
            help="make headers unique by deleting all but the final entry with a given header (the sequence is ignored, so order matters, you may want to sort by sequence first for reproducibility)",
            action="store_true",
            default=False,
        )
        parser.set_defaults(func=self.func)

    @staticmethod
    def standard_generator(args, gen):

        seqs = OrderedDict()
        for seq in gen.next():
            try:
                seqs[seq] += 1
            except KeyError:
                seqs[seq] = 1

        if args.repeated:
            sgen = ((k, v) for k, v in seqs.items() if v > 1)
        elif args.uniq:
            sgen = ((k, v) for k, v in seqs.items() if v == 1)
        else:
            sgen = seqs.items()

        if args.count:
            for k, v in sgen:
                yield "{}\t{}".format(v, k.header)
        else:
            for k, v in sgen:
                yield k

    @staticmethod
    def pack_generator(args, gen):
        seqs = OrderedDict()
        for seq in gen.next():
            if seq.seq in seqs:
                seqs[seq.seq].append(seq.header)
            else:
                seqs[seq.seq] = [seq.header]
        for q, h in seqs.items():
            seq = FSeq(header=args.pack_sep.join(h), seq=q)
            yield seq

    @staticmethod
    def unpack_generator(args, gen):
        for seq in gen.next():
            headers = seq.header.split(args.pack_sep)
            for header in headers:
                yield FSeq(header=header, seq=seq.seq)

    @staticmethod
    def final_header_generator(args, gen):
        seqs = OrderedDict()
        for seq in gen.next():
            seqs[seq.header] = seq.seq
        for header, sequence in seqs.items():
            seq = FSeq(header=header, seq=sequence)
            yield seq

    def generator(self, args, gen):
        if args.final_header:
            g = self.final_header_generator
        elif args.pack:
            g = self.pack_generator
        elif args.unpack:
            g = self.unpack_generator
        else:
            g = self.standard_generator

        for a in g(args, gen):
            yield a


class Wc(Subcommand):
    def _parse(self):
        cmd_name = "wc"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="roughly emulates the UNIX wc command",
            description="""Outputs the total number of entries and the total
            sequence length (TAB delimited).""",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "-m",
            "--chars",
            help="writes the summed length of all sequences",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "-l",
            "--lines",
            help="writes the total number of sequences",
            action="store_true",
            default=False,
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
            out.write("%s\n" % nchars)
        elif args.lines and not args.chars:
            out.write("%s\n" % nseqs)
        else:
            out.write("{}\t{}\n".format(nseqs, nchars))


class Tail(Subcommand):
    def _parse(self):
        cmd_name = "tail"
        parser = self.subparsers.add_parser(
            cmd_name,
            usage=self.usage.format(cmd_name),
            help="writes the last sequences in a file",
            description="""`smof tail` is modeled after GNU tail and follows
            the same basic conventions except it is entry-based rather than
            line-based. `smof tail` will output ONE sequence (rather than the
            10 line default for `tail`)""",
        )
        parser.add_argument(
            "nseqs", help="-K print last K entries", metavar="K", nargs="?"
        )
        parser.add_argument(
            "-n",
            "--entries",
            metavar="K",
            help="print last K entries; or use -n +K to output starting with the Kth",
        )
        parser.add_argument(
            "fh",
            help="input fasta sequence (default = stdin)",
            metavar="INPUT",
            nargs="*",
        )
        parser.add_argument(
            "-f",
            "--first",
            help="print first K letters of each sequence",
            metavar="K",
            type=counting_number,
        )
        parser.add_argument(
            "-l",
            "--last",
            help="print last K letters of each sequence",
            metavar="K",
            type=counting_number,
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        fromtop = False
        nstring = None

        try:
            # if the first positional is a file, treat is as a file
            # otherwise, it may be for example, '-5', or whatever
            if os.access(args.nseqs, os.R_OK):
                args.fh = [args.nseqs] + args.fh
            else:
                nstring = args.nseqs
        except TypeError:
            pass

        if args.nseqs and args.entries:
            err("ERROR: do not use -n along with positional number argument")

        if not nstring:
            if args.entries:
                nstring = args.entries
            else:
                nstring = "1"

        try:
            m = re.match(r"([-+])(\d+)", nstring)
            if m.group(1) == "+":
                fromtop = True
            nstring = int(m.group(2))
        except AttributeError:
            try:
                nstring = int(nstring)
            except ValueError:
                err("N must be formatted as '[+-]12'")

        if fromtop:
            i = 1
            for seq in gen.next():
                if i >= nstring:
                    yield headtailtrunk(seq, args.first, args.last)
                i += 1
        else:
            from collections import deque

            try:
                lastseqs = deque(maxlen=nstring)
            except ValueError:
                err("--nseqs argument must be positive")
            for seq in gen.next():
                lastseqs.append(seq)

            for s in lastseqs:
                yield headtailtrunk(s, args.first, args.last)


# =======
# EXECUTE
# =======


def main():
    if os.name is not "nt":
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    if os.name is "nt":
        os.system("color")  # allows ANSI color on windows
    args = parse()
    gen = FSeqGenerator(args)
    args.func(args, gen, out=sys.stdout)
