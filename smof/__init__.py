import math
import re
import sys
import string
from hashlib import md5
from collections import Counter

# ========
# Commands
# ========


def clean(
    gen,
    seq_type=None,
    toupper=False,
    tolower=False,
    toseq=False,
    standardize=False,
    reduce_header=False,
    mask_lowercase=False,
    mask_irregular=False,
):

    # Make translator, this is only possible if type is given
    trans = None
    if seq_type:
        if seq_type.lower()[0] in ["a", "p"]:
            # irregulars proteins include selenocysteine (U) and protein
            # ambiguous characters
            irr = "".join(Alphabet.PROT_AMB) + "U"
            unk = "X"
            standard_trans = str.maketrans("._", "--")
        elif seq_type.lower()[0] in ["n", "d"]:
            # irregular nucleotides are ambiguous characters
            irr = "".join(Alphabet.DNA_AMB)
            unk = "N"
            standard_trans = str.maketrans("Xx._", "Nn--")
        else:
            err("Type not recognized")

        a = ""
        # Get irregular characters
        if mask_irregular:
            a += irr.upper() + irr.lower()

        # convert lowercase to unknown
        if mask_lowercase:
            a += string.ascii_lowercase

        a = "".join(set(a))

        b = unk * len(a)
        if irr:
            trans = str.maketrans(a, b)

    for seq in gen.next(purge_color=True):
        if reduce_header:
            seq.header = ParseHeader.firstword(seq.header, delimiter=" \t|")

        if standardize:
            try:
                seq.seq = seq.seq.translate(standard_trans)
            except UnboundLocalError:
                err("Please provide a type argument (-t)")

        # WARNING: order is important here, don't swap thoughtlesly
        # Remove all nonletters or wanted, otherwise just remove space
        if toseq:
            seq.seq = re.sub(r"[^A-Za-z]", "", seq.seq)
        else:
            seq.seq = re.sub(r"[^\S]", "", seq.seq)

        # Irregular or lowercase to unknown
        if trans:
            seq.seq = seq.seq.translate(trans)

        # Change everything to desired case
        if toupper:
            seq.seq = seq.seq.upper()
        elif tolower:
            seq.seq = seq.seq.lower()

        yield seq


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
    return "".join(aa)


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
