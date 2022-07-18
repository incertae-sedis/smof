import math
import re
import sys
import string
import os
import io
import hashlib
import collections
import itertools
from smof.version import __version__


def to_pair(seq):
    """
    Convert a FastaEntry object to a simple header/sequence pair
    """
    return (seq.header, seq.seq)


def open_fasta(xs):
    """
    Given a single fasta file or a list of fasta files, return a generator that
    will yield individual entries. The returned object is the expected input to
    all fasta processing functions in smof.
    """
    return _stream_entries(xs)


def print_fasta(xs, *args, **kwargs):
    """
    Print all entries in a fasta stream
    """
    for seq in _stream_entries(xs):
        seq.print(*args, **kwargs)


def read_fasta_str(lines, *args, **kwargs):
    """
    text: an iterator of strings
    """
    seq_list = []
    header = None

    for line in lines:
        line = line.strip()
        if line == "" or line[0] == "#":
            continue
        if line[0] == ">":
            if seq_list:
                yield FastaEntry(header, "".join(seq_list), *args, **kwargs)
            elif header:
                # NOTE: yields an empty sequence! This is usually
                # a BAD THING, but it can happen in the wild
                yield FastaEntry(header, "", *args, **kwargs)
            seq_list = []
            header = line[1:]
        # '' is valid for a header
        elif header is not None:
            seq_list.append(line)
        else:
            _err("First fasta line must begin with '>'")

    # process the last sequence
    if header is not None:
        if seq_list:
            yield FastaEntry(header, "".join(seq_list), *args, **kwargs)
        else:
            # NOTE: yields empty sequence!
            yield FastaEntry(header, "", *args, **kwargs)


def read_fasta(fastafile, *args, **kwargs):
    """
    fastafile may be a filename or a file object
    """
    if isinstance(fastafile, str):
        f = open(fastafile, "r")
    else:
        f = fastafile

    for seq in read_fasta_str(f, *args, **kwargs):
        yield seq

    f.close()


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
    delimiter="[ |]",
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
            _err("Type not recognized")

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

    for seq in gen:

        if reduce_header:
            seq.header = _parse_header_firstword(seq.header, delimiter=delimiter)

        if standardize:
            try:
                seq.seq = seq.seq.translate(standard_trans)
            except UnboundLocalError:
                _err("Please provide a type argument (-t)")

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


def cut(gen, indices, complement=False):
    i = 0
    if complement:
        for seq in gen:
            if not i in indices:
                yield seq
            i += 1
    else:
        m = max(indices)
        for seq in gen:
            if i > m:
                break
            if i in indices:
                yield seq
            i += 1


def consensus(gen, table=False):
    seqs = [s for s in gen]
    imax = max([len(s.seq) for s in seqs])
    try:
        transpose = [[s.seq[i] for s in seqs] for i in range(0, imax)]
    except IndexError:
        _err("All sequences must be of equivalent length")

    characters = list(set(("").join([s.seq for s in seqs])))

    if table:
        rows = []
        for column in transpose:
            c = collections.Counter(column)
            rows.append([c[x] for x in characters])
        return (characters, rows)
    else:
        consensus = [collections.Counter(c).most_common()[0][0] for c in transpose]
        header = "Consensus"
        return FastaEntry(header, "".join(consensus))


class GrepOptions:
    def __init__(
        self,
        pattern=None,
        match_sequence=False,
        file=None,
        files_without_match=False,
        files_with_matches=False,
        wrap=None,
        perl_regexp=False,
        ambiguous_nucl=False,
        case_sensitive=False,
        invert_match=False,
        only_matching=False,
        before_context=False,
        after_context=False,
        context=False,
        count=False,
        count_matches=False,
        line_regexp=False,
        exact=False,
        gapped=False,
        both_strands=False,
        reverse_only=False,
        no_color=False,
        force_color=False,
        preserve_color=False,
        color="bold_red",
        gff=False,
        gff_type="regex_match",
        fastain=None,
    ):
        self.pattern = pattern
        self.match_sequence = match_sequence
        self.file = file
        self.files_without_match = files_without_match
        self.files_with_matches = files_with_matches
        self.wrap = wrap
        self.perl_regexp = perl_regexp
        self.ambiguous_nucl = ambiguous_nucl
        self.case_sensitive = case_sensitive
        self.invert_match = invert_match
        self.only_matching = only_matching
        self.before_context = before_context
        self.after_context = after_context
        self.context = context
        self.count = count
        self.count_matches = count_matches
        self.line_regexp = line_regexp
        self.exact = exact
        self.gapped = gapped
        self.both_strands = both_strands
        self.reverse_only = reverse_only
        self.no_color = no_color
        self.force_color = force_color
        self.preserve_color = preserve_color
        self.color = color
        self.gff = gff
        self.gff_type = gff_type
        self.fastain = fastain


def grep(gen, **kwargs):
    args = GrepOptions(**kwargs)
    src = GrepSearch(args)
    return src.search(gen)


def head(gen, nseqs, first, last, allbut):
    i = 1
    if allbut:
        seqs = list()
        for seq in gen:
            if i > nseqs:
                yield _headtailtrunk(seqs.pop(0), first, last)
            seqs.append(seq)
            i += 1
    else:
        for seq in gen:
            yield _headtailtrunk(seq, first, last)
            if i == nseqs:
                break
            i += 1


def permute(gen, seed=42, word_size=1, start_offset=0, end_offset=0):
    import random

    if seed:
        random.seed(seed)
    w = word_size
    start = start_offset
    end = end_offset
    for seq in gen:
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
        header = _parse_header_permute(seq.header, start, end, w)
        yield FastaEntry(header, out)


def reverse(gen, complement=False, no_validate=False):
    """ Reverse each sequence """
    if complement:

        def f(s):
            return FastaEntry.getrevcomp(s)

    else:

        def f(s):
            s.reverse()
            return s

    if complement and not no_validate:

        def func(s):
            if s.get_moltype() == "dna":
                return f(s)
            else:
                msg = "Cannot take reverse complement of the sequence '%s' since it does not appear to DNA"
                _err(msg % _parse_header_firstword(s.header))

    else:
        func = f

    for seq in gen:
        yield func(seq)


def sniff(gen):
    seqsum = FastaDescription()
    for seq in gen:
        seqsum.add_seq(seq)
    return seqsum


def stat_seq(
    gen,
    length=False,
    counts=False,
    proportion=False,
    case_sensitive=False,
    count_lower=False,
):
    if not (counts or proportion or case_sensitive or count_lower):
        length = True

    seqlist = []
    charset = set()
    if length and not (counts or proportion):
        for seq in gen:
            seqid = _parse_header_firstword(seq.header)
            yield [seqid, len(seq.seq)]
    else:
        for seq in gen:
            seqstat = FastaEntryStat(seq)
            seqlist.append(seqstat)
            charset.update(seqstat.counts)

        ignorecase = not case_sensitive
        kwargs = {"masked": count_lower, "length": length, "ignorecase": ignorecase}

        yield FastaEntryStat.getheader(charset, **kwargs)

        count_offset = length + count_lower + 1
        for q in seqlist:
            line = q.aslist(
                charset=charset, header_fun=_parse_header_firstword, **kwargs
            )
            if counts:
                out = line
            elif proportion:
                total = sum(line[count_offset:])
                props = [c / total for c in line[count_offset:]]
                out = itertools.chain(line[0:count_offset], props)
            yield out


def stat_file(gen, count_characters=False):
    g = FastaStat()
    for seq in gen:
        g.add_seq(FastaEntryStat(seq, count=count_characters))
    return g


def subseq(gen, a, b, color=None, annotate=False):
    for seq in gen:
        start, end = sorted([a, b])
        end = min(end, len(seq.seq))

        # Check boundaries
        if start > len(seq.seq):
            _err("Start position must be less than seq length")

        if color:
            c = Colors.COLORS[color]
            seq.color_seq(start - 1, end, c)
            outseq = seq
        else:
            outseq = seq.subseq(start - 1, end, annotate=annotate)
            if (a > b) and seq.get_moltype() == "dna":
                outseq = FastaEntry.getrevcomp(outseq)
        yield outseq


def gff_subseq(gen, gff_file, keep=False, color=None):
    subseqs = defaultdict(list)
    for line in gff_file:
        row = line.split("\t")
        try:
            a, b = int(row[3]), int(row[4])
            if row[6] == "-":
                subseqs[row[0]].append({"start": max(a, b), "end": min(a, b)})
            else:
                subseqs[row[0]].append({"start": a, "end": b})
        except IndexError:
            _err("Improper gff3 file")
        except ValueError:
            _err("gff bounds must be integers")

    for seq in gen:
        seqid = _parse_header_firstword(seq.header)
        try:
            if seqid not in subseqs.keys():
                raise KeyError
        except KeyError:
            if keep:
                yield seq
            continue

        if color:
            for s in subseqs[seqid]:
                seq = subseq(seq, s["start"], s["end"], color)
            yield seq
        else:
            for s in subseqs[seqid]:
                yield subseq(seq, s["start"], s["end"])


def find_max_orf(dna, from_start=False):
    dna = dna.translate(FastaEntry.ungapper).upper()
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


def _get_orf(dna, all_frames=False, from_start=False, translate=True):
    dna = dna.translate(FastaEntry.ungapper).upper()
    if all_frames:
        cds_start, cds_length = find_max_orf(dna, from_start=from_start)
        if cds_start is None:
            return ""
        cds = dna[cds_start : cds_start + cds_length]
        if translate:
            return translate_dna(cds)
        else:
            return cds
    else:
        aa = translate_dna(dna)
        if from_start and not aa[0] == "M":
            aa = ""
        return aa


def translate(gen, all_frames=False, from_start=False, cds=False):
    for seq in gen:
        orf = _get_orf(
            seq.seq, all_frames=all_frames, from_start=from_start, translate=not cds
        )
        yield FastaEntry(header=seq.header, seq=orf)


def uniq(gen, repeated=False, uniq=False, count=False):
    seqs = collections.OrderedDict()
    for seq in gen:
        try:
            seqs[seq] += 1
        except KeyError:
            seqs[seq] = 1

    if repeated:
        sgen = ((k, v) for k, v in seqs.items() if v > 1)
    elif uniq:
        sgen = ((k, v) for k, v in seqs.items() if v == 1)
    else:
        sgen = seqs.items()

    if count:
        for k, v in sgen:
            yield "{}\t{}".format(v, k.header)
    else:
        for k, v in sgen:
            yield k


def pack(gen, sep):
    seqs = collections.OrderedDict()
    for seq in gen:
        if seq.seq in seqs:
            seqs[seq.seq].append(seq.header)
        else:
            seqs[seq.seq] = [seq.header]
    for q, h in seqs.items():
        seq = FastaEntry(header=sep.join(h), seq=q)
        yield seq


def unpack(gen, sep):
    for seq in gen:
        headers = seq.header.split(sep)
        for header in headers:
            yield FastaEntry(header=header, seq=seq.seq)


def uniq_headers(gen, removed=False):
    seqs = collections.OrderedDict()
    for seq in gen:
        if seq.header in seqs:
            if removed:
                seq.print(color=False, out=args.removed)
        else:
            seqs[seq.header] = seq.seq
    for header, sequence in seqs.items():
        seq = FastaEntry(header=header, seq=sequence)
        yield seq


# =================
# UTILITY FUNCTIONS
# =================


def _counter_caser(counter, lower=False):
    """
    Sums cases in Collections.Counter object
    """
    if lower:
        out = counter + collections.Counter(
            {k.lower(): v for k, v in counter.items() if k.isupper()}
        )
        out = out - collections.Counter(
            {k: v for k, v in counter.items() if k.isupper()}
        )
    else:
        out = counter + collections.Counter(
            {k.upper(): v for k, v in counter.items() if k.islower()}
        )
        out = out - collections.Counter(
            {k: v for k, v in counter.items() if k.islower()}
        )
    return out


def _sum_lower(counter):
    lc = [v for k, v in counter.items() if k in string.ascii_lowercase]
    return sum(lc)


def _guess_type(counts):
    """
    Predict sequence type from character counts (dna|rna|prot|ambiguous|illegal)
    """
    if isinstance(counts, str):
        counts = collections.Counter(counts)
    elif isinstance(counts, FastaEntry):
        counts = collections.Counter(counts.seq)

    # Convert all to upper case
    counts = _counter_caser(counts)
    # Remove gaps from Counter
    counts = collections.Counter(
        {k: n for k, n in counts.items() if k not in Alphabet.GAP}
    )

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


def _headtailtrunk(seq, first=None, last=None):
    """
    This function is used by the Head and Tail classes to portray partial
    sections of sequences.
    """
    outseq = FastaEntry(seq.header, seq.seq)
    if first and last:
        if first + last < len(seq.seq):
            outseq.header = _parse_header_firstword(
                seq.header
            ) + "|TRUNCATED:first-{}_last-{}".format(first, last)
            outseq.seq = "{}{}{}".format(seq.seq[0:first], "...", seq.seq[-last:])
    elif first:
        outseq.header = _parse_header_firstword(
            seq.header
        ) + "|TRUNCATED:first-{}".format(first)
        outseq.seq = seq.seq[0:first]
    elif last:
        outseq.header = _parse_header_firstword(
            seq.header
        ) + "|TRUNCATED:last-{}".format(last)
        outseq.seq = seq.seq[-last:]
    elif first == 0 and last == 0:
        _err("Illegal empty sequence, dying ...")
    return outseq


def _ascii_histchar(dif, chars=" .~*O"):
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


def _err(msg):
    sys.exit(msg)


def ambiguous2perl(pattern):
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
    perlpat = []
    in_bracket = False
    escaped = False
    for c in pattern:
        amb = c in DNA_AMB
        if c == "\\":
            escaped = True
            continue
        elif escaped:
            c = c if amb else "\\" + c
            escaped = False
        elif amb:
            v = DNA_AMB[c]
            c = v if in_bracket else "[%s]" % v
        elif c == "[":
            in_bracket = True
        elif c == "]":
            in_bracket = False
        perlpat.append(c)
    return "".join(perlpat)


def translate_dna(dna):
    # remove gaps
    dna = dna.translate(FastaEntry.ungapper).upper()
    aa = []
    for i in range(2, len(dna), 3):
        codon = dna[i - 2 : i + 1]
        if codon in FastaEntry.codon_table:
            aa.append(FastaEntry.codon_table[codon])
        else:
            aa.append("X")
    return "".join(aa)


def _parse_header_firstword(h, delimiter="[ \t]"):
    return re.sub("%s.*" % delimiter, "", h)


def _parse_header_description(h):
    return re.sub(r"^\S+\s*", "", h)


def _parse_header_add_suffix(h, suffix):
    return re.sub(r"^(\S+)(.*)", "\\1|%s\\2" % suffix, h)


def _parse_header_add_tag(h, tag, value):
    return re.sub(r"^(\S+)(.*)", "\\1 %s=%s\\2" % (tag, value), h)


def _parse_header_subseq(h, a, b):
    header = "%s|subseq(%d..%d) %s" % (
        _parse_header_firstword(h),
        a,
        b,
        _parse_header_description(h),
    )
    return header.strip()


def _parse_header_permute(h, start, end, wordsize):
    header = "%s|permutation:start=%d;end=%d;word_size=%d %s" % (
        _parse_header_firstword(h),
        start,
        end,
        wordsize,
        _parse_header_description(h),
    )
    return header.strip()


def _parse_header_ncbi_format(h, fields):
    raise NotImplementedError


def _parse_header_regex_group(h, regex):
    raise NotImplementedError


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
        if isinstance(thing, FastaEntry):
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
            _err(
                "ColorString can only append strings, FastaEntry, or ColorString objects"
            )
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


class FastaDescription:
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
        @type seq: FastaEntry object
        """
        # Add md5 hashes of sequences and headers to respective sets
        self.seqs.update([hashlib.md5(bytes(seq.seq, "ascii")).digest()])
        self.headers.update([hashlib.md5(bytes(seq.header, "ascii")).digest()])

        counts = collections.Counter(seq.seq)

        # Ungapped the sequence is required for downstream analysis
        is_gapped = self._handle_gaps(seq, counts)
        if is_gapped:
            seq.ungap()
            counts = collections.Counter(seq.seq)

        scase = self._handle_case(counts)
        # Sum upper and lowercase
        if scase not in "upper":
            counts = _counter_caser(counts)

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
        stype = _guess_type(counts)
        self.ntype[stype] += 1
        return stype

    def __str__(self):
        # Total number of sequences
        nseqs = self.get_nseqs()

        result = []

        # Print number of uniq and total sequences
        if self.count_degenerate_seqs():
            uniq = nseqs - self.count_degenerate_seqs()
            result.append("{} uniq sequences ({} total)".format(uniq, nseqs))
        else:
            result.append("Total sequences: {}".format(nseqs))

        # Warn if there are any duplicate headers
        if self.count_degenerate_headers():
            uniq = nseqs - self.count_degenerate_headers()
            result.append("WARNING: headers are not unique ({}/{})".format(uniq, nseqs))

        # Warn if there are any illegal characters
        if self.ntype["illegal"]:
            result.append("WARNING: illegal characters found")

        def write_dict(d, name, N):
            # Print keys if value is greater than 0
            uniq = [[k, v] for k, v in d.items() if v > 0]
            # E.g. If all of the sequences are proteins, print 'All prot'
            if len(uniq) == 1:
                result.append("All {}".format(uniq[0][0]))
            # Otherwise print the count and proportion of each represented type
            else:
                result.append("{}:".format(name))
                for k, v in sorted(uniq, key=lambda x: -x[1]):
                    result.append("  {:<20} {:<10} {:>7.4%}".format(k + ":", v, v / N))

        def write_feat(d, text, N, drop=False):
            # If no sequences are of this type (e.g. 'prot'), do nothing
            if N == 0:
                return
            result.append("%s" % text)
            # Sort the dictionary by value
            for k, v in sorted(list(d.items()), key=lambda x: -x[1]):
                # If the key is represented, print its count and proportion
                if (drop and v != 0) or not drop:
                    result.append("  {:<20} {:<10} {:>7.4%}".format(k + ":", v, v / N))

        write_dict(self.ntype, "Sequence types", nseqs)
        write_dict(self.ncase, "Sequences cases", nseqs)

        nnucl = self.ntype["dna"] + self.ntype["rna"]
        nprot = self.ntype["prot"]
        write_feat(self.nfeat, "Nucleotide Features", nnucl, drop=True)
        write_feat(self.pfeat, "Protein Features:", nprot)
        write_feat(self.ufeat, "Universal Features:", nseqs)

        result.append("")

        return "\n".join(result)


class FastaStat:
    def __init__(self):
        self.counts = collections.Counter()
        self.nseqs = 0
        self.lengths = []

    def add_seq(self, stat):
        if stat.counts:
            self.counts += stat.counts
        self.nseqs += 1
        self.lengths.append(stat.length)

    def get_length(self):
        lines = []
        total = sum(self.lengths)
        N = len(self.lengths)
        if N > 1:
            s = _summary(self.lengths)

            # Yield total number of sequences
            lines.append("{:10s} {}".format("nseq:", len(self.lengths)))

            # lines.append totla number of letters
            lines.append("{:10s} {}".format("nchars:", sum(self.lengths)))

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
            lstr = ", ".join([str(x) for x in sorted(self.lengths)])
            lines.append("nchars: {}".format(lstr))
        return "\n".join(lines)

    def get_hist(self, title=None, height=10, width=60, log=False):
        lines = []
        try:
            import numpy
        except ImportError:
            _err("Please install numpy (needed for histograms)")

        if title:
            lines.append("")
            lines.append(title)

        if log:
            lengths = [math.log(x, 2) for x in self.lengths]
        else:
            lengths = self.lengths

        y = numpy.histogram(lengths, bins=width)[0]
        y = [height * x / max(y) for x in y]

        for row in reversed(range(height)):
            out = "".join([_ascii_histchar(h - row) for h in y])
            lines.append("|{}|".format(out))
        return "\n".join(lines)

    def get_aaprofile(self, title=None, height=10, case_sensitive=False):
        lines = []
        if title:
            lines.append("")
            lines.append(title)

        colorAA = ColorAA()
        aacols = []
        for chars, group, color in colorAA.group:
            for c in chars:
                if not case_sensitive and c.islower():
                    continue
                cheight = height * self.counts[c] / max(self.counts.values())
                aacols.append([c, cheight, color])
        # Draw histogram
        for row in reversed(range(height)):
            out = "".join(
                [c + _ascii_histchar(y - row, chars=" .:'|") for l, y, c in aacols]
            )
            out = "{}{}".format(out, Colors.OFF)
            lines.append(out)
        names = "".join([l for l, y, c in aacols])
        lines.append(names + Colors.OFF)
        return "\n".join(lines)

    def get_count(
        self,
        count_lower=False,
        case_sensitive=False,
        type=False,
        counts=False,
        proportion=False,
    ):
        lines = []
        lower = _sum_lower(self.counts) if count_lower else None
        if not case_sensitive:
            self.counts = _counter_caser(self.counts)

        if type:
            lines.append(_guess_type(self.counts))

        N = sum(self.lengths)
        slen = str(len(str(max(self.counts.values()))) + 2)
        count_iter = sorted(self.counts.items(), key=lambda x: -x[1])
        if counts ^ proportion:
            for k, v in count_iter:
                val = v / N if proportion else v
                if counts:
                    exp = "{}{:>%sd}" % slen
                else:
                    exp = "{}{:>11.5%}"
                lines.append(exp.format(k, val))
        elif counts and proportion:
            for k, v in count_iter:
                outstr = "{}{:>" + slen + "d}{:>11.5%}"
                lines.append(outstr.format(k, v, v / N))

        if count_lower:
            lines.append("{:10s} {} ({:.1%})".format("lower:", lower, lower / N))
        return "\n".join(lines)


class FastaEntry:
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
        self.seq = self.seq.translate(FastaEntry.ungapper)
        self.header = _parse_header_add_suffix(self.header, "ungapped")

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
            self.moltype = _guess_type(self)

    def get_moltype(self):
        if not self.moltype:
            self.set_moltype()
        return self.moltype

    def header_upper(self):
        self.header = self.header.upper()

    def subseq(self, a, b, annotate=True):
        if annotate:
            header = _parse_header_subseq(self.header, a + 1, b)
        else:
            header = self.header
        newseq = FastaEntry(header, self.seq[a:b])
        if self.colseq:
            newseq.colseq = self.colseq.copy()
            newseq.colseq.subseq(a, b)
        return newseq

    def add_filename(self):
        if self.filename:
            self.header = _parse_header_add_tag(
                h=self.header, tag="filename", value=self.filename
            )
            self.colheader = None

    def reverse(self):
        self.seq = self.seq[::-1]
        if self.handle_color:
            self.colseq.reverse(len(self.seq))
        self.header = _parse_header_add_suffix(self.header, "reverse")

    @classmethod
    def getrevcomp(cls, seq):
        trans = lambda s: s[::-1].translate(FastaEntry.revtrans)
        if isinstance(seq, str):
            return trans(seq)
        elif isinstance(seq, FastaEntry):
            newheader = _parse_header_add_suffix(seq.header, "revcom")
            newseq = FastaEntry(newheader, trans(seq.seq))
            if seq.colseq:
                newseq.colseq = seq.colseq
                newseq.colseq.reverse(len(seq.seq))
            if seq.colheader:
                # TODO implement this
                pass
            return newseq


def _stream_entries(entries, *args, **kwargs):
    if (
        not hasattr(entries, "__iter__")
        or isinstance(entries, str)
        or isinstance(entries, io.TextIOWrapper)
    ):
        entries = [entries]
    else:
        entries = entries

    for entry in entries:
        if isinstance(entry, tuple):
            # if this is a pair, create a FastaEntry object
            yield FastaEntry(entry[0], entry[1], filename=None, *args, *kwargs)

        # if this is a single entry, yield it
        elif isinstance(entry, FastaEntry):
            yield entry

        # maybe it is a fasta file?
        elif isinstance(entry, str) or isinstance(entry, io.TextIOWrapper):
            for seq in read_fasta(entry, *args, **kwargs):
                yield seq

        else:
            print("Can't handle this type:" + str(type(entry)), file=sys.stderr)
            raise ShitInput


class FastaEntryStat:
    def __init__(self, seq, count=True):
        self.counts = collections.Counter(seq.seq) if count else None
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
                _err("Cannot process header: '{}'".format(self.header))
        if length:
            line.append(self.length)

        if masked:
            line.append(_sum_lower(self.counts))

        if ignorecase:
            charset = set("".join(charset).upper())
            self.counts = _counter_caser(self.counts)

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


def _N50(xs, issorted=False):
    xs = sorted(xs) if not issorted else xs
    N = sum(xs)
    total = 0
    for i in range(len(xs) - 1, -1, -1):
        total += xs[i]
        if total > N / 2:
            return xs[i]


def _mean(xs):
    if not xs:
        mu = float("nan")
    else:
        mu = sum(xs) / len(xs)
    return mu


def _median(xs, issorted=False):
    return _quantile(xs, 0.5, issorted=issorted)


def _sd(xs):
    if len(xs) < 2:
        stdev = float("nan")
    else:
        mean = sum(xs) / len(xs)
        stdev = (sum((y - mean) ** 2 for y in xs) / (len(xs) - 1)) ** 0.5
    return stdev


def _quantile(xs, q, issorted=False):
    """
    Calculates quantile as the weighted average between indices
    """
    # Die if out of bounds
    if not 0 <= q <= 1:
        _err("quantile must be between 0 and 1")

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


def _summary(xs):
    xs = sorted(xs)
    out = {
        "min": xs[0],
        "max": xs[-1],
        "1st_qu": _quantile(xs, 0.25, issorted=True),
        "median": _quantile(xs, 0.50, issorted=True),
        "3rd_qu": _quantile(xs, 0.75, issorted=True),
        "mean": _mean(xs),
        "sd": _sd(xs),
        "N50": _N50(xs, issorted=True),
    }
    return out


class GrepSearch:
    def __init__(self, args):
        self.clean_args = self._process_arguments(args)
        pat, wrapper = self._get_pattern(self.clean_args)
        self.matcher = self._create_matcher(self.clean_args, pat, wrapper)
        self.generator = self._makegen(self.clean_args)

    def search(self, gen):
        for item in self.generator(gen, self.matcher):
            yield item

    @staticmethod
    def _process_arguments(args):
        # Stop if there are any incompatible options
        if args.count_matches and args.invert_match:
            _err("--count-matches argument is incompatible with --invert-matches")

        if args.line_regexp and args.wrap:
            _err("--line_regexp is incompatible with --wrap")

        if args.gff and (args.exact or args.line_regexp):
            _err("--gff is incompatible with --exact and --line_regexp")

        if args.gff and (args.files_without_match or args.files_with_matches):
            _err("--gff is incompatible with -l and -L options")

        if args.fastain:
            args.match_sequence = True

        if not args.match_sequence and (
            args.both_strands or args.ambiguous_nucl or args.reverse_only
        ):
            _err("If you want to search sequence, set -q flag")

        if args.wrap and args.perl_regexp:
            _err(
                "Patterns found in --wrap captures must be literal (-P and -w incompatible)"
            )

        if args.ambiguous_nucl:
            args.perl_regexp = True

        if args.force_color and args.no_color:
            _err("WTF? --force-color AND --no-color?")

        if args.only_matching and (
            args.exact or args.gff or args.count or args.count_matches
        ):
            args.only_matching = False

        if args.only_matching and args.invert_match:
            _err("--only-matching is incompatible with --inver-match")

        if (args.perl_regexp or args.ambiguous_nucl) and args.exact:
            _err("--exact works only with literal strings (incompatible with -P or -G")

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
    def _get_pattern(args):
        pat = set()
        if args.fastain:
            # read patterns from a fasta file
            pat.update((s.seq for s in read_fasta(args.fastain)))
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
            _err("Please provide a pattern")

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
                    for d in matcher(FastaEntry.getrevcomp(seq), strand="-"):
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
                r = lambda x: matcher(FastaEntry.getrevcomp(x), strand="-")
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
                for seq in gen:
                    row[0] = _parse_header_firstword(seq.header)
                    matches = list(matcher(seq))
                    for m in matches:
                        row[3] = m["pos"][0] + 1
                        row[4] = m["pos"][1]
                        row[6] = m["strand"]
                        yield "\t".join([str(s) for s in row])

        elif args.count or args.count_matches:

            def sgen(gen, matcher):
                seqcount = collections.OrderedDict()
                seqmatch = collections.OrderedDict()
                for seq in gen:
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
            print("C")

            def sgen(gen, matcher):
                seqmat = collections.OrderedDict()
                for seq in gen:
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
                for seq in gen:
                    matches = matcher(seq)
                    text = seq.seq if args.match_sequence else seq.header
                    for m in matches:
                        match = text[m["pos"][0] : m["pos"][1]]
                        if args.match_sequence:
                            header = "%s|subseq(%d..%d) %s" % (
                                _parse_header_firstword(seq.header),
                                m["pos"][0],
                                m["pos"][1],
                                _parse_header_description(seq.header),
                            )
                            yield FastaEntry(header, match)
                        else:
                            yield match

        else:

            def sgen(gen, matcher):
                for seq in gen:
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


def md5sum(
    gen,
    ignore_case=False,
    each_sequence=False,
    all_sequences=False,
    all_headers=False,
    replace_header=False,
):
    md5hash = hashlib.md5()
    # Hash the sequences only (in input order)
    if all_sequences:
        fun = lambda s, h: md5hash.update(s)
    # Hash the headers only (in input order)
    elif all_headers:
        fun = lambda s, h: md5hash.update(h)
    # DEFAULT: Hash headers and sequences (concatenated)
    # Equivalent to:
    # $ tr -d '\n>' < myfile.fa | md5sum
    else:
        fun = lambda s, h: md5hash.update(h + s)

    for seq in gen:
        if ignore_case:
            seq.header_upper()
            seq.seq_upper()
        s = seq.seq.encode("ascii")
        h = seq.header.encode("ascii")
        # Write <header>\t<sequence hash> for each sequence
        if replace_header:
            yield FastaEntry(hashlib.md5(s).hexdigest(), seq.seq)
        elif each_sequence:
            yield "{}\t{}".format(
                _parse_header_firstword(seq.header), hashlib.md5(s).hexdigest()
            )
        else:
            fun(s, h)

    # Print output hash for cumulative options
    if not (each_sequence or replace_header):
        yield md5hash.hexdigest()
