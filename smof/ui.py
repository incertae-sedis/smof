import argparse
import math
import re
import sys
import string
import os
import signal
import textwrap
from collections import Counter
from collections import defaultdict
from collections import OrderedDict

from smof.functions import *
from smof.functions import _headtailtrunk
from smof.functions import _stream_entries
from smof.functions import _err
from smof.version import __version__

# =============
# argpase types
# =============


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
        _err("{} is deprecated".format(argv[0]))
    if argv[0] in ["idsearch", "retrieve", "search", "rmfields"]:
        _err("{} is deprecated, use 'smof grep'".format(argv[0]))
    if argv[0] == "winnow":
        _err("`winnow` is deprecated, use `smof filter`")
    if argv[0] == "chksum":
        _err("`winnow` is deprecated, use `smof md5sum`")
    if argv[0] == "perm":
        _err("`perm` is deprecated, use `smof permute`")

    # parse arguments, the default function that will ultimately be run is
    # selected according to tthe user-chosen subcommand
    args = parser.parser.parse_args(argv)

    return args


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
            if isinstance(output, FastaEntry):
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
            help="Remove all text from header that follows the first word (delimited by the value of the --delimiter argument, '[ |]' by default)",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "--delimiter",
            help="The regex delimiter used by --reduce-header",
            default="[ |]",
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
            _err("Please provide sequence type (--type)")

        if args.tolower and args.toupper:
            _err("Err, you want me to convert to lower AND upper?")

    def generator(self, args, gen):
        # Catches illegal combinations of arguments
        self._process_args(args)

        args_dict = vars(args)

        return clean(
            gen,
            seq_type=args.type,
            toupper=args.toupper,
            tolower=args.tolower,
            toseq=args.toseq,
            standardize=args.standardize,
            reduce_header=args.reduce_header,
            delimiter=args.delimiter,
            mask_lowercase=args.mask_lowercase,
            mask_irregular=args.mask_irregular,
        )

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
                _err(
                    'The argument for --composition must be three space separated values, e.g. "GC > .5"'
                )
            legal_signs = ("<", "<=", ">=", ">", "==", "=", "!=")
            if not sign in legal_signs:
                _err(
                    "Middle term must be a comparison symbol ('<', '<=', '>=', '>', '==', '=', '!=')"
                )
            if sign == "=":
                sign = "=="
            try:
                per = float(per)
            except ValueError:
                _err("Third value must be a float")
            if not 0 <= per <= 1:
                _err("Third value must be between 0 and 1")
            ch = set(str(ch))

            def evaluate(s):
                c = Counter(s)
                p = sum([c[x] for x in ch]) / len(s)
                return eval("p {} {}".format(sign, per))

            tests.append(evaluate)

        for seq in gen:
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
        return md5sum(
            gen,
            ignore_case=args.ignore_case,
            each_sequence=args.each_sequence,
            all_sequences=args.all_sequences,
            all_headers=args.all_headers,
            replace_header=args.replace_header,
        )


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
        return permute(
            gen,
            seed=args.seed,
            word_size=args.word_size,
            start_offset=args.start_offset,
            end_offset=args.end_offset,
        )


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
        self.force_color = args.force_color
        return reverse(gen, complement=args.complement, no_validate=args.no_validate)


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
        NotImplemented

    def write(self, args, gen, out=sys.stdout):
        """
        This function basically just formats and prints the information in a
        FastaDescription object
        """
        seqsum = sniff(gen)

        out.write(str(seqsum))


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

    def generator(self, args, gen):
        args = self._process_args(args)
        if args.byseq:
            g = stat_seq(
                gen,
                length=args.length,
                counts=args.counts,
                proportion=args.proportion,
                case_sensitive=args.case_sensitive,
                count_lower=args.count_lower,
            )

            for item in g:
                yield args.delimiter.join([str(i) for i in item])
        else:
            # Do I need to count the characters? (much faster if I don't)
            need_count = any(
                (
                    args.counts,
                    args.proportion,
                    args.count_lower,
                    args.type,
                    args.aa_profile,
                )
            )

            g = stat_file(gen, count_characters=need_count)

            if need_count:
                yield g.get_count(
                    count_lower=args.count_lower,
                    case_sensitive=args.case_sensitive,
                    type=args.type,
                    counts=args.counts,
                    proportion=args.proportion,
                )

            if args.length:
                yield g.get_length()

            if args.hist:
                if args.log_hist:
                    yield g.get_hist(title="Flat histogram")
                else:
                    yield g.get_hist()

            if args.log_hist:
                if args.hist:
                    yield g.get_hist(title="Log2 histogram", log=True)
                else:
                    yield g.get_hist(log=True)

            if args.aa_profile:
                if args.hist or args.log_hist:
                    yield g.get_aaprofile(title="AA profile")
                else:
                    yield g.get_aaprofile()


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
        for s in gen:
            yield s

    def write(self, args, gen, out=None):
        p = args.prefix
        N = args.number
        used = set()
        for i, seq in enumerate(self.generator(args, gen)):
            fnum = i // N if args.seqs else i % N
            outfile = "%s%s.fasta" % (p, str(fnum))
            if not outfile in used and os.path.isfile(outfile):
                _err('Split refuses to overwrite "%s"' % outfile)
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
            "-a",
            "--annotate",
            help="Append the subsequence interval to the defline",
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
        self.force_color = args.force_color
        if args.gff:
            sgen = gff_subseq(gen, gff_file=args.gff, color=args.color)
        else:
            sgen = subseq(
                gen,
                a=args.bounds[0],
                b=args.bounds[1],
                color=args.color,
                annotate=args.annotate,
            )

        return sgen


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
        return translate(
            gen, all_frames=args.all_frames, from_start=args.from_start, cds=args.cds
        )


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
        seqs = [s for s in gen]
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
        seqs = [s for s in gen]

        if args.numeric_sort and not args.regex:
            _err("--numeric does nothing unless with --regex")

        # Set type of order determining variable
        if args.numeric_sort:

            def typer(x):
                try:
                    return float(x)
                except ValueError:
                    _err("'{}' cannot be numerically sorted".format(x))

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
                    _err("No match for regex '{}'".format(args.regex))
                except IndexError:
                    _err("Nothing was captured in regex '{}'".format(args.regex))

        elif args.key:
            try:
                key = int(args.key) - 1

                def sortterm(x):
                    try:
                        return typer(x.header.split(args.field_separator)[key])
                    except IndexError:
                        _err("Cannot sort by column '{}', two few columns".format(key))

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
                            _err(str(e))
                    try:
                        return typer(d[key])
                    except KeyError:
                        _err("Could not find key '{}'".format(key))
                    except Exception as e:
                        _err(str(e))

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
        result = consensus(gen, table=args.table)

        if args.table:
            out.write("\t".join(result[0]))
            out.write("\n")
            transpose = [[s.seq[i] for s in seqs] for i in range(0, imax)]
            for row in transpose:
                out.write("\t".join([str(x) for x in row]))
                out.write("\n")
        else:
            result.print()


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
                        _err(
                            "'{}-{}' does not specify a valid range".format(
                                str(s), str(t)
                            )
                        )
                except ValueError:
                    _err("Cannot parse '{}'".format(args.fields))

        return cut(gen, indices=indices, complement=args.complement)


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
        allbut = False
        if args.entries:
            if args.nseqs:
                _err("Please don't use nseqs with --entries")
            try:
                if args.entries[0] == "-":
                    allbut = True
                    nseqs = int(args.entries[1:])
                else:
                    nseqs = int(args.entries)
            except AttributeError:
                _err("-n (--entries) must be a number")
        elif args.nseqs:
            try:
                nseqs = int(re.match(r"-(\d+)", args.nseqs).group(1))
            except AttributeError:
                _err("N must be formatted as '-12'")
        else:
            nseqs = 1

        return head(gen, nseqs=nseqs, first=args.first, last=args.last, allbut=allbut)


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

    def generator(self, args, gen):
        searcher = GrepSearch(args)
        self.force_color = args.force_color
        return searcher.search(gen)


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
            "--first-header",
            help="remove entries with duplicate headers (keep only the first)",
            action="store_true",
            default=False,
        )
        parser.add_argument(
            "--removed",
            metavar="FILE",
            help="With -f, store removed sequences in FILE",
            type=argparse.FileType("w"),
        )
        parser.set_defaults(func=self.func)

    def generator(self, args, gen):
        if args.first_header:
            return uniq_headers(gen, removed=args.removed)
        elif args.pack:
            return pack(gen, sep=args.pack_sep)
        elif args.unpack:
            return unpack(gen, sep=args.pack_sep)
        else:
            return uniq(gen, repeated=args.repeated, uniq=args.uniq, count=args.count)


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
        for seq in gen:
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
            _err("ERROR: do not use -n along with positional number argument")

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
                _err("N must be formatted as '[+-]12'")

        if fromtop:
            i = 1
            for seq in gen:
                if i >= nstring:
                    yield _headtailtrunk(seq, args.first, args.last)
                i += 1
        else:
            from collections import deque

            try:
                lastseqs = deque(maxlen=nstring)
            except ValueError:
                _err("--nseqs argument must be positive")
            for seq in gen:
                lastseqs.append(seq)

            for s in lastseqs:
                yield _headtailtrunk(s, args.first, args.last)


# =======
# EXECUTE
# =======


def main():
    if os.name == "nt":
        os.system("color")  # allows ANSI color on windows
    else:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    args = parse()

    # This is only relevant to Head and Tail
    if "nseqs" in args:
        # This resolve cases where there is a positional filename and no
        # nseqs given.
        # If the first positional argument is a readable filename, treat
        # it as input. Otherwise, try to interpret it as a number
        if os.access(args.nseqs, os.R_OK):
            args.fh = [args.nseqs] + args.fh
            args.nseqs = None

    if "pattern" in args:        
        # If the pattern is readable, it is probably meant to be an input, not
        # a pattern
        if args.file or args.fastain:
            args.fh = [args.pattern] + args.fh
            args.pattern = None

    # If no input is given,
    # and if smof is not reading user input from stdin,
    # assume piped input is from STDIN
    try:
        if not args.fh:
            files = [sys.stdin]
        else:
            files = args.fh
    # If args does not have a .fh argument, then try treating args itself
    # as the input
    except AttributeError:
        files = [args]

    handle_color = ("preserve_color" in args) and bool(args.preserve_color)

    gen = _stream_entries(files, handle_color=handle_color)

    args.func(args, gen, out=sys.stdout)
