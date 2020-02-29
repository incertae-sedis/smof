#!/usr/bin/env python3

import smof.main as smof
import unittest
import argparse
import sys
import tempfile
import math
import os
from collections import Counter
from io import StringIO


class dummy:
    def __init__(self, fh):
        self.fh = [fh]


def get_output(seq, argv):
    argv = [str(s) for s in argv]
    out = StringIO()
    args = smof.parse(argv)
    args.fh = [seq]
    gen = smof.FSeqGenerator(args)
    args.func(args, gen, out=out)
    return out.getvalue().strip().split("\n")


class TestParseHeader(unittest.TestCase):
    def test_firstword(self):
        self.assertEqual(smof.ParseHeader.firstword("abc xyz"), "abc")

    def test_description_present(self):
        self.assertEqual(smof.ParseHeader.description("abc xyz"), "xyz")

    def test_description_absent(self):
        self.assertEqual(smof.ParseHeader.description("abc"), "")

    def test_add_suffix_with_desc(self):
        self.assertEqual(smof.ParseHeader.add_suffix("abc", "s"), "abc|s")

    def test_add_suffix_without_desc(self):
        self.assertEqual(smof.ParseHeader.add_suffix("abc xyz", "s"), "abc|s xyz")

    def test_subseq_with_desc(self):
        self.assertEqual(smof.ParseHeader.subseq("abc", 1, 10), "abc|subseq(1..10)")

    def test_subseq_without_desc(self):
        self.assertEqual(
            smof.ParseHeader.subseq("abc xyz", 1, 10), "abc|subseq(1..10) xyz"
        )

    def test_permute_with_desc(self):
        self.assertEqual(
            smof.ParseHeader.permute("abc", 1, 2, 3),
            "abc|permutation:start=1;end=2;word_size=3",
        )

    def test_permute_without_desc(self):
        self.assertEqual(
            smof.ParseHeader.permute("abc xyz", 1, 2, 3),
            "abc|permutation:start=1;end=2;word_size=3 xyz",
        )


class TestFSeq(unittest.TestCase):
    def test_subseq(self):
        header = "seq1 description"
        seq = "AcGGNttt"
        seqobj = smof.FSeq(header, seq)
        sq = seqobj.subseq(1, 5)
        self.assertEqual(sq.seq, "cGGN")
        self.assertEqual(sq.header, "seq1|subseq(2..5) description")

    def test_getrevcomp_fromStringInput(self):
        seq = "ACGTT"
        self.assertEqual(smof.FSeq.getrevcomp(seq), "AACGT")

    def test_getrevcomp_fromFSeqInput(self):
        header = "seq1"
        seq = "ACGTT"
        seqobj = smof.FSeq(header, seq)
        rc = smof.FSeq.getrevcomp(seqobj)
        self.assertEqual(rc.seq, "AACGT")
        self.assertEqual(rc.header, "seq1|revcom")

    def test_getrevcomp_extended_alphabet(self):
        # test uracil and unknown
        f = "ACGTUNacgtun"
        r = "naacgtNAACGT"
        self.assertEqual(smof.FSeq.getrevcomp(f), r)

        # W = [AT] <--> S = [GC]
        # M = [AC] <--> K = [GT]
        # R = [AG] <--> Y = [CT]
        f = "wmrWMRskySKY"
        r = "RMWrmwYKSyks"
        self.assertEqual(smof.FSeq.getrevcomp(f), r)

        # B = [GTC] <--> V = [ACG]
        # D = [AGT] <--> H = [ACT]
        f = "BDVHbdvh"
        r = "dbhvDBHV"
        self.assertEqual(smof.FSeq.getrevcomp(f), r)

    def test_ungap(self):
        header = "seq1"
        seq = "A.C-G_T"
        seqobj = smof.FSeq(header, seq)
        seqobj.ungap()
        self.assertEqual(seqobj.seq, "ACGT")
        self.assertEqual(seqobj.header, "seq1|ungapped")

    def test_reverse(self):
        header = "seq1"
        seq = "ACGTT"
        seqobj = smof.FSeq(header, seq)
        seqobj.reverse()
        self.assertEqual(seqobj.seq, "TTGCA")
        self.assertEqual(seqobj.header, "seq1|reverse")


class TestStatFun(unittest.TestCase):
    def test_N50(self):
        self.assertEqual(smof.StatFun.N50([1, 2, 3.1]), 3.1)

    def test_mean(self):
        self.assertEqual(smof.StatFun.mean([1, 2, 3]), 2)

    def test_median(self):
        self.assertEqual(smof.StatFun.median([1, 2, 3]), 2)

    def test_sd(self):
        self.assertEqual(smof.StatFun.sd([1, 2, 3]), 1)
        self.assertTrue(math.isnan(smof.StatFun.sd([1])))

    def test_quantile(self):
        self.assertEqual(smof.StatFun.quantile([1, 2, 3], 0.5), 2)

    def test_summary_values(self):
        x = list(range(-5, 11))
        o = smof.StatFun.summary(x)
        self.assertEqual(o["min"], -5)
        self.assertEqual(o["1st_qu"], -1.25)
        self.assertEqual(o["median"], 2.5)
        self.assertEqual(o["3rd_qu"], 6.25)
        self.assertEqual(o["max"], 10)


class TestFileDescription(unittest.TestCase):
    def setUp(self):
        self.seqs = {
            "p-normal": "SMIF",
            "p-selenocysteine": "SMUF",
            "p-unknown": "SMXF",
            "p-ambiguous": "SMBF",
            "p-illegal": "SMOF",
            "p-terminal-stop": "SMIF*",
            "p-internal-stop": "SM*F",
            "p-initial-Met": "MSMIF",
            "p-lowercase": "smif",
            "p-mixedcase": "sMiF",
            "0000": "GGTtaagGCCGGT",
            "0001": "GGTgGCCGGT",
            "0010": "GGTtaaGCCGGT",
            "0011": "GGTGCCGGT",
            "0100": "GGTtaagGCCTAA",
            "0101": "GGTgGCCTAA",
            "0110": "GGTtaaGCCTAA",
            "0111": "GGTGCCTAA",
            "1000": "ATGtaagGCCGGT",
            "1001": "ATGgGCCGGT",
            "1010": "ATGtaaGCCGGT",
            "1011": "ATGGCCGGT",
            "1100": "ATGtaagGCCTAA",
            "1101": "ATGgGCCTAA",
            "1110": "ATGtaaGCCTAA",
            "1111": "ATGGCCTAA",
            "n-gapped": "-GATACA_CAT.TAG.",
            "p-gapped": ".YELL-AT_FAT-PHIL.",
            "illegal": "CHOCOLATE",
        }

    def _prep_fd(self, keys):
        fd = smof.FileDescription()
        for key in keys:
            seq = smof.FSeq(header=key, seq=self.seqs[key])
            fd.add_seq(seq)
        return fd

    def _equal_counts(self, test, true):
        return (
            test.ntype == true.ntype
            and test.ncase == true.ncase
            and test.pfeat == true.pfeat
            and test.nfeat == true.nfeat
            and test.ufeat == true.ufeat
        )

    def test_prot_normal(self):
        fd_test = self._prep_fd(["p-normal"])

        fd_true = smof.FileDescription()
        fd_true.ntype["prot"] += 1
        fd_true.ncase["uppercase"] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_selenocysteine(self):
        fd_test = self._prep_fd(["p-selenocysteine"])

        fd_true = smof.FileDescription()
        fd_true.ntype["prot"] += 1
        fd_true.ncase["uppercase"] += 1
        fd_true.pfeat["selenocysteine"] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_unknown(self):
        fd_test = self._prep_fd(["p-unknown"])

        fd_true = smof.FileDescription()
        fd_true.ntype["prot"] += 1
        fd_true.ncase["uppercase"] += 1
        fd_true.ufeat["unknown"] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_ambiguous(self):
        fd_test = self._prep_fd(["p-ambiguous"])

        fd_true = smof.FileDescription()
        fd_true.ntype["prot"] += 1
        fd_true.ncase["uppercase"] += 1
        fd_true.ufeat["ambiguous"] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_illegal(self):
        fd_test = self._prep_fd(["p-illegal"])

        fd_true = smof.FileDescription()
        fd_true.ntype["illegal"] += 1
        fd_true.ncase["uppercase"] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_terminal_stop(self):
        fd_test = self._prep_fd(["p-terminal-stop"])

        fd_true = smof.FileDescription()
        fd_true.ntype["prot"] += 1
        fd_true.ncase["uppercase"] += 1
        fd_true.pfeat["terminal-stop"] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_internal_stop(self):
        fd_test = self._prep_fd(["p-internal-stop"])

        fd_true = smof.FileDescription()
        fd_true.ntype["prot"] += 1
        fd_true.ncase["uppercase"] += 1
        fd_true.pfeat["internal-stop"] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_initial_met(self):
        fd_test = self._prep_fd(["p-initial-Met"])

        fd_true = smof.FileDescription()
        fd_true.ntype["prot"] += 1
        fd_true.ncase["uppercase"] += 1
        fd_true.pfeat["initial-Met"] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_lowercase(self):
        fd_test = self._prep_fd(["p-lowercase"])

        fd_true = smof.FileDescription()
        fd_true.ntype["prot"] += 1
        fd_true.ncase["lowercase"] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_mixedcase(self):
        fd_test = self._prep_fd(["p-mixedcase"])

        fd_true = smof.FileDescription()
        fd_true.ntype["prot"] += 1
        fd_true.ncase["mixedcase"] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_gapped(self):
        fd_test = self._prep_fd(["p-gapped"])

        fd_true = smof.FileDescription()
        fd_true.ntype["prot"] += 1
        fd_true.ncase["uppercase"] += 1
        fd_true.ufeat["gapped"] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_nucl_gapped(self):
        fd_test = self._prep_fd(["n-gapped"])

        fd_true = smof.FileDescription()
        fd_true.ntype["dna"] += 1
        fd_true.ncase["uppercase"] += 1
        fd_true.ufeat["gapped"] += 1
        fd_true.nfeat["0111"] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_illegal(self):
        fd_test = self._prep_fd(["illegal"])

        fd_true = smof.FileDescription()
        fd_true.ntype["illegal"] += 1
        fd_true.ncase["uppercase"] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_has_start(self):
        start = lambda s: smof.FileDescription._has_start(s)
        self.assertTrue(start("ATGGGT"))
        self.assertFalse(start("GGTATG"))
        self.assertFalse(start("GATGG"))
        self.assertFalse(start("GGGGG"))

    def test_has_stop(self):
        stop = lambda s: smof.FileDescription._has_stop(s)
        for codon in ("TAA", "TAG", "TGA"):
            self.assertTrue(stop("GGG%s" % codon))
            self.assertTrue(stop("GG%s" % codon))
            self.assertFalse(stop("%sGG" % codon))
            self.assertFalse(stop("G%sGG" % codon))
        self.assertTrue(stop("TAATGA"))
        self.assertTrue(stop("TAAgTGA"))

    def test_is_sense(self):
        sense = lambda s: smof.FileDescription._is_sense(s)
        self.assertTrue(sense("ATGGGGCCCTAA"))
        self.assertTrue(sense("CCCGGGCCCAAA"))
        self.assertTrue(sense("CCCGGGCCCTAA"))
        self.assertTrue(sense("CCCGTAAAAAAGG"))
        self.assertFalse(sense("CCCTAACCCAAA"))
        self.assertFalse(sense("ATGTAACCCAAA"))
        self.assertFalse(sense("ATGCCCTAATAA"))

    def test_is_triple(self):
        triple = lambda s: smof.FileDescription._is_triple(s)
        self.assertTrue(triple("ATG"))
        self.assertTrue(triple("ATGGGG"))
        self.assertFalse(triple("ATGGGGg"))
        self.assertFalse(triple("ATGGGGgg"))

    def test_nfeat(self):
        from itertools import product as product

        def test_profile(prof):
            fd = self._prep_fd([prof])
            return bool(fd.nfeat[prof])

        for prof in ["".join(c) for c in product("01", repeat=4)]:
            self.assertTrue(test_profile(prof), "nfeat profile: %s, failed" % prof)


class TestUtilities(unittest.TestCase):
    def setUp(self):
        self.seq = smof.FSeq("seq", "ACDEFSTVWY")

    def test_counter_caser(self):
        self.assertEqual(smof.counter_caser(Counter("Aaa")), {"A": 3})
        self.assertEqual(smof.counter_caser(Counter("Aaa"), True), {"a": 3})

    def test_sum_lower(self):
        self.assertEqual(smof.sum_lower(Counter("AaaFf")), 3)
        self.assertEqual(smof.sum_lower(Counter("AAAFF")), 0)
        self.assertEqual(smof.sum_lower(Counter("AaaF.{")), 2)

    def test_guess_type_input(self):
        # String input
        self.assertEqual(smof.guess_type("FFFF"), "prot")
        # Counter object input
        self.assertEqual(smof.guess_type(Counter("FFFF")), "prot")
        # FSeq object input
        self.assertEqual(smof.guess_type(smof.FSeq("s1", "FFFF")), "prot")
        # Gaps should be ignored
        self.assertEqual(smof.guess_type(smof.FSeq("s1", "F-F_F.F")), "prot")
        # Case should be ignored
        self.assertEqual(smof.guess_type(smof.FSeq("s1", "ffff")), "prot")

    def test_guess_type_dna(self):
        self.assertEqual(smof.guess_type("GATACA"), "dna")
        self.assertEqual(smof.guess_type("GATACANNN"), "dna")
        self.assertEqual(smof.guess_type("NNNNNN"), "dna")

    def test_guess_type_rna(self):
        self.assertEqual(smof.guess_type("GAUACA"), "rna")

    def test_guess_type_prot(self):
        self.assertEqual(smof.guess_type("FAMNX"), "prot")
        self.assertEqual(smof.guess_type("XXXXX"), "prot")

    def test_guess_type_illegal(self):
        self.assertEqual(smof.guess_type("DAMO"), "illegal")
        self.assertEqual(smof.guess_type("DAM!"), "illegal")
        # A nucleotide sequence can't have both U and T
        self.assertEqual(smof.guess_type("GATU"), "illegal")
        # Space is illegal
        self.assertEqual(smof.guess_type("DAM "), "illegal")
        # Gaps should NOT be counted as illegal
        self.assertNotEqual(smof.guess_type("D.A-M_"), "illegal")
        # * should not be illegal (indicates STOP in protein sequence)
        self.assertNotEqual(smof.guess_type("DA*M*"), "illegal")

    def test_guess_type_ambiguous(self):
        self.assertEqual(smof.guess_type("A"), "ambiguous")
        self.assertEqual(smof.guess_type("AT"), "ambiguous")
        self.assertNotEqual(smof.guess_type("ATG"), "ambiguous")
        self.assertNotEqual(smof.guess_type("AUG"), "ambiguous")
        # Greater than 80% nucleotide characters with ambiguous is dna
        self.assertEqual(smof.guess_type("ATGGR"), "ambiguous")
        self.assertEqual(smof.guess_type("ATGGGR"), "dna")
        # Sequences containing only ambiguous nucleotides (could be dna or
        # protein) are counted as ambiguous regardless of lenght
        self.assertEqual(smof.guess_type("WYS"), "ambiguous")
        self.assertEqual(smof.guess_type("RYSWKMDBHV"), "ambiguous")
        # But if one unambiguous aa is added ('F')
        self.assertEqual(smof.guess_type("FRYSWKMDBHV"), "prot")

    def test_counting_number(self):
        import argparse

        self.assertRaises(argparse.ArgumentTypeError, smof.counting_number, 0)
        self.assertRaises(argparse.ArgumentTypeError, smof.counting_number, -1)

    def test_positive_int(self):
        import argparse

        self.assertRaises(argparse.ArgumentTypeError, smof.positive_int, -1)
        try:
            smof.positive_int(0)
            zero_calls_exception = False
        except argparse.ArgumentTypeError:
            zero_calls_exception = True
        self.assertFalse(zero_calls_exception)

    def test_headtailtrunk_first(self):
        # Note: argparse will only allow positive integers
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=1, last=0).seq, "A")
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=5, last=0).seq, "ACDEF")
        self.assertEqual(
            smof.headtailtrunk(seq=self.seq, first=20, last=0).seq, "ACDEFSTVWY"
        )

    def test_headtailtrunk_last(self):
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=0, last=1).seq, "Y")
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=0, last=5).seq, "STVWY")
        self.assertEqual(
            smof.headtailtrunk(seq=self.seq, first=0, last=20).seq, "ACDEFSTVWY"
        )

    def test_headtailtrunk_firstandlast(self):
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=1, last=1).seq, "A...Y")
        self.assertEqual(
            smof.headtailtrunk(seq=self.seq, first=2, last=3).seq, "AC...VWY"
        )
        self.assertEqual(
            smof.headtailtrunk(seq=self.seq, first=5, last=5).seq, "ACDEFSTVWY"
        )
        self.assertEqual(
            smof.headtailtrunk(seq=self.seq, first=6, last=6).seq, "ACDEFSTVWY"
        )

    def test_headtailtrunk_doublezero(self):
        self.assertRaises(SystemExit, smof.headtailtrunk, seq=self.seq, first=0, last=0)

    def test_headtailtrunk_doublenone(self):
        self.assertEqual(smof.headtailtrunk(seq=self.seq).seq, "ACDEFSTVWY")


class TestAmbiguous2Perl(unittest.TestCase):
    def test_simple(self):
        self.assertEqual(smof.ambiguous2perl("Y"), "[CT]")

    def test_escaped(self):
        self.assertEqual(smof.ambiguous2perl("\Y"), "Y")
        self.assertEqual(smof.ambiguous2perl("\YY"), "Y[CT]")
        self.assertEqual(smof.ambiguous2perl("Y\Y+Y\\{2"), "[CT]Y+[CT]\\{2")

    def test_bracketed(self):
        self.assertEqual(smof.ambiguous2perl("[Y]"), "[CT]")
        self.assertEqual(smof.ambiguous2perl("[Y\]]"), "[CT\]]")
        # Of course, '[' and ']' should NEVER be in a DNA sequence, but ...
        self.assertEqual(smof.ambiguous2perl("[\[Y\]]"), "[\[CT\]]")


class TestFSeqGenerator(unittest.TestCase):
    def setUp(self):
        self.seq1 = smof.FSeq(header="seq1", seq="ACGTA")
        self.seq2 = smof.FSeq(header="seq2", seq="GGTT")
        self.seq1_spaced = smof.FSeq(header="seq1", seq="AC GTA")
        self.seq2_spaced = smof.FSeq(header="seq2", seq="GGTT")
        self.seq1_weird = smof.FSeq(header="seq1 >weirdness", seq="ACGTA")
        self.seq1_funky = smof.FSeq(
            header="seq1|asdf:!@(*#& !@#$%^&*())_+", seq="ACGTA"
        )

        self.good = [">seq1\n", "ACGT\n", "A\n", ">seq2\n", "GGT\n", "T\n"]
        self.good_empty_lines = [">seq1", "ACGT", "A", "\n", ">seq2", "GGT", "T", "\n"]
        self.weird_empty_lines = [
            "\n",
            ">seq1",
            "ACGT",
            "\n",
            "A",
            "\n",
            ">seq2",
            "GGT",
            "T",
            "\n",
        ]
        self.spaced = [" >seq1", "AC GT", "A", " >seq2 ", " GGT", "T "]
        self.well_commented = [
            "# this is a comment",
            "# so is this",
            ">seq1",
            "ACGT",
            "A",
            ">seq2",
            "GGT",
            "T",
        ]
        self.interspersed_comments = [
            "# this is a comment",
            ">seq1",
            "ACGT",
            "A",
            "# so is this",
            ">seq2",
            "GGT",
            "T",
        ]
        self.bad_first = ["A", ">seq1", "ACGT", "A", ">seq2", "GGT", "T"]
        self.empty_seq = [">seq1", ">seq2", "GGT", "T"]
        self.empty_last_seq = [">seq1", "ACGT", "A", ">seq2"]
        self.internal_gt = [">seq1 >weirdness", "ACGT", "A"]
        self.funky_header = [">seq1|asdf:!@(*#& !@#$%^&*())_+", "ACGT", "A"]
        self.no_sequence = []

    def cmp_seqs(self, fh, exp_seqs):
        seq = dummy(fh)
        g = smof.FSeqGenerator(seq)
        obs_seqs = [s for s in g.next()]
        for obs, exp in zip(obs_seqs, exp_seqs):
            if (obs.header != exp.header) or (obs.seq != exp.seq):
                print([obs.header, exp.header])
                return False
        return True

    def is_valid(self, fh):
        seq = dummy(fh)
        try:
            g = smof.FSeqGenerator(seq)
            out = [s for s in g.next()]
            return True
        except BaseException:
            return False

    def test_good(self):
        self.assertTrue(self.cmp_seqs(self.good, (self.seq1, self.seq2)))

    def test_good_empty_lines(self):
        self.assertTrue(self.cmp_seqs(self.good_empty_lines, (self.seq1, self.seq2)))

    def test_weird_empty_lines(self):
        self.assertTrue(self.cmp_seqs(self.weird_empty_lines, (self.seq1, self.seq2)))

    def test_spaced(self):
        self.assertTrue(
            self.cmp_seqs(self.spaced, (self.seq1_spaced, self.seq2_spaced))
        )

    def test_well_commented(self):
        self.assertTrue(self.cmp_seqs(self.well_commented, (self.seq1, self.seq2)))

    def test_interspersed_comments(self):
        self.assertTrue(
            self.cmp_seqs(self.interspersed_comments, (self.seq1, self.seq2))
        )

    def test_funky_header(self):
        self.assertTrue(self.cmp_seqs(self.funky_header, [self.seq1_funky]))

    def test_internal_gt(self):
        self.assertTrue(self.cmp_seqs(self.internal_gt, [self.seq1_weird]))

    def test_bad_first(self):
        self.assertFalse(self.is_valid(self.bad_first))

    def test_empty_seq(self):
        # empty sequences are now supported
        self.assertTrue(self.is_valid(self.empty_seq))

    def test_empty_last_seq(self):
        # empty sequences are now supported
        self.assertTrue(self.is_valid(self.empty_last_seq))

    def test_no_sequence(self):
        self.assertTrue(self.is_valid(self.no_sequence))


class TestMd5sum(unittest.TestCase):
    def setUp(self):
        self.seqs = [">asdf", "ASDF", ">qwer", "TYUI"]

    def test_default(self):
        self.assertEqual(
            get_output(self.seqs, ["md5sum"]), ["28fd532b933aaa89d2188b98241a8b46"]
        )

    def test_eachseq(self):
        self.assertEqual(
            get_output(self.seqs, ["md5sum", "-q"]),
            [
                "asdf\t6d87a19f011653459575ceb722db3b69",
                "qwer\t6e9758614cca89162b2d19922de103bb",
            ],
        )

    def test_replaceseq(self):
        self.assertEqual(
            get_output(self.seqs, ["md5sum", "-r"]),
            [
                ">6d87a19f011653459575ceb722db3b69",
                "ASDF",
                ">6e9758614cca89162b2d19922de103bb",
                "TYUI",
            ],
        )

    def test_headers(self):
        self.assertEqual(
            get_output(self.seqs, ["md5sum", "-d"]),
            ["c69874b898abb180ac71bd99bc16f8fb"],
        )

    def test_seqs(self):
        self.assertEqual(
            get_output(self.seqs, ["md5sum", "-s"]),
            ["ed9b124094bc93e7f611da252d06f628"],
        )


class TestClean(unittest.TestCase):
    def setUp(self):
        self.seq = [">a", " gAtA cA-NY "]
        self.aaseq = [">p", " gAtA cA-NB "]
        self.longseq = [">l", "A" * 91]
        self.header = [">l a", "A", ">m|a", "A"]
        self.gapamb = [">a", "yy--_.ATT"]
        self.nonstandard_dna = [">a", ".-_XxCAT"]
        self.nonstandard_pro = [">a", ".-_GANDALF"]

    def test_default(self):
        self.assertEqual(get_output(self.seq, ["clean"])[1], "gAtAcA-NY")

    def test_case(self):
        self.assertEqual(get_output(self.seq, ["clean", "-u"])[1], "GATACA-NY")
        self.assertEqual(get_output(self.seq, ["clean", "-l"])[1], "gataca-ny")

    def test_masking(self):
        self.assertEqual(
            get_output(self.seq, ["clean", "-t", "nucl", "-m"])[1], "NANANA-NY"
        )
        self.assertEqual(
            get_output(self.seq, ["clean", "-t", "nucl", "-mr"])[1], "NANANA-NN"
        )

    def test_type(self):
        for d in ["n", "nu", "nuc", "nucl", "dna"]:
            self.assertEqual(
                get_output(self.seq, ["clean", "-t", d, "-r"])[1], "gAtAcA-NN"
            )
        for d in ["p", "pro", "prot", "protein", "aa", "pep"]:
            self.assertEqual(
                get_output(self.aaseq, ["clean", "-t", d, "-r"])[1], "gAtAcA-NX"
            )

    def test_toseq(self):
        self.assertEqual(get_output([">a", "ASD!@(#*& D"], ["clean", "-x"])[1], "ASDD")

    def test_irregulars(self):
        self.assertEqual(
            get_output([">p", "YbJuZ"], ["clean", "-t", "p", "-r"])[1], "YXXXX"
        )
        self.assertEqual(
            get_output([">n", "ATRySWkMDbHVG"], ["clean", "-t", "n", "-r"])[1],
            "ATNNNNNNNNNNG",
        )

        # Unambiguously illegal characters are not masked
        self.assertEqual(
            get_output([">p", "YOU]"], ["clean", "-t", "p", "-r"])[1], "YOX]"
        )
        self.assertEqual(
            get_output([">n", "ATryjG*"], ["clean", "-t", "n", "-r"])[1], "ATNNjG*"
        )

    def test_nonstandard_dna(self):
      self.assertEqual(
          get_output(self.nonstandard_dna, ["clean", "-d", "-t", "n"])[1], "---NnCAT"
      )

    def test_nonstandard_pro(self):
      self.assertEqual(
          get_output(self.nonstandard_pro, ["clean", "-d", "-t", "p"])[1], "---GANDALF"
      )

    def test_wrap(self):
        self.assertEqual(get_output(self.longseq, ["clean", "-w", "30"])[1], "A" * 30)
        self.assertEqual(get_output(self.longseq, ["clean", "-w", "30"])[4], "A")
        # test no-wrap
        self.assertEqual(get_output(self.longseq, ["clean", "-w", "0"])[1], "A" * 91)

    def test_reduce_header(self):
        self.assertEqual(get_output(self.header, ["clean", "-s"])[0], ">l")
        self.assertEqual(get_output(self.header, ["clean", "-s"])[2], ">m")

    def test_reduce_and_mask(self):
        self.assertEqual(get_output(self.gapamb, ["clean", "-uxrt", "n"])[1], "NNATT")


class TestFilter(unittest.TestCase):
    def setUp(self):
        self.seq = [">a", "ASDFX", ">b", "ASDF", ">c", "ASD", ">d", ""]

    def test_shorter_than(self):
        self.assertEqual(
            get_output(self.seq, ["filter", "-s", 5])[0::2], [">a", ">b", ">c", ">d"]
        )
        self.assertEqual(
            get_output(self.seq, ["filter", "-s", 4])[0::2], [">b", ">c", ">d"]
        )
        self.assertEqual(get_output(self.seq, ["filter", "-s", 3])[0::2], [">c", ">d"])
        self.assertEqual(get_output(self.seq, ["filter", "-s", 0])[0::2], [">d"])

    def test_longer_than(self):
        self.assertEqual(
            get_output(self.seq, ["filter", "-l", 0])[0::2], [">a", ">b", ">c", ">d"]
        )
        self.assertEqual(
            get_output(self.seq, ["filter", "-l", 3])[0::2], [">a", ">b", ">c"]
        )
        self.assertEqual(get_output(self.seq, ["filter", "-l", 4])[0::2], [">a", ">b"])
        self.assertEqual(get_output(self.seq, ["filter", "-l", 5])[0::2], [">a"])

    def test_composition(self):
        comp = [">a", "AAAAG.....", ">b", "AG........", ">c", "AAX......."]
        self.assertEqual(
            get_output(comp, ["filter", "-c", "X == 0"])[0::2], [">a", ">b"]
        )
        self.assertEqual(
            get_output(comp, ["filter", "-c", "X = 0"])[0::2], [">a", ">b"]
        )
        self.assertEqual(get_output(comp, ["filter", "-c", "X != 0"])[0::2], [">c"])
        self.assertEqual(get_output(comp, ["filter", "-c", "AG > .3"])[0::2], [">a"])
        self.assertEqual(
            get_output(comp, ["filter", "-c", "AG < .3"])[0::2], [">b", ">c"]
        )
        self.assertEqual(get_output(comp, ["filter", "-c", "AG >= .5"])[0::2], [">a"])
        self.assertEqual(
            get_output(comp, ["filter", "-c", "AG <= .5"])[0::2], [">a", ">b", ">c"]
        )


class TestHeaderGrep(unittest.TestCase):
    def setUp(self):
        self.headers = [
            ">gg sco 12",
            "A",
            ">gg bob 48a",
            "A",
            ">gl har 61",
            "A",
            ">aL har 61",
            "A",
        ]

    def test_default(self):
        self.assertEqual(
            get_output(self.headers, ["grep", "-y", "bob"]), [">gg bob 48a", "A"]
        )
        self.assertEqual(
            get_output(self.headers, ["grep", "-y", "gg"]),
            [">gg sco 12", "A", ">gg bob 48a", "A"],
        )

    def test_perl(self):
        self.assertEqual(
            get_output(self.headers, ["grep", "-yP", ".g"]),
            [">gg sco 12", "A", ">gg bob 48a", "A"],
        )
        self.assertEqual(
            get_output(self.headers, ["grep", "-yP", "bob|sco"]),
            [">gg sco 12", "A", ">gg bob 48a", "A"],
        )
        self.assertEqual(
            get_output(self.headers, ["grep", "-yP", "\d+[a-z]"]), [">gg bob 48a", "A"]
        )

    def test_invert(self):
        self.assertEqual(
            get_output(self.headers, ["grep", "-yvP", "^g"]), [">aL har 61", "A"]
        )

    def test_case_sensitive(self):
        self.assertEqual(
            get_output(self.headers, ["grep", "-yI", "aL"]), [">aL har 61", "A"]
        )
        self.assertEqual(get_output(self.headers, ["grep", "-yI", "al"]), [""])
        self.assertEqual(
            get_output(self.headers, ["grep", "-y", "al"]), [">aL har 61", "A"]
        )
        self.assertEqual(
            get_output(self.headers, ["grep", "-y", "aL"]), [">aL har 61", "A"]
        )

    def test_count(self):
        self.assertEqual(get_output(self.headers, ["grep", "-cP", "gg"]), ["2"])

    def test_only_matching(self):
        self.assertEqual(
            get_output([">a;glob.", "GACFADE"], ["grep", "-oP", "g..b\."]), ["glob."]
        )

    def test_only_matching_wrap(self):
        self.assertEqual(
            get_output(
                [">a;glob.", "GACFADE"], ["grep", "-w", "a;([^.]+)", "-o", "glob"]
            ),
            ["glob"],
        )

    def test_line_regexp(self):
        self.assertEqual(get_output(self.headers, ["grep", "-x", "gg"]), [""])
        self.assertEqual(
            get_output(self.headers, ["grep", "-x", "gg sco 12"]), [">gg sco 12", "A"]
        )

    def test_exact(self):
        self.assertEqual(get_output(self.headers, ["grep", "-X", "gg"]), [""])
        self.assertEqual(
            get_output(self.headers, ["grep", "-X", "gg sco 12"]), [">gg sco 12", "A"]
        )


class TestSequenceGrep(unittest.TestCase):
    def setUp(self):
        self.seqs = [
            ">a",
            "AAGATACA",
            ">b",
            "GAACATAACAT",
            ">c",
            "aaaaa",
            ">d",
            "aaaaaa",
            ">e",
            "A",
        ]
        self.revseqs = [">a", "AAG", ">b", "CTT"]

    def test_default(self):
        self.assertEqual(
            get_output(self.seqs, ["grep", "-qy", "gataca"]), [">a", "AAGATACA"]
        )

    def test_ambiguous_nucl_encodings(self):
        for h, q in [
            ("M", "AC"),
            ("R", "AG"),
            ("W", "AT"),
            ("S", "CG"),
            ("Y", "CT"),
            ("K", "GT"),
            ("V", "ACG"),
            ("H", "ACT"),
            ("D", "AGT"),
            ("B", "CGT"),
            ("N", "ACGT"),
        ]:
            self.assertNotEqual(
                get_output([">{}".format(h), q], ["grep", "-qyG", "^{}+$".format(h)]),
                [""],
            )
            compl = "".join(set("ACGT") - set(q))
            if compl:
                self.assertEqual(
                    get_output([">{}".format(h), compl], ["grep", "-qyG", h]), [""]
                )

    def test_ambiguous_nucl_regex(self):
        self.assertEqual(
            get_output(self.seqs, ["grep", "-qyG", "R{4}Y"]), [">a", "AAGATACA"]
        )
        self.assertEqual(
            get_output(self.seqs, ["grep", "-qyG", "[^Y]{4}Y"]), [">a", "AAGATACA"]
        )

    def test_count(self):
        self.assertEqual(get_output(self.seqs, ["grep", "-cq", "aa"]), ["4"])

    def test_matches(self):
        self.assertEqual(get_output(self.seqs, ["grep", "-qm", "aa"]), ["8"])

    def test_count_matches(self):
        self.assertEqual(get_output(self.seqs, ["grep", "-qcm", "aa"]), ["4\t8"])

    def test_both_strands(self):
        self.assertEqual(
            get_output(self.revseqs, ["grep", "-qy", "AA"]), self.revseqs[0:2]
        )
        self.assertEqual(get_output(self.revseqs, ["grep", "-qby", "AA"]), self.revseqs)
        self.assertEqual(
            get_output(self.revseqs, ["grep", "-qy", "AG"]), self.revseqs[0:2]
        )
        self.assertEqual(get_output(self.revseqs, ["grep", "-qby", "AG"]), self.revseqs)

    def test_gff(self):
        self.assertEqual(
            get_output(self.seqs, ["grep", "--gff", "CAT"]),
            [
                "b\tsmof-{}\tregex_match\t4\t6\t.\t.\t.\t.".format(smof.__version__),
                "b\tsmof-{}\tregex_match\t9\t11\t.\t.\t.\t.".format(smof.__version__),
            ],
        )

    def test_gff_context(self):
        self.assertEqual(
            get_output(self.seqs, ["grep", "--gff", "-A 1", "CAT"]),
            [
                "b\tsmof-{}\tregex_match\t4\t7\t.\t.\t.\t.".format(smof.__version__),
                "b\tsmof-{}\tregex_match\t9\t11\t.\t.\t.\t.".format(smof.__version__),
            ],
        )
        self.assertEqual(
            get_output(self.seqs, ["grep", "--gff", "-B 1", "CAT"]),
            [
                "b\tsmof-{}\tregex_match\t3\t6\t.\t.\t.\t.".format(smof.__version__),
                "b\tsmof-{}\tregex_match\t8\t11\t.\t.\t.\t.".format(smof.__version__),
            ],
        )

    def test_gff_seqid(self):
        self.assertEqual(
            get_output([">a|b.1 desc", "AC"], ["grep", "--gff", "A"]),
            ["a|b.1\tsmof-{}\tregex_match\t1\t1\t.\t.\t.\t.".format(smof.__version__)],
        )

    def test_only_matching(self):
        self.assertEqual(
            get_output([">a", "GACFADE"], ["grep", "-qoP", "A."])[1::2], ["AC", "AD"]
        )
        self.assertEqual(
            get_output([">a", "GAACFADE"], ["grep", "-qoP", "A.*?D"])[1::2], ["AACFAD"]
        )

    def test_only_matching_context(self):
        self.assertEqual(
            get_output([">a", "GACFADE"], ["grep", "-qoP", "-A 1", "A."])[1::2],
            ["ACF", "ADE"],
        )
        self.assertEqual(
            get_output([">a", "GACFADE"], ["grep", "-qoP", "-A 2", "A."])[1::2],
            ["ACFA", "ADE"],
        )
        self.assertEqual(
            get_output([">a", "GACFADE"], ["grep", "-qoP", "-B 1", "A."])[1::2],
            ["GAC", "FAD"],
        )
        self.assertEqual(
            get_output([">a", "GACFADE"], ["grep", "-qoP", "-B 2", "A."])[1::2],
            ["GAC", "CFAD"],
        )

    def test_only_matching_context_both(self):
        self.assertEqual(
            get_output([">a", "GAAGGGTTA"], ["grep", "-qoPb", "-A 1", "AA"])[1::2],
            ["AAG", "GTT"],
        )

    def test_only_matching_context_reverse(self):
        self.assertEqual(
            get_output([">a", "GAAGGGTTA"], ["grep", "-qoPr", "-A 1", "AA"])[1], "GTT"
        )

    def test_only_matching_wrap(self):
        self.assertEqual(
            get_output([">a", "GACFADE"], ["grep", "-qw", "CF(..)", "-o", "AD"])[1],
            "AD",
        )

    def test_only_matching_wrap_reverse(self):
        self.assertEqual(
            get_output([">a", "GACFADE"], ["grep", "-qw", "CF(..)", "-o", "AD"])[1],
            "AD",
        )
        self.assertEqual(
            get_output([">a", "GAAGGGTTA"], ["grep", "-qbw", "AAC(..)", "-o", "CC"])[1],
            "GG",
        )

    def test_gapped_search(self):
        self.assertEqual(
            get_output([">a", "GA-CFADE"], ["grep", "-qgy", "AC"])[1], "GA-CFADE"
        )

    def test_gapped_search_only(self):
        self.assertEqual(
            get_output([">a", "GA-CFADE"], ["grep", "-qgyo", "AC"])[1], "A-C"
        )
        self.assertEqual(
            get_output([">a", "--GA--C-F-ADE"], ["grep", "-qgyo", "ACF"])[1], "A--C-F"
        )
        self.assertEqual(
            get_output([">a", "G--ACF-ADE"], ["grep", "-qgyo", "ACF"])[1], "ACF"
        )

    def test_gapped_search_only_revcom(self):
        self.assertEqual(
            get_output([">a", "GATA-CA"], ["grep", "-qyorg", "GTA"])[1], "TA-C"
        )
        self.assertEqual(
            get_output([">a", "GATA-CA"], ["grep", "-qyobg", "GTA"])[1], "TA-C"
        )

    def test_line_regexp(self):
        self.assertEqual(get_output(self.seqs, ["grep", "-qx", "GAA"]), [""])
        self.assertEqual(
            get_output(self.seqs, ["grep", "-qx", "GAACATAACAT"]), [">b", "GAACATAACAT"]
        )
        self.assertEqual(
            get_output(self.seqs, ["grep", "-Pqx", "GAA.*"]), [">b", "GAACATAACAT"]
        )

    def test_exact(self):
        # Partial exact matches return nothing
        self.assertEqual(get_output(self.seqs, ["grep", "-qX", "GAA"]), [""])
        # Full exact matches return everything
        self.assertEqual(
            get_output(self.seqs, ["grep", "-qX", "GAACATAACAT"]), [">b", "GAACATAACAT"]
        )

    def test_fastain(self):
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b">a\nGAT")
        f.close()
        self.assertEqual(
            get_output(self.seqs, ["grep", "-y", "--fastain", f.name]),
            [">a", "AAGATACA"],
        )
        self.assertEqual(
            get_output(self.seqs, ["grep", "-yo", "--fastain", f.name])[1], "GAT"
        )
        self.assertEqual(
            get_output(self.seqs, ["grep", "--gff", "--fastain", f.name])[0].split(
                "\t"
            )[3:5],
            ["3", "5"],
        )
        os.unlink(f.name)


class TestGrepBadCombinations(unittest.TestCase):
    def setUp(self):
        self.seq = [">a", "A"]

    def test_wrap_incompatible_options(self):
        self.assertRaises(
            SystemExit, get_output, self.seq, ["grep", "-Pw", "a(b)", "a"]
        )
        self.assertRaises(
            SystemExit, get_output, self.seq, ["grep", "-Gw", "a(b)", "a"]
        )

    def test_mv(self):
        self.assertRaises(SystemExit, get_output, self.seq, ["grep", "-mv", "a"])

    def test_seq_only_matches_against_header(self):
        self.assertRaises(SystemExit, get_output, self.seq, ["grep", "-b", "a"])
        self.assertRaises(SystemExit, get_output, self.seq, ["grep", "-r", "a"])
        self.assertRaises(SystemExit, get_output, self.seq, ["grep", "-G", "a"])

    def test_contradictory_options(self):
        self.assertRaises(SystemExit, get_output, self.seq, ["grep", "-yY", "a"])

    def test_only_matching_incompatible_options(self):
        self.assertRaises(SystemExit, get_output, self.seq, ["grep", "-ov", "a"])

    def test_exact_incompatible_options(self):
        self.assertRaises(SystemExit, get_output, self.seq, ["grep", "-PX", "a"])
        self.assertRaises(SystemExit, get_output, self.seq, ["grep", "-GX", "a"])
        self.assertRaises(
            SystemExit, get_output, self.seq, ["grep", "--gff", "-x", "a"]
        )
        self.assertRaises(
            SystemExit, get_output, self.seq, ["grep", "--gff", "-X", "a"]
        )


class TestHeadandTail(unittest.TestCase):
    def setUp(self):
        self.seq = [">a", "GATACA", ">b", "GALLIF", ">c", "SPARTA", ">d", "ATHENS"]

    def test_defaults(self):
        self.assertEqual(get_output(self.seq, ["head"]), [">a", "GATACA"])
        self.assertEqual(
            get_output(self.seq, ["head", "-2"]), [">a", "GATACA", ">b", "GALLIF"]
        )

    def test_head_n(self):
        self.assertEqual(get_output(self.seq, ["head", "-n", "1"]), [">a", "GATACA"])
        self.assertEqual(get_output(self.seq, ["head", "-n", "+1"]), [">a", "GATACA"])
        self.assertEqual(get_output(self.seq, ["head", "-n", "-1"]), self.seq[0:-2])

    def test_tail_n(self):
        self.assertEqual(get_output(self.seq, ["tail", "-n", "1"]), [">d", "ATHENS"])
        self.assertEqual(get_output(self.seq, ["tail", "-n", "-1"]), [">d", "ATHENS"])
        self.assertEqual(get_output(self.seq, ["tail", "-n", "+2"]), self.seq[2:])

    def test_tail(self):
        self.assertEqual(get_output(self.seq, ["tail"]), [">d", "ATHENS"])
        self.assertEqual(
            get_output(self.seq, ["tail", "-2"]), [">c", "SPARTA", ">d", "ATHENS"]
        )
        self.assertEqual(
            get_output(self.seq, ["tail", "+2"]),
            [">b", "GALLIF", ">c", "SPARTA", ">d", "ATHENS"],
        )

    def test_fl(self):
        self.assertEqual(get_output(self.seq, ["head", "-f", 2, "-l", 1])[1], "GA...A")
        self.assertEqual(get_output(self.seq, ["tail", "-f", 2, "-l", 1])[1], "AT...S")


class TestPermute(unittest.TestCase):
    def setUp(self):
        self.seq = [">a", "WHEREISMYTARDIS"]

    def test_default(self):
        self.assertEqual(
            get_output(self.seq, ["permute", "--seed", 42]),
            [">a|permutation:start=0;end=0;word_size=1", "MTISSADYHEIERWR"],
        )

    def test_word_size(self):
        self.assertEqual(
            get_output(self.seq, ["permute", "--seed", 42, "-w", 3]),
            [">a|permutation:start=0;end=0;word_size=3", "TARREISMYDISWHE"],
        )
        self.assertEqual(
            get_output(self.seq, ["permute", "--seed", 42, "-w", 5]),
            [">a|permutation:start=0;end=0;word_size=5", "ARDISISMYTWHERE"],
        )

    def test_offsets(self):
        self.assertEqual(
            get_output(self.seq, ["permute", "--seed", 42, "-w", 4, "-s", 3]),
            [">a|permutation:start=3;end=0;word_size=4", "WHERDISMYTAREIS"],
        )
        self.assertEqual(
            get_output(self.seq, ["permute", "--seed", 123, "-w", 4, "-s", 5]),
            [">a|permutation:start=5;end=0;word_size=4", "WHEREISTARDISMY"],
        )
        self.assertEqual(
            get_output(self.seq, ["permute", "--seed", 123, "-w", 4, "-e", 3]),
            [">a|permutation:start=0;end=3;word_size=4", "YTAREISMWHERDIS"],
        )


class TestReverse(unittest.TestCase):
    def setUp(self):
        self.seq = [">a1", "LIVED", ">b2", "MILLER"]
        self.reverse = [">a1|reverse", "DEVIL", ">b2|reverse", "RELLIM"]
        # This test sequence is adapted from the Sequence Manipulation Suite
        # (http://www.bioinformatics.org/sms2/rev_comp.html)
        self.seq2 = [
            ">s1 sample1",
            "garkbdctymvhu",
            ">s2 sample2",
            "ctymvhgarkbda",
            ">s3 sample3",
            "ccccccccccga",
        ]
        self.seq2_revcomp = [
            ">s1|revcom sample1",
            "adbkraghvmytc",
            ">s2|revcom sample2",
            "thvmytcdbkrag",
            ">s3|revcom sample3",
            "tcgggggggggg",
        ]

    def test_default(self):
        self.assertEqual(get_output(self.seq, ["reverse"]), self.reverse)

    def test_reverse_complement(self):
        self.assertEqual(get_output(self.seq2, ["reverse", "-cV"]), self.seq2_revcomp)


class TestSample(unittest.TestCase):
    def setUp(self):
        self.seqs = [">1", "A", ">2", "A", ">3", "A", ">4", "A", ">5", "A"]

    def test_default(self):
        self.assertEqual(get_output(self.seqs, ["sample", "--seed", "5"]), [">5", "A"])
        self.assertEqual(
            get_output(self.seqs, ["sample", "--seed", "5", "-n", "2"]),
            [">5", "A", ">3", "A"],
        )
        self.assertEqual(
            get_output(self.seqs, ["sample", "-n", "2", "--seed", "123"]),
            [">1", "A", ">3", "A"],
        )


class TestSort(unittest.TestCase):
    def setUp(self):
        self.unsorted = [
            ">g=c;d=100",
            "AAA",
            ">g=d;d=30",
            "AA",
            ">g=b;d=9",
            "AAAA",
            ">g=a;d=200",
            "A",
        ]
        self.default = [
            ">g=a;d=200",
            "A",
            ">g=b;d=9",
            "AAAA",
            ">g=c;d=100",
            "AAA",
            ">g=d;d=30",
            "AA",
        ]
        self.default_reverse = [
            ">g=d;d=30",
            "AA",
            ">g=c;d=100",
            "AAA",
            ">g=b;d=9",
            "AAAA",
            ">g=a;d=200",
            "A",
        ]
        self.length = [
            ">g=a;d=200",
            "A",
            ">g=d;d=30",
            "AA",
            ">g=c;d=100",
            "AAA",
            ">g=b;d=9",
            "AAAA",
        ]
        self.regex = [
            ">g=c;d=100",
            "AAA",
            ">g=a;d=200",
            "A",
            ">g=d;d=30",
            "AA",
            ">g=b;d=9",
            "AAAA",
        ]
        self.regex_numeric = [
            ">g=b;d=9",
            "AAAA",
            ">g=d;d=30",
            "AA",
            ">g=c;d=100",
            "AAA",
            ">g=a;d=200",
            "A",
        ]

    def test_default(self):
        self.assertEqual(get_output(self.unsorted, ["sort"]), self.default)

    def test_default(self):
        self.assertEqual(
            get_output(self.unsorted, ["sort", "-r"]), self.default_reverse
        )

    def test_length_sort(self):
        self.assertEqual(get_output(self.unsorted, ["sort", "-l"]), self.length)

    def test_regex_sort(self):
        self.assertEqual(
            get_output(self.unsorted, ["sort", "-x", "d=(\d+)"]), self.regex
        )

    def test_numeric_sort(self):
        self.assertEqual(
            get_output(self.unsorted, ["sort", "-x", "d=(\d+)", "-n"]),
            self.regex_numeric,
        )


class TestSubseq(unittest.TestCase):
    def setUp(self):
        self.seq = [">a", "GATACA"]
        self.aaseq = [">p", "PICKLE"]

    def test_default(self):
        self.assertEqual(get_output(self.seq, ["subseq", "-b", 1, 1])[1], "G")
        self.assertEqual(get_output(self.seq, ["subseq", "-b", 5, 6])[1], "CA")

    def test_overbounds(self):
        self.assertEqual(get_output(self.seq, ["subseq", "-b", 1, 100])[1], "GATACA")
        self.assertRaises(SystemExit, get_output, self.seq, ["subseq", "-b", 7, 7])

    def test_revcomp(self):
        self.assertEqual(get_output(self.seq, ["subseq", "-b", 3, 6])[1], "TACA")
        self.assertEqual(get_output(self.aaseq, ["subseq", "-b", 1, 3])[1], "PIC")
        # guess_type function adientifies GATACA as dna, then takes revcomp
        self.assertEqual(get_output(self.seq, ["subseq", "-b", 6, 3])[1], "TGTA")
        # PICKLE however doesn't appear to be dna, so reversing does nothing
        self.assertEqual(get_output(self.aaseq, ["subseq", "-b", 3, 1])[1], "PIC")


class TestTranslateDNA(unittest.TestCase):
    def test_find_max_orf(self):
        self.assertEqual(smof.find_max_orf("", from_start=True), (None, 0))
        self.assertEqual(smof.find_max_orf("T", from_start=True), (None, 0))
        self.assertEqual(smof.find_max_orf("ATG", from_start=True), (0, 3))
        self.assertEqual(smof.find_max_orf("TAAATG", from_start=True), (3, 3))
        self.assertEqual(smof.find_max_orf("TAAATGTAG", from_start=True), (3, 3))
        self.assertEqual(smof.find_max_orf("TAAATGATGTAG", from_start=True), (3, 6))
        self.assertEqual(smof.find_max_orf("taaatgatgtag", from_start=True), (3, 6))
        self.assertEqual(smof.find_max_orf("AATG", from_start=True), (1, 3))
        self.assertEqual(smof.find_max_orf("AAAAATGATGTTTTAA", from_start=True), (4, 9))
        self.assertEqual(smof.find_max_orf("AAAAATGATGTTTTAA", from_start=False), (0, 15))
        self.assertEqual(smof.find_max_orf("aaaaatgatgttttaa", from_start=False), (0, 15))
        self.assertEqual(smof.find_max_orf("ATG", from_start=False), (0, 3))
        self.assertEqual(smof.find_max_orf("TAAATG", from_start=False), (1, 3))
        self.assertEqual(smof.find_max_orf("TAAATGTAG", from_start=False), (1, 6))
        self.assertEqual(smof.find_max_orf("TAAATGATGTAG", from_start=False), (2, 9))

    def test_simple(self):
        self.assertEqual(smof.translate_dna(""), "")
        self.assertEqual(smof.translate_dna("T"), "")
        self.assertEqual(smof.translate_dna("ATG"), "M")
        self.assertEqual(smof.translate_dna("atG"), "M")
        self.assertEqual(smof.translate_dna("ATGT"), "M")
        self.assertEqual(smof.translate_dna("AT"), "")
        self.assertEqual(smof.translate_dna(""), "")
        self.assertEqual(smof.translate_dna("TAA"), "*")
        self.assertEqual(smof.translate_dna("TAAATG"), "*M")
        self.assertEqual(smof.translate_dna("taAatG"), "*M")
        self.assertEqual(
            smof.translate_dna(
                "TTTTCTTATTGTTTCTCCTACTGCTTATCATAATGATTGTCGTAGTGGCTTCCTCATCGTCTCCCCCACCGCCTACCACAACGACTGCCGCAGCGGATTACTAATAGTATCACCAACAGCATAACAAAAAGAATGACGAAGAGGGTTGCTGATGGTGTCGCCGACGGCGTAGCAGAAGGAGTGGCGGAGGGG"
            ),
            "FSYCFSYCLS**LS*WLPHRLPHRLPQRLPQRITNSITNSITKRMTKRVADGVADGVAEGVAEG",
        )

    def test_ambiguous_handling(self):
        self.assertEqual(smof.translate_dna("TTTATGYNNATG"), "FMXM")
        self.assertEqual(smof.translate_dna("YTT"), "X")
        self.assertEqual(smof.translate_dna("-AT---GT_..TT-"), "MF")

    def test_all_frames(self):
        self.assertEqual(smof.get_orf("TTTATGT", all_frames=True), "FM")
        self.assertEqual(smof.get_orf("ATGTGA", all_frames=False), "M*")
        self.assertEqual(smof.get_orf("aTgtGa", all_frames=False), "M*")
        self.assertEqual(smof.get_orf("ATGTGA", all_frames=True), "M")
        self.assertEqual(
            smof.get_orf("ATGTGAT", all_frames=True), "CD"
        )  # match 2nd frame
        self.assertEqual(
            smof.get_orf("ATGTGATTAATTTATGTTTTTTTTT", all_frames=True), "VINLCFF"
        )  # match in 3rd frame

    def test_require_start_and_all_frames(self):
        self.assertEqual(
            smof.get_orf(
                "ATGTGATTAATTTATGTTTTTTTTT", all_frames=True, from_start=True
            ),
            "MFFF",
        )  # match in 3rd frame

    def test_get_cds(self):
        self.assertEqual( smof.get_orf(""), "")
        self.assertEqual( smof.get_orf("T"), "")
        self.assertEqual( smof.get_orf("TT"), "")
        self.assertEqual( smof.get_orf("TTTTTT", from_start=True), "")

        self.assertEqual(
            smof.get_orf(
                "ATG", all_frames=True, from_start=True, translate=False
            ),
            "ATG",
        )
        self.assertEqual(
            smof.get_orf(
                "TATGTTTTGA", all_frames=True, from_start=True, translate=False
            ),
            "ATGTTT",
        )
        self.assertEqual(
            smof.get_orf(
                "TATGTTTTGA", all_frames=True, from_start=True, translate=False
            ),
            "ATGTTT",
        )


class TestUniq(unittest.TestCase):
    def setUp(self):
        self.all_uniq = [">a", "CAT", ">b", "HAT", ">c", "A"]
        self.unsorted = [">b", "HAT", ">a", "CAT", ">c", "A"]
        self.repeated = [">a", "CAT", ">b", "HAT", ">c", "A", ">b", "HAT", ">c", "A"]
        self.redundant = [">a", "CAT", ">b", "CAT", ">c", "CAT", ">d", "HAT"]
        self.packed = [">a|b|c", "CAT", ">d", "HAT"]
        self.repeated_header = [">a", "CAT", ">d", "HAT", ">a", "MATH"]

    def test_default(self):
        self.assertEqual(get_output(self.all_uniq, ["uniq"]), self.all_uniq)
        self.assertEqual(get_output(self.repeated, ["uniq"]), self.all_uniq)
        self.assertEqual(get_output(self.unsorted, ["uniq"]), self.unsorted)

    def test_uniq(self):
        self.assertEqual(get_output(self.all_uniq, ["uniq", "-u"]), self.all_uniq)
        self.assertEqual(get_output(self.repeated, ["uniq", "-u"]), [">a", "CAT"])

    def test_duplicated(self):
        self.assertEqual(get_output(self.all_uniq, ["uniq", "-d"]), [""])
        self.assertEqual(
            get_output(self.repeated, ["uniq", "-d"]), [">b", "HAT", ">c", "A"]
        )

    def test_duplicated(self):
        self.assertEqual(
            get_output(self.repeated_header, ["uniq", "-f"]),
            [">a", "MATH", ">d", "HAT"],
        )

    def test_count(self):
        self.assertEqual(
            get_output(self.all_uniq, ["uniq", "-c"]), ["1\ta", "1\tb", "1\tc"]
        )
        self.assertEqual(
            get_output(self.repeated, ["uniq", "-c"]), ["1\ta", "2\tb", "2\tc"]
        )

    def test_pack(self):
        # pack does not guarantee conservation of order, so I compare sets
        self.assertEqual(
            set(get_output(self.redundant, ["uniq", "-p", "-z", "|"])),
            set([">a|b|c", "CAT", ">d", "HAT"]),
        )

    def test_unpack(self):
        self.assertEqual(
            get_output(self.packed, ["uniq", "-P", "-z", "|"]), self.redundant
        )


class TestWc(unittest.TestCase):
    def setUp(self):
        self.seq = [">a", "CAT", ">b", "HAT", ">c", "A"]

    def test_default(self):
        self.assertEqual(get_output(self.seq, ["wc"]), ["3\t7"])

    def test_nseqs(self):
        self.assertEqual(get_output(self.seq, ["wc", "-l"]), ["3"])

    def test_nchars(self):
        self.assertEqual(get_output(self.seq, ["wc", "-m"]), ["7"])


if __name__ == "__main__":
    unittest.main()
