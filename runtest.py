#!/usr/bin/env python3

import smof
import unittest
import argparse
import sys
from tempfile import TemporaryFile
from collections import Counter
from io import StringIO


class TestFSeq(unittest.TestCase):
    def test_subseq(self):
        header = 'seq1'
        seq = 'AcGGNttt'
        seqobj = smof.FSeq(header, seq)
        sq = seqobj.subseq(1, 5)
        self.assertEqual(sq.seq, 'cGGN')
        self.assertEqual(sq.header, 'seq1|SUBSEQ(2..5)')

    def test_getrevcomp_fromStringInput(self):
        seq = 'ACGTT'
        self.assertEqual(smof.FSeq.getrevcomp(seq), 'AACGT')

    def test_getrevcomp_fromFSeqInput(self):
        header = 'seq1'
        seq = 'ACGTT'
        seqobj = smof.FSeq(header,seq)
        rc = smof.FSeq.getrevcomp(seqobj)
        self.assertEqual(rc.seq, 'AACGT')
        self.assertEqual(rc.header, 'seq1|REVCOM')

    def test_ungap(self):
        header = 'seq1'
        seq = 'A.C-G_T'
        seqobj = smof.FSeq(header,seq)
        seqobj.ungap()
        self.assertEqual(seqobj.seq, 'ACGT')
        self.assertEqual(seqobj.header, 'seq1|UNGAPPED')

    def test_reverse(self):
        header = 'seq1'
        seq = 'ACGTT'
        seqobj = smof.FSeq(header,seq)
        seqobj.reverse()
        self.assertEqual(seqobj.seq, 'TTGCA')
        self.assertEqual(seqobj.header, 'seq1|REVERSE')

class TestStatFun(unittest.TestCase):
    def test_N50(self):
        self.assertEqual(smof.StatFun.N50([1,2,3.1]), 3.1)

    def test_mean(self):
        self.assertEqual(smof.StatFun.mean([1,2,3]), 2)

    def test_median(self):
        self.assertEqual(smof.StatFun.median([1,2,3]), 2)

    def test_sd(self):
        self.assertEqual(smof.StatFun.sd([1,2,3]), 1)

    def test_quantile(self):
        self.assertEqual(smof.StatFun.quantile([1,2,3], 0.5), 2)

    def test_summary_values(self):
        x = list(range(-5, 11))
        o = smof.StatFun.summary(x)
        self.assertEqual(o['min'], -5)
        self.assertEqual(o['1st_qu'], -1.25)
        self.assertEqual(o['median'], 2.5)
        self.assertEqual(o['3rd_qu'], 6.25)
        self.assertEqual(o['max'], 10)

class TestFileDescription(unittest.TestCase):
    def setUp(self):
        self.seqs = {
            'p-normal':'SMIF',
            'p-selenocysteine':'SMUF',
            'p-unknown':'SMXF',
            'p-ambiguous':'SMBF',
            'p-illegal':'SMOF',
            'p-terminal-stop':'SMIF*',
            'p-internal-stop':'SM*F',
            'p-initial-Met':'MSMIF',
            'p-lowercase':'smif',
            'p-mixedcase':'sMiF',
            '0000':'GGTtaagGCCGGT',
            '0001':'GGTgGCCGGT',
            '0010':'GGTtaaGCCGGT',
            '0011':'GGTGCCGGT',
            '0100':'GGTtaagGCCTAA',
            '0101':'GGTgGCCTAA',
            '0110':'GGTtaaGCCTAA',
            '0111':'GGTGCCTAA',
            '1000':'ATGtaagGCCGGT',
            '1001':'ATGgGCCGGT',
            '1010':'ATGtaaGCCGGT',
            '1011':'ATGGCCGGT',
            '1100':'ATGtaagGCCTAA',
            '1101':'ATGgGCCTAA',
            '1110':'ATGtaaGCCTAA',
            '1111':'ATGGCCTAA',
            'n-gapped':'-GATACA_CAT.TAG.',
            'p-gapped':'.YELL-AT_FAT-PHIL.',
            'illegal':'CHOCOLATE'
        }

    def _prep_fd(self, keys):
        fd = smof.FileDescription()
        for key in keys:
            seq = smof.FSeq(header=key, seq=self.seqs[key])
            fd.add_seq(seq)
        return(fd)

    def _equal_counts(self, test, true):
        return(
            test.ntype == true.ntype and
            test.ncase == true.ncase and
            test.pfeat == true.pfeat and
            test.nfeat == true.nfeat and
            test.ufeat == true.ufeat
        )

    def test_prot_normal(self):
        fd_test = self._prep_fd(['p-normal'])

        fd_true = smof.FileDescription()
        fd_true.ntype['prot'] += 1
        fd_true.ncase['uppercase'] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_selenocysteine(self):
        fd_test = self._prep_fd(['p-selenocysteine'])

        fd_true = smof.FileDescription()
        fd_true.ntype['prot'] += 1
        fd_true.ncase['uppercase'] += 1
        fd_true.pfeat['selenocysteine'] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_unknown(self):
        fd_test = self._prep_fd(['p-unknown'])

        fd_true = smof.FileDescription()
        fd_true.ntype['prot'] += 1
        fd_true.ncase['uppercase'] += 1
        fd_true.ufeat['unknown'] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_ambiguous(self):
        fd_test = self._prep_fd(['p-ambiguous'])

        fd_true = smof.FileDescription()
        fd_true.ntype['prot'] += 1
        fd_true.ncase['uppercase'] += 1
        fd_true.ufeat['ambiguous'] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_illegal(self):
        fd_test = self._prep_fd(['p-illegal'])

        fd_true = smof.FileDescription()
        fd_true.ntype['illegal'] += 1
        fd_true.ncase['uppercase'] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_terminal_stop(self):
        fd_test = self._prep_fd(['p-terminal-stop'])

        fd_true = smof.FileDescription()
        fd_true.ntype['prot'] += 1
        fd_true.ncase['uppercase'] += 1
        fd_true.pfeat['terminal-stop'] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_internal_stop(self):
        fd_test = self._prep_fd(['p-internal-stop'])

        fd_true = smof.FileDescription()
        fd_true.ntype['prot'] += 1
        fd_true.ncase['uppercase'] += 1
        fd_true.pfeat['internal-stop'] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_initial_met(self):
        fd_test = self._prep_fd(['p-initial-Met'])

        fd_true = smof.FileDescription()
        fd_true.ntype['prot'] += 1
        fd_true.ncase['uppercase'] += 1
        fd_true.pfeat['initial-Met'] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_lowercase(self):
        fd_test = self._prep_fd(['p-lowercase'])

        fd_true = smof.FileDescription()
        fd_true.ntype['prot'] += 1
        fd_true.ncase['lowercase'] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_mixedcase(self):
        fd_test = self._prep_fd(['p-mixedcase'])

        fd_true = smof.FileDescription()
        fd_true.ntype['prot'] += 1
        fd_true.ncase['mixedcase'] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_prot_gapped(self):
        fd_test = self._prep_fd(['p-gapped'])

        fd_true = smof.FileDescription()
        fd_true.ntype['prot'] += 1
        fd_true.ncase['uppercase'] += 1
        fd_true.ufeat['gapped'] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_nucl_gapped(self):
        fd_test = self._prep_fd(['n-gapped'])

        fd_true = smof.FileDescription()
        fd_true.ntype['dna'] += 1
        fd_true.ncase['uppercase'] += 1
        fd_true.ufeat['gapped'] += 1
        fd_true.nfeat['0111'] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_illegal(self):
        fd_test = self._prep_fd(['illegal'])

        fd_true = smof.FileDescription()
        fd_true.ntype['illegal'] += 1
        fd_true.ncase['uppercase'] += 1

        self.assertTrue(self._equal_counts(fd_true, fd_test))

    def test_has_start(self):
        start = lambda s: smof.FileDescription._has_start(s)
        self.assertTrue(start('ATGGGT'))
        self.assertFalse(start('GGTATG'))
        self.assertFalse(start('GATGG'))
        self.assertFalse(start('GGGGG'))

    def test_has_stop(self):
        stop = lambda s: smof.FileDescription._has_stop(s)
        for codon in ('TAA', 'TAG', 'TGA'):
            self.assertTrue(stop('GGG%s' % codon))
            self.assertTrue(stop('GG%s' % codon))
            self.assertFalse(stop('%sGG' % codon))
            self.assertFalse(stop('G%sGG' % codon))
        self.assertTrue(stop('TAATGA'))
        self.assertTrue(stop('TAAgTGA'))

    def test_is_sense(self):
        sense = lambda s: smof.FileDescription._is_sense(s)
        self.assertTrue(sense('ATGGGGCCCTAA'))
        self.assertTrue(sense('CCCGGGCCCAAA'))
        self.assertTrue(sense('CCCGGGCCCTAA'))
        self.assertTrue(sense('CCCGTAAAAAAGG'))
        self.assertFalse(sense('CCCTAACCCAAA'))
        self.assertFalse(sense('ATGTAACCCAAA'))
        self.assertFalse(sense('ATGCCCTAATAA'))

    def test_is_triple(self):
        triple = lambda s: smof.FileDescription._is_triple(s)
        self.assertTrue(triple('ATG'))
        self.assertTrue(triple('ATGGGG'))
        self.assertFalse(triple('ATGGGGg'))
        self.assertFalse(triple('ATGGGGgg'))

    def test_nfeat(self):
        from itertools import product as product
        def test_profile(prof):
            fd = self._prep_fd([prof])
            return(bool(fd.nfeat[prof]))
        for prof in [''.join(c) for c in product('01', repeat=4)]:
            self.assertTrue(test_profile(prof), 'nfeat profile: %s, failed' % prof)

class TestUtilities(unittest.TestCase):
    def setUp(self):
        self.seq = smof.FSeq('seq', 'ACDEFSTVWY')

    def test_counter_caser(self):
        self.assertEqual(smof.counter_caser(Counter('Aaa')), {'A':3})
        self.assertEqual(smof.counter_caser(Counter('Aaa'), True), {'a':3})

    def test_sum_lower(self):
        self.assertEqual(smof.sum_lower(Counter('AaaFf')), 3)
        self.assertEqual(smof.sum_lower(Counter('AAAFF')), 0)
        self.assertEqual(smof.sum_lower(Counter('AaaF.{')), 2)

    def test_guess_type_input(self):
        # String input
        self.assertEqual(smof.guess_type('FFFF'), 'prot')
        # Counter object input
        self.assertEqual(smof.guess_type(Counter('FFFF')), 'prot')
        # FSeq object input
        self.assertEqual(smof.guess_type(smof.FSeq('s1', 'FFFF')), 'prot')
        # Gaps should be ignored
        self.assertEqual(smof.guess_type(smof.FSeq('s1', 'F-F_F.F')), 'prot')
        # Case should be ignored
        self.assertEqual(smof.guess_type(smof.FSeq('s1', 'ffff')), 'prot')

    def test_guess_type_dna(self):
        self.assertEqual(smof.guess_type('GATACA'), 'dna')
        self.assertEqual(smof.guess_type('GATACANNN'), 'dna')
        self.assertEqual(smof.guess_type('NNNNNN'), 'dna')

    def test_guess_type_rna(self):
        self.assertEqual(smof.guess_type('GAUACA'), 'rna')

    def test_guess_type_prot(self):
        self.assertEqual(smof.guess_type('FAMNX'), 'prot')
        self.assertEqual(smof.guess_type('XXXXX'), 'prot')

    def test_guess_type_illegal(self):
        self.assertEqual(smof.guess_type('DAMO'), 'illegal')
        self.assertEqual(smof.guess_type('DAM!'), 'illegal')
        # A nucleotide sequence can't have both U and T
        self.assertEqual(smof.guess_type('GATU'), 'illegal')
        # Space is illegal
        self.assertEqual(smof.guess_type('DAM '), 'illegal')
        # Gaps should NOT be counted as illegal
        self.assertNotEqual(smof.guess_type('D.A-M_'), 'illegal')
        # * should not be illegal (indicates STOP in protein sequence)
        self.assertNotEqual(smof.guess_type('DA*M*'), 'illegal')

    def test_guess_type_ambiguous(self):
        self.assertEqual(smof.guess_type('A'), 'ambiguous')
        self.assertEqual(smof.guess_type('AT'), 'ambiguous')
        self.assertNotEqual(smof.guess_type('ATG'), 'ambiguous')
        self.assertNotEqual(smof.guess_type('AUG'), 'ambiguous')
        # Greater than 80% nucleotide characters with ambiguous is dna
        self.assertEqual(smof.guess_type('ATGGR'), 'ambiguous')
        self.assertEqual(smof.guess_type('ATGGGR'), 'dna')
        # Sequences containing only ambiguous nucleotides (could be dna or
        # protein) are counted as ambiguous regardless of lenght
        self.assertEqual(smof.guess_type('WYS'), 'ambiguous')
        self.assertEqual(smof.guess_type('RYSWKMDBHV'), 'ambiguous')
        # But if one unambiguous aa is added ('F')
        self.assertEqual(smof.guess_type('FRYSWKMDBHV'), 'prot')

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
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=1, last=0).seq, 'A')
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=5, last=0).seq, 'ACDEF')
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=20, last=0).seq, 'ACDEFSTVWY')

    def test_headtailtrunk_last(self):
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=0, last=1).seq, 'Y')
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=0, last=5).seq, 'STVWY')
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=0, last=20).seq, 'ACDEFSTVWY')

    def test_headtailtrunk_firstandlast(self):
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=1, last=1).seq, 'A...Y')
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=2, last=3).seq, 'AC...VWY')
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=5, last=5).seq, 'ACDEFSTVWY')
        self.assertEqual(smof.headtailtrunk(seq=self.seq, first=6, last=6).seq, 'ACDEFSTVWY')

    def test_headtailtrunk_doublezero(self):
        self.assertRaises(SystemExit, smof.headtailtrunk, seq=self.seq, first=0, last=0)

    def test_headtailtrunk_doublenone(self):
        self.assertEqual(smof.headtailtrunk(seq=self.seq).seq, 'ACDEFSTVWY')

class TestAmbiguous2Perl(unittest.TestCase):
    def test_simple(self):
        self.assertEqual(smof.ambiguous2perl('Y'), '[CT]')

    def test_escaped(self):
        self.assertEqual(smof.ambiguous2perl('\Y'), 'Y')
        self.assertEqual(smof.ambiguous2perl('\YY'), 'Y[CT]')
        self.assertEqual(smof.ambiguous2perl('Y\Y+Y\\{2'), '[CT]Y+[CT]\\{2')

    def test_bracketed(self):
        self.assertEqual(smof.ambiguous2perl('[Y]'), '[CT]')
        self.assertEqual(smof.ambiguous2perl('[Y\]]'), '[CT\]]')
        # Of course, '[' and ']' should NEVER be in a DNA sequence, but ...
        self.assertEqual(smof.ambiguous2perl('[\[Y\]]'), '[\[CT\]]')

class TestFSeqGenerator(unittest.TestCase):
    def setUp(self):
        self.seq1 = smof.FSeq(header='seq1', seq='ACGTA')
        self.seq2 = smof.FSeq(header='seq2', seq='GGTT')
        self.seq1_spaced = smof.FSeq(header='seq1', seq='AC GTA')
        self.seq2_spaced = smof.FSeq(header='seq2', seq='GGTT')
        self.seq1_weird = smof.FSeq(header="seq1 >weirdness", seq='ACGTA')
        self.seq1_funky = smof.FSeq(header="seq1|asdf:!@(*#& !@#$%^&*())_+", seq="ACGTA")

        self.good = [
            ">seq1\n", "ACGT\n", "A\n",
            ">seq2\n", "GGT\n", "T\n"
        ]
        self.good_empty_lines = [
            ">seq1", "ACGT", "A", "\n",
            ">seq2", "GGT", "T", "\n"
        ]
        self.weird_empty_lines = [
            "\n", ">seq1", "ACGT", "\n", "A", "\n",
            ">seq2", "GGT", "T", "\n"
        ]
        self.spaced = [
            " >seq1", "AC GT", "A",
            " >seq2 ", " GGT", "T "
        ]
        self.well_commented = [
            "# this is a comment",
            "# so is this",
            ">seq1", "ACGT", "A",
            ">seq2", "GGT", "T"
        ]
        self.interspersed_comments = [
            "# this is a comment",
            ">seq1", "ACGT", "A",
            "# so is this",
            ">seq2", "GGT", "T"
        ]
        self.bad_first = [
            "A",
            ">seq1", "ACGT", "A",
            ">seq2", "GGT", "T"
        ]
        self.empty_seq = [
            ">seq1",
            ">seq2", "GGT", "T"
        ]
        self.empty_last_seq = [
            ">seq1", "ACGT", "A",
            ">seq2"
        ]
        self.internal_gt = [
            ">seq1 >weirdness", "ACGT", "A"
        ]
        self.funky_header = [
            ">seq1|asdf:!@(*#& !@#$%^&*())_+", "ACGT", "A"
        ]
        self.no_sequence = []

    def cmp_seqs(self, fh, exp_seqs):
        g = smof.FSeqGenerator(fh)
        obs_seqs = [s for s in g.next()]
        for obs, exp in zip(obs_seqs, exp_seqs):
            if (obs.header != exp.header) or (obs.seq != exp.seq):
                print([obs.header, exp.header])
                return(False)
        return(True)

    def is_valid(self, fh):
        try:
            g = smof.FSeqGenerator(fh)
            out = [s for s in g.next()]
            return(True)
        except BaseException:
            return(False)

    def test_good(self):
        self.assertTrue(self.cmp_seqs(self.good, (self.seq1, self.seq2)))

    def test_good_empty_lines(self):
        self.assertTrue(self.cmp_seqs(self.good_empty_lines, (self.seq1, self.seq2)))

    def test_weird_empty_lines(self):
        self.assertTrue(self.cmp_seqs(self.weird_empty_lines, (self.seq1, self.seq2)))

    def test_spaced(self):
        self.assertTrue(self.cmp_seqs(self.spaced, (self.seq1_spaced, self.seq2_spaced)))

    def test_well_commented(self):
        self.assertTrue(self.cmp_seqs(self.well_commented, (self.seq1, self.seq2)))

    def test_interspersed_comments(self):
        self.assertTrue(self.cmp_seqs(self.interspersed_comments, (self.seq1, self.seq2)))

    def test_funky_header(self):
        self.assertTrue(self.cmp_seqs(self.funky_header, [self.seq1_funky]))

    def test_internal_gt(self):
        self.assertTrue(self.cmp_seqs(self.internal_gt, [self.seq1_weird]))

    def test_bad_first(self):
        self.assertFalse(self.is_valid(self.bad_first))

    def test_empty_seq(self):
        self.assertFalse(self.is_valid(self.empty_seq))

    def test_empty_last_seq(self):
        self.assertFalse(self.is_valid(self.empty_last_seq))

    def test_no_sequence(self):
        self.assertFalse(self.is_valid(self.no_sequence))


def get_output(seq, argv):
    argv = [str(s) for s in argv]
    out = StringIO()
    gen = smof.FSeqGenerator(seq)
    args = smof.parse(argv)
    args.func(args, gen, out=out)
    return(out.getvalue().strip().split("\n"))

class TestChksum(unittest.TestCase):
    def setUp(self):
        self.seqs = [
            '>asdf', 'ASDF',
            '>qwer', 'TYUI'
        ]
    def test_default(self):
        self.assertEqual(
            get_output(self.seqs, ['chksum']),
            ['28fd532b933aaa89d2188b98241a8b46'])
    def test_eachseq(self):
        self.assertEqual(
            get_output(self.seqs, ['chksum', '-q']),
                ['asdf\t6d87a19f011653459575ceb722db3b69',
                 'qwer\t6e9758614cca89162b2d19922de103bb']
            )
    def test_headers(self):
        self.assertEqual(
            get_output(self.seqs, ['chksum', '-d']),
            ['c69874b898abb180ac71bd99bc16f8fb'])
    def test_seqs(self):
        self.assertEqual(
            get_output(self.seqs, ['chksum', '-s']),
            ['ed9b124094bc93e7f611da252d06f628'])

class TestClean(unittest.TestCase):
    def setUp(self):
        self.seq = ['>a', ' gAtA cA-NY ']
        self.aaseq = ['>p', ' gAtA cA-NB ']

    def test_default(self):
        self.assertEqual(get_output(self.seq, ['clean'])[1], 'gAtAcA-NY')

    def test_case(self):
        self.assertEqual(get_output(self.seq, ['clean', '-u'])[1], 'GATACA-NY')
        self.assertEqual(get_output(self.seq, ['clean', '-l'])[1], 'gataca-ny')

    def test_masking(self):
        self.assertEqual(get_output(self.seq, ['clean', '-t', 'nucl', '-m'])[1], 'NANANA-NY')
        self.assertEqual(get_output(self.seq, ['clean', '-t', 'nucl', '-mr'])[1], 'NANANA-NN')

    def test_type(self):
        for d in ['n', 'nu', 'nuc', 'nucl', 'dna']:
            self.assertEqual(get_output(self.seq, ['clean', '-t', d, '-r'])[1], 'gAtAcA-NN')
        for d in ['p', 'pro', 'prot', 'protein', 'aa', 'pep']:
            self.assertEqual(get_output(self.aaseq, ['clean', '-t', d, '-r'])[1], 'gAtAcA-NX')

    def test_toseq(self):
        self.assertEqual(get_output(['>a', 'ASD!@(#*& D'], ['clean', '-x'])[1], 'ASDD')

    def test_irregulars(self):
        self.assertEqual(get_output(['>p', 'YbJuZ'], ['clean', '-t', 'p', '-r'])[1], 'YXXXX')
        self.assertEqual(get_output(['>n', 'ATRySWkMDbHVG'], ['clean', '-t', 'n', '-r'])[1], 'ATNNNNNNNNNNG')

        # Unambiguously illegal characters are not masked
        self.assertEqual(get_output(['>p', 'YOU]'], ['clean', '-t', 'p', '-r'])[1], 'YOX]')
        self.assertEqual(get_output(['>n', 'ATryjG*'], ['clean', '-t', 'n', '-r'])[1], 'ATNNjG*')

class TestFasta2Csv(unittest.TestCase):
    def setUp(self):
        self.seq = [
            '>a1', 'AAA',
            '>b2', 'CCC'
        ]
        self.tsv = [
            'a1\tAAA',
            'b2\tCCC'
        ]
        self.csv = [
            'a1,AAA',
            'b2,CCC'
        ]
        self.headers = [
            'seqid\tseq',
            'a1\tAAA',
            'b2\tCCC'
        ]
    def test_default(self):
        self.assertEqual(get_output(self.seq, ['fasta2csv']), self.tsv)

    def test_delimiter(self):
        self.assertEqual(get_output(self.seq, ['fasta2csv', '-d', ',']), self.csv)

    def test_header(self):
        self.assertEqual(get_output(self.seq, ['fasta2csv', '-r']), self.headers)

class TestHeaderGrep(unittest.TestCase):
    def setUp(self):
        self.headers = [
            '>gg sco 12', 'A',
            '>gg bob 48a', 'A',
            '>gl har 61', 'A',
            '>aL har 61', 'A'
        ]

    def test_default(self):
        self.assertEqual(get_output(self.headers,
            ['grep', '-y', 'bob']),
            ['>gg bob 48a', 'A'])
        self.assertEqual(get_output(self.headers,
            ['grep', '-y', 'gg']),
            ['>gg sco 12', 'A', '>gg bob 48a', 'A'])

    def test_perl(self):
        self.assertEqual(get_output(self.headers,
            ['grep', '-yP', '.g']),
            ['>gg sco 12', 'A', '>gg bob 48a', 'A'])
        self.assertEqual(get_output(self.headers,
            ['grep', '-yP', 'bob|sco']),
            ['>gg sco 12', 'A', '>gg bob 48a', 'A'])
        self.assertEqual(get_output(self.headers,
            ['grep', '-yP', '\d+[a-z]']),
            ['>gg bob 48a', 'A'])

    def test_invert(self):
        self.assertEqual(get_output(self.headers,
            ['grep', '-yvP', '^g']),
            ['>aL har 61', 'A'])

    def test_case_sensitive(self):
        self.assertEqual(get_output(self.headers,
            ['grep', '-yI', 'aL']),
            ['>aL har 61', 'A'])
        self.assertEqual(get_output(self.headers,
            ['grep', '-yI', 'al']),
            [''])
        self.assertEqual(get_output(self.headers,
            ['grep', '-y', 'al']),
            ['>aL har 61', 'A'])
        self.assertEqual(get_output(self.headers,
            ['grep', '-y', 'aL']),
            ['>aL har 61', 'A'])

    def test_count(self):
        self.assertEqual(get_output(self.headers, ['grep', '-cP', 'gg']), ['2'])

class TestSequenceGrep(unittest.TestCase):
    def setUp(self):
        self.seqs = [
            '>a', 'AAGATACA',
            '>b', 'GAACATAACAT',
            '>c', 'aaaaa',
            '>d', 'aaaaaa',
            '>e', 'A'
        ]
        self.revseqs = [
            '>a', 'AAG',
            '>b', 'CTT']

    def test_default(self):
        self.assertEqual(get_output(self.seqs,
            ['grep', '-qy', 'gataca']),
            ['>a', 'AAGATACA'])

    def test_ambiguous_nucl_encodings(self):
            for h,q in [
             ('M', 'AC'),
             ('R', 'AG'),
             ('W', 'AT'),
             ('S', 'CG'),
             ('Y', 'CT'),
             ('K', 'GT'),
             ('V', 'ACG'),
             ('H', 'ACT'),
             ('D', 'AGT'),
             ('B', 'CGT'),
             ('N', 'ACGT')]:
                self.assertNotEqual(get_output(['>{}'.format(h), q], ['grep', '-qyB', '^{}+$'.format(h)]), [''])
                compl = ''.join(set('ACGT') - set(q))
                if compl:
                    self.assertEqual(get_output(['>{}'.format(h), compl], ['grep', '-qyB', h]), [''])

    def test_ambiguous_nucl_regex(self):
        self.assertEqual(get_output(self.seqs, ['grep', '-qyB', 'R{4}Y']), ['>a', 'AAGATACA'])
        self.assertEqual(get_output(self.seqs, ['grep', '-qyB', '[^Y]{4}Y']), ['>a', 'AAGATACA'])

    def test_count(self):
        self.assertEqual(get_output(self.seqs, ['grep', '-cq', 'aa']), ['4'])

    def test_matches(self):
        self.assertEqual(get_output(self.seqs, ['grep', '-qm', 'aa']), ['8'])

    def test_count_matches(self):
        self.assertEqual(get_output(self.seqs, ['grep', '-qcm', 'aa']), ['4\t8'])

    def test_both_strands(self):
        self.assertEqual(get_output(self.revseqs, ['grep', '-qy', 'AA']), self.revseqs[0:2])
        self.assertEqual(get_output(self.revseqs, ['grep', '-qby', 'AA']), self.revseqs)
        self.assertEqual(get_output(self.revseqs, ['grep', '-qy', 'AG']), self.revseqs[0:2])
        self.assertEqual(get_output(self.revseqs, ['grep', '-qby', 'AG']), self.revseqs)

    def test_gff(self):
        self.assertEqual(
            get_output(self.seqs,
            ['grep', '--gff', 'CAT']),
            ['b	smof-1.9.0	regex_match	4	6	.	.	.	.',
            'b	smof-1.9.0	regex_match	9	11	.	.	.	.'])


class TestHeadandTail(unittest.TestCase):
    def setUp(self):
        self.seq=['>a','GATACA',
                  '>b','GALLIF',
                  '>c','SPARTA']

    def test_defaults(self):
        self.assertEqual(get_output(self.seq, ['head']), ['>a','GATACA'])
        self.assertEqual(get_output(self.seq, ['head', '-n', '2']), ['>a','GATACA','>b','GALLIF'])

    def test_n(self):
        self.assertEqual(get_output(self.seq, ['tail']), ['>c','SPARTA'])
        self.assertEqual(get_output(self.seq, ['tail', '-n', '2']), ['>b','GALLIF','>c','SPARTA'])

    def test_fl(self):
        self.assertEqual(get_output(self.seq, ['head', '-f', 2, '-l', 1])[1], 'GA...A')
        self.assertEqual(get_output(self.seq, ['tail', '-f', 2, '-l', 1])[1], 'SP...A')

class TestPerm(unittest.TestCase):
    def setUp(self):
        self.seq = ['>a', 'WHEREISMYTARDIS']

    def test_default(self):
        self.assertEqual(get_output(
            self.seq,
            ['perm', '--seed', 42]),
            ['>a|PERMUTATION:start=0;end=0;word_size=1', 'MTISSADYHEIERWR'])

    def test_word_size(self):
        self.assertEqual(get_output(
            self.seq,
            ['perm', '--seed', 42, '-w', 3]),
            ['>a|PERMUTATION:start=0;end=0;word_size=3', 'TARREISMYDISWHE'])
        self.assertEqual(get_output(
            self.seq,
            ['perm', '--seed', 42, '-w', 5]),
            ['>a|PERMUTATION:start=0;end=0;word_size=5', 'ARDISISMYTWHERE'])

    def test_offsets(self):
        self.assertEqual(get_output(
            self.seq,
            ['perm', '--seed', 42, '-w', 4, '-s', 3]),
            ['>a|PERMUTATION:start=3;end=0;word_size=4', 'WHERDISMYTAREIS'])
        self.assertEqual(get_output(
            self.seq,
            ['perm', '--seed', 123, '-w', 4, '-s', 5]),
            ['>a|PERMUTATION:start=5;end=0;word_size=4', 'WHEREISTARDISMY'])
        self.assertEqual(get_output(
            self.seq,
            ['perm', '--seed', 123, '-w', 4, '-e', 3]),
            ['>a|PERMUTATION:start=0;end=3;word_size=4', 'YTAREISMWHERDIS'])

class TestRename(unittest.TestCase):
    def setUp(self):
        self.seq=[
            '>freddy','A',
            '>fred','A',
            '>bob','A']

    def test_replace(self):
        self.assertEqual(get_output(self.seq, ['rename', 'fred', 'a']), ['>ady', 'A', '>a', 'A', '>bob', 'A'])
        self.assertEqual(get_output(self.seq, ['rename', '[frb]', '']), ['>eddy', 'A', '>ed', 'A', '>o', 'A'])

    def test_replace_where(self):
        self.assertEqual(get_output(self.seq, ['rename', '$', '~', 'fred']), ['>freddy~', 'A', '>fred~', 'A', '>bob', 'A'])

class TestReverse(unittest.TestCase):
    def setUp(self):
        self.seq = [
            '>a1', 'LIVED',
            '>b2', 'MILLER'
        ]
        self.reverse = [
            '>a1|REVERSE', 'DEVIL',
            '>b2|REVERSE', 'RELLIM'
        ]
    def test_default(self):
        self.assertEqual(get_output(self.seq, ['reverse']), self.reverse)

class TestSample(unittest.TestCase):
    def setUp(self):
        self.seqs=[
            '>1', 'A',
            '>2', 'A',
            '>3', 'A',
            '>4', 'A',
            '>5', 'A']

    def test_default(self):
        self.assertEqual(get_output(self.seqs, ['sample', '--seed', '5']), ['>5', 'A'])
        self.assertEqual(get_output(self.seqs, ['sample', '--seed', '5', '2']), ['>5', 'A', '>3', 'A'])
        self.assertEqual(get_output(self.seqs, ['sample', '2', '--seed', '123']), ['>1', 'A', '>3', 'A'])

class TestSort(unittest.TestCase):
    def setUp(self):
        self.unsorted=[
            '>g=c;d=100','AAA',
            '>g=d;d=30','AA',
            '>g=b;d=9','AAAA',
            '>g=a;d=200','A',
        ]
        self.default=[
            '>g=a;d=200','A',
            '>g=b;d=9','AAAA',
            '>g=c;d=100','AAA',
            '>g=d;d=30','AA',
        ]
        self.default_reverse=[
            '>g=d;d=30','AA',
            '>g=c;d=100','AAA',
            '>g=b;d=9','AAAA',
            '>g=a;d=200','A',
        ]
        self.length=[
            '>g=a;d=200','A',
            '>g=d;d=30','AA',
            '>g=c;d=100','AAA',
            '>g=b;d=9','AAAA',
        ]
        self.regex=[
            '>g=c;d=100','AAA',
            '>g=a;d=200','A',
            '>g=d;d=30','AA',
            '>g=b;d=9','AAAA',
        ]
        self.regex_numeric=[
            '>g=b;d=9','AAAA',
            '>g=d;d=30','AA',
            '>g=c;d=100','AAA',
            '>g=a;d=200','A',
        ]

    def test_default(self):
        self.assertEqual(get_output(self.unsorted, ['sort']), self.default)

    def test_default(self):
        self.assertEqual(get_output(self.unsorted, ['sort', '-r']), self.default_reverse)

    def test_length_sort(self):
        self.assertEqual(get_output(self.unsorted, ['sort', '-l']), self.length)

    def test_regex_sort(self):
        self.assertEqual(get_output(self.unsorted, ['sort', '-x', 'd=(\d+)']), self.regex)

    def test_numeric_sort(self):
        self.assertEqual(get_output(self.unsorted, ['sort', '-x', 'd=(\d+)', '-n']), self.regex_numeric)

class TestSubseq(unittest.TestCase):
    def setUp(self):
        self.seq=['>a', 'GATACA']
        self.aaseq=['>p', 'PICKLE']

    def test_default(self):
        self.assertEqual(get_output(self.seq, ['subseq', '-b', 1, 1])[1], 'G')
        self.assertEqual(get_output(self.seq, ['subseq', '-b', 5, 6])[1], 'CA')

    def test_overbounds(self):
        self.assertEqual(get_output(self.seq, ['subseq', '-b', 1, 100])[1], 'GATACA')
        self.assertRaises(SystemExit, get_output, self.seq, ['subseq', '-b', 7, 7])

    def test_revcomp(self):
        self.assertEqual(get_output(self.seq, ['subseq', '-b', 3, 6])[1], 'TACA')
        self.assertEqual(get_output(self.aaseq, ['subseq', '-b', 1, 3])[1], 'PIC')
        # guess_type function adientifies GATACA as dna, then takes revcomp
        self.assertEqual(get_output(self.seq, ['subseq', '-b', 6, 3])[1], 'TGTA')
        # PICKLE however doesn't appear to be dna, so reversing does nothing
        self.assertEqual(get_output(self.aaseq, ['subseq', '-b', 3, 1])[1], 'PIC')

class TestUniq(unittest.TestCase):
    def setUp(self):
        self.all_uniq=[
            '>a','CAT',
            '>b','HAT',
            '>c','A']
        self.unsorted=[
            '>b','HAT',
            '>a','CAT',
            '>c','A']
        self.repeated=[
            '>a','CAT',
            '>b','HAT',
            '>c','A',
            '>b','HAT',
            '>c','A'
        ]

    def test_default(self):
        self.assertEqual(get_output(self.all_uniq, ['uniq']), self.all_uniq)
        self.assertEqual(get_output(self.repeated, ['uniq']), self.all_uniq)
        self.assertEqual(get_output(self.unsorted, ['uniq']), self.unsorted)

    def test_uniq(self):
        self.assertEqual(get_output(self.all_uniq, ['uniq', '-u']), self.all_uniq)
        self.assertEqual(get_output(self.repeated, ['uniq', '-u']), ['>a', 'CAT'])

    def test_duplicated(self):
        self.assertEqual(get_output(self.all_uniq, ['uniq', '-d']), [''])
        self.assertEqual(get_output(self.repeated, ['uniq', '-d']), ['>b', 'HAT',  '>c', 'A'])

    def test_count(self):
        self.assertEqual(get_output(self.all_uniq, ['uniq', '-c']), ['1\ta', '1\tb', '1\tc'])
        self.assertEqual(get_output(self.repeated, ['uniq', '-c']), ['1\ta', '2\tb',  '2\tc'])

class TestWc(unittest.TestCase):
    def setUp(self):
        self.seq=['>a','CAT',
                  '>b','HAT',
                  '>c','A']
    def test_default(self):
        self.assertEqual(get_output(self.seq, ['wc']), ['3\t7'])

    def test_nseqs(self):
        self.assertEqual(get_output(self.seq, ['wc', '-l']), ['3'])

    def test_nchars(self):
        self.assertEqual(get_output(self.seq, ['wc', '-m']), ['7'])

class TestWinnow(unittest.TestCase):
    def setUp(self):
        self.seq = [
            '>a', 'ASDFX',
            '>b', 'ASDF',
            '>c', 'ASD'
        ]

    def test_contain(self):
        self.assertEqual(get_output(self.seq, ['winnow', '-c', 'X'])[0::2], ['>b', '>c'])
        self.assertEqual(get_output(self.seq, ['winnow', '-c', 'XF'])[0::2], ['>c'])

    def test_not_contain(self):
        self.assertEqual(get_output(self.seq, ['winnow', '-vc', 'X'])[0::2], ['>a'])
        self.assertEqual(get_output(self.seq, ['winnow', '-vc', 'XF'])[0::2], ['>a', '>b'])

    def test_shorter_than(self):
        self.assertEqual(get_output(self.seq, ['winnow', '-s', 3])[0::2], ['>a', '>b', '>c'])
        self.assertEqual(get_output(self.seq, ['winnow', '-s', 4])[0::2], ['>a', '>b'])
        self.assertEqual(get_output(self.seq, ['winnow', '-s', 5])[0::2], ['>a'])

    def test_longer_than(self):
        self.assertEqual(get_output(self.seq, ['winnow', '-vs', 3])[0::2], [''])
        self.assertEqual(get_output(self.seq, ['winnow', '-vs', 4])[0::2], ['>c'])
        self.assertEqual(get_output(self.seq, ['winnow', '-vs', 5])[0::2], ['>b', '>c'])

    def test_composition(self):
        comp = [
            '>a', 'AAAAG.....',
            '>b', 'AG........'
        ]
        self.assertEqual(get_output(comp, ['winnow', '-p', 'AG < 1'])[0::2], [''])
        self.assertEqual(get_output(comp, ['winnow', '-p', 'AG <= 0.5'])[0::2], [''])
        self.assertEqual(get_output(comp, ['winnow', '-p', 'AG < 0.5'])[0::2], ['>a'])
        self.assertEqual(get_output(comp, ['winnow', '-p', 'AG < 0.2'])[0::2], ['>a', '>b'])


if __name__ == '__main__':
    unittest.main()
