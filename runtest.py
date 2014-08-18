#!/usr/bin/env python3

import smof
import unittest
import argparse
import sys
from tempfile import TemporaryFile

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



if __name__ == '__main__':
    # a = ['>seq1','AcGGNttt','AGGGG','>seq2','atgat']
    # gen = getgen(a)
    unittest.main()
