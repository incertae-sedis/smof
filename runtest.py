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


if __name__ == '__main__':
    # a = ['>seq1','AcGGNttt','AGGGG','>seq2','atgat']
    # gen = getgen(a)
    unittest.main()
