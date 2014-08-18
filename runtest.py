#!/usr/bin/env python3

import smof
import unittest
import argparse
import sys
from tempfile import TemporaryFile

class TestFSeq(unittest.TestCase):
    def test_fseq_subseq(self):
        header = 'seq1'
        seq = 'AcGGNttt'
        seqobj = smof.FSeq(header, seq)
        sq = seqobj.subseq(1, 5)
        self.assertEqual(sq.seq, 'cGGN')
        self.assertEqual(sq.header, 'seq1|SUBSEQ(2..5)')

    def test_fseq_getrevcomp_fromStringInput(self):
        seq = 'ACGTT'
        self.assertEqual(smof.FSeq.getrevcomp(seq), 'AACGT')

    def test_fseq_getrevcomp_fromFSeqInput(self):
        header = 'seq1'
        seq = 'ACGTT'
        seqobj = smof.FSeq(header,seq)
        rc = smof.FSeq.getrevcomp(seqobj)
        self.assertEqual(rc.seq, 'AACGT')
        self.assertEqual(rc.header, 'seq1|REVCOM')

    def test_fseq_ungap(self):
        header = 'seq1'
        seq = 'A.C-G_T'
        seqobj = smof.FSeq(header,seq)
        seqobj.ungap()
        self.assertEqual(seqobj.seq, 'ACGT')
        self.assertEqual(seqobj.header, 'seq1|UNGAPPED')

    def test_fseq_reverse(self):
        header = 'seq1'
        seq = 'ACGTT'
        seqobj = smof.FSeq(header,seq)
        seqobj.reverse()
        self.assertEqual(seqobj.seq, 'TTGCA')
        self.assertEqual(seqobj.header, 'seq1|REVERSE')


if __name__ == '__main__':
    # a = ['>seq1','AcGGNttt','AGGGG','>seq2','atgat']
    # gen = getgen(a)
    unittest.main()
