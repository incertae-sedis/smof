#!/usr/bin/env python3

import smof
import unittest
import argparse
import sys
from tempfile import TemporaryFile

class TestSmof(unittest.TestCase):
    def test_fseq(self):
        header = 'seq1'
        seq = 'AcGGNttt'
        seqobj = smof.FSeq(header,seq)

        self.assertEqual(seqobj.countmasked(), 4)
        self.assertEqual(seqobj.subseq(1,5),'cGGN')
        self.assertEqual(seqobj.getvalue('header'), 'seq1')
        self.assertEqual(smof.FSeq.getrevcomp(seqobj.seq), 'aaaNCCgT')


if __name__ == '__main__':
    # a = ['>seq1','AcGGNttt','AGGGG','>seq2','atgat']
    # gen = getgen(a)
    unittest.main()
