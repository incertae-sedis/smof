#!/usr/bin/env python3

import smof
import unittest
import argparse
import sys
from tempfile import TemporaryFile

def getgen(t):
    f1 = TemporaryFile()
    f1.write(bytes('\n'.join(t), 'ascii'))
    gen = smof.FSeqGenerator(fh=open(f1.name, 'rb'))
    return(gen)

def gettuple(f):
    with open(f, 'r') as f:
        return([x for x in f])

class TestSmof(unittest.TestCase):
    def test_fseq(self):
        header = 'seq1'
        seq = 'AcGGNttt'
        seqobj = smof.FSeq(header,seq)

        self.assertEqual(seqobj.countmasked(), 4)
        self.assertEqual(seqobj.subseq(1,5),'cGGN')
        self.assertEqual(seqobj.getvalue('header'), 'seq1')
        self.assertEqual(smof.FSeq.getrevcomp(seqobj.seq), 'aaaNCCgT')

    def test_fasta2csv(self):
        a = ['>seq1','AcGGNttt','AGGGG','>seq2','atgat']
        gen = getgen(a)
        args = smof.parse(['fasta2csv'])


if __name__ == '__main__':
    # a = ['>seq1','AcGGNttt','AGGGG','>seq2','atgat']
    # gen = getgen(a)
    unittest.main()
