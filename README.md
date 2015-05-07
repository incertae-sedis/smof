smof
====

A collection FASTA sequence functions with convenient command-line interface

Getting Help
============

Detailed instructions on how to use each command in `smof` is available via the
'-h' option.

To list subcommands

``` bash
$ smof -h
```

Get help on a specific subcommand

``` bash
$ smof grep -h
```

UNIX-like commands
==================

This group of subcommands include commands based off UNIX builtins.

## `smof grep`

Whereas GNU grep searches lines for matches, `smof grep` searches either the
FASTA headers OR the fasta sequence.

Extract the entry with the id AT5G49640

``` bash
smof grep AT1G10960 < at.faa
```

Extract entries containing the sequence 'FFQQ'

```bash
smof grep -q FFQQ < at.faa
```

Extract only the FFQQ motif and the 3 amino acids flanking it.
```bash
smof grep -qoC3 FFQQ < at.faa
```

Write the output in gff format
```bash
smof grep -qC3 --gff FFQQ < at.faa
```

Count occurences (on both strands) of a DNA pattern using IUPAC extended
nucleotide alphabet.
``` bash
$ smof grep -qrG YYNCTATAWAWASM < myfile.fna
692
```

Find non-overlapping open reading frames of length greater than 100 codons.
This is meant as an example of regex searching. This will NOT give you a great
answer. smof does not consider frames (nor will it ever). It will not find the
set of longest possible ORFs. If you want to identify ORFs, you should use a
specialized program. That said:

``` bash
$ smof grep -qPr --gff 'ATG(.{3}){99,}?(TAA|TGA|TAG)' < myfile.fna
chr3    smof-1.19.0   regex_match   357   668   .  +  .  .
chr3    smof-1.19.0   regex_match   823   1152  .  +  .  .
chr3    smof-1.19.0   regex_match   1230  1568  .  +  .  .
```

## `smof md5sum`

## `smof head` and `tail`

## `smof sort`

## `smof split`

## `smof sample`

## `smof uniq`

## `smof wc`

String manipulation commands
============================

## `smof permute`

## `smof reverse`

## `smof subseq`

Get the sequences with headers matching 'chr3', get the subsequence 357,668

``` bash
$ cat myfile.fna | smof grep chr3 | smof subseq 357 668 
>chr3
atggtcctttctcttgtttcttctctgtgttgttgagattagtttgtttaggtttgatagcgttgattttggcctgcgtt
tggtgactcatatggtttgattggagtttgtttctgggttttatggttttggttgaagcgacatttttttgtggaatatg
gtttttgcaaaatattttgttccggatgagtaatatctacggtgctgctgtgagaattatgctattgttttgcaggtcct
gttcttaatctttcatcgcttttgtgcttattgtctccttgtcgtttatgttgagtggtgtttgggctttag
```

If you want to do extract many sequences from a fasta file using a gff file as
a guide (or other gff/bed manipulations), consider using a specialized tools
such as 'bedtools'.

Biological sequence tools
=========================

## clean

## filter

## sniff 

## stat
