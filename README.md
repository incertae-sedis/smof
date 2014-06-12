smof
====

Collection FASTA sequence functions with convenient command-line interface

Examples
========

Get help on all subcommands

``` bash
$ smof -h
```

Get help on a specific subcommand

``` bash
$ smof grep -h
```

Count occurences of a DNA pattern using IUPAC extended nucleotide alphabet
(note short-option clustering is possible, so rather than write '-m -B', we can
write '-mB').

``` bash
$ cat myfile.fna | smof grep -mB YYNCTATAWAWASM
692
```

Write positions of these patterns (on both strands) to a gff3 file

``` bash
$ cat myfile.fna | smof grep -Br --gff YYNCTATAWAWASM > TATA-box.gff3
```

Find non-overlapping open reading frames of length greater than 100 codons.
This is meant as an example of regex searching. This will NOT give you a great
answer. smof does not consider frames (yet). It will not find the set of
longest possible ORFs. If you want to identify ORFs, you should use a
specialized program. That said:

``` bash
$ cat myfile.fna |
  smof grep -P -r --gff 'ATG(.{3}){99,}?(TAA|TGA|TAG)' | head 3
chr3    smof-1.2.1   regex_match   357   668   .  +  .  .
chr3    smof-1.2.1   regex_match   823   1152  .  +  .  .
chr3    smof-1.2.1   regex_match   1230  1568  .  +  .  .
```

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
