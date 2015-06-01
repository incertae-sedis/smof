smof
====

A collection FASTA sequence functions with convenient command-line interface

Getting Help
============

Detailed instructions on how to use each command in `smof` is available via the
'-h' option.

To list subcommands

``` bash
smof -h
```

Get help on a specific subcommand

``` bash
smof grep -h
```

Installation
============

Currently smof has only been tested in a UNIX environment.

``` bash
git clone https://github.com/zbwrnz/smof
cd smof
./runtest.py
ln -s smof.py /usr/local/bin/smof
```

UNIX-like commands
==================

This group of subcommands include commands based off UNIX builtins.

## `smof grep`

Whereas GNU grep searches lines for matches, `smof grep` searches either the
FASTA headers OR the fasta sequence.

Extract the entry with the id AT5G49640

``` bash
smof grep AT1G10960 at.faa
```

Extract entries containing the sequence 'FFQQ'

```bash
smof grep -q FFQQ at.faa
```

Extract only the FFQQ motif and the 3 amino acids flanking it.
```bash
smof grep -qoC3 FFQQ at.faa
```

Write the output in gff format
```bash
smof grep -qC3 --gff FFQQ at.faa
```

Count occurences (on both strands) of a DNA pattern using IUPAC extended
nucleotide alphabet.
```bash
smof grep -qrG YYNCTATAWAWASM myfile.fna
  692
```

Find non-overlapping open reading frames of length greater than 100 codons.
This is meant as an example of regex searching. This will NOT give you a great
answer. smof does not consider frames (nor will it ever). It will not find the
set of longest possible ORFs. If you want to identify ORFs, you should use a
specialized program. That said:

``` bash
smof grep -qPr --gff 'ATG(.{3}){99,}?(TAA|TGA|TAG)' myfile.fna
  chr3    smof-1.19.0   regex_match   357   668   .  +  .  .
  chr3    smof-1.19.0   regex_match   823   1152  .  +  .  .
  chr3    smof-1.19.0   regex_match   1230  1568  .  +  .  .
```

A particularly powerful function of `smof grep` is the ability to read a whole
file of patterns and match them against a regex capture. This allows O(n),
rather than O(mn) as in GNU grep, extraction of entries containing a particular
pattern. For example if your headers are formatted like
'locus|xxx|taxon|yyy|gi|zzz' and you have a file of thousands of gi numbers,
you can easily extract all the sequences in the FASTA file matching one of
these gi numbers with the following command:

```bash
smof grep -w 'gi\|(\d+)' -f gi_numbers.txt seq.fa
```

## `smof md5sum`

This tool is useful if you want a checksum for a FASTA file that is independent
of format (e.g. column width).

## `smof head` and `tail`

These functions mimic their GNU counterparts but on the entry, rather than
line, level. For example `smof head` prints the first entry in a file and `smof
-5` prints the first 5. Similarly for `smof tail`. 

## `smof sort`

`smof sort` can be used to simply sort sequences alphabetically by header. It
can also sort by sequence length. One useful feature with no homolog in GNU
sort is the ability to sort by regex capture. For example, if the FASTA headers
are formated like 'locus|xxx|taxon|yyy|gi|zzz', you can sort them numerically
by taxon with the command `smof sort -nx 'taxon\|(\d+)'`.

## `smof sample`

`smof sample` allows extraction of a random sample of entries. With no
arguments, it reads the entire file into memory and outputs a random one.

## `smof split`

This command allows easily splitting of a large file into many smaller files.

## `smof uniq`

This is currently a pretty useless command.

## `smof wc`

Outputs the number of characters and entries in the fasta file.

String manipulation commands
============================

## `smof permute`

Permutes the letters of a sequence

## `smof reverse`

Reverses a sequence (does NOT take the reverse complement)

## `smof subseq`

Get the sequences with headers matching 'chr3', get the subsequence 357,668

``` bash
cat myfile.fna | smof grep chr3 | smof subseq -b 357 668 
  >chr3
  atggtcctttctcttgtttcttctctgtgttgttgagattagtttgtttaggtttgatagcgttgattttggcctgcgtt
  tggtgactcatatggtttgattggagtttgtttctgggttttatggttttggttgaagcgacatttttttgtggaatatg
  gtttttgcaaaatattttgttccggatgagtaatatctacggtgctgctgtgagaattatgctattgttttgcaggtcct
  gttcttaatctttcatcgcttttgtgcttattgtctccttgtcgtttatgttgagtggtgtttgggctttag
```

If the start is higher than the end, and the sequence appears to be a DNA
sequence, then smof will take the reverse complement.

`smof subseq` can also read from a gff file. However, if you want to do extract
many sequences from a fasta file using a gff file as a guide (or other gff/bed
manipulations), consider using a specialized tools such as 'bedtools'.


Biological sequence tools
=========================

## `smof clean`

This command can be used to tidy a sequence. You can change the column width,
remove gaps and stops, change sequence to letters, etc.

## `smof filter`

Output only sequence that meet a set of conditions.

If you want to only keep sequences that are longer than 100 letters

```bash
smof filter -l 100 myfile.fa
```

Or shorter than 100 letters

```bash
smof filter -s 100 myfile.fa
```

Or that have greater than 60% AFILMVW content (hydrophobic amino acids)

```bash
smof filter -c 'AFILMVW > .6' myfile.fa
```

## `smof sniff`

This command runs a number of checks on a FASTA file and is useful in
diagnosing problems. For details, run `smof sniff -h`.

## `smof stat`

Just run `smof stat -h`. This thing is awesome and can do all sorts of crazy
shit.
