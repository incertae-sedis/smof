[![stable](http://badges.github.io/stability-badges/dist/stable.svg)](http://github.com/badges/stability-badges)
[![Build Status](https://travis-ci.org/incertae-sedis/smof.svg?branch=master)](https://travis-ci.org/incertae-sedis/smof)
[![Docker Docker build](https://img.shields.io/docker/cloud/build/incertaesedis/smof.svg)](https://hub.docker.com/r/incertaesedis/smof/) [![docker pulls](https://img.shields.io/docker/pulls/incertaesedis/smof.svg)](https://hub.docker.com/r/incertaesedis/smof/)
![PyPI](https://img.shields.io/pypi/v/smof.svg)
[![DOI](https://zenodo.org/badge/19203682.svg)](https://zenodo.org/badge/latestdoi/19203682)

smof - Simple Manipulation Of FASTA
====

UNIX-style FASTA tools

Installation
============

```
pip install smof
```

Functions
=========

`smof` is divided into the following subcommands:

 | subcommand  | description                                           |
 | ----------  | ----------------------------------------------------- |
 | `cut`       | emulates UNIX cut command, where fields are entries   |
 | `clean`     | cleans fasta files                                    |
 | `consensus` | finds the consensus sequence for aligned sequence     |
 | `filter`    | extracts sequences meeting the given conditions       |
 | `grep`      | roughly emulates the UNIX grep command                |
 | `md5sum`    | calculate an md5 checksum for the input sequences     |
 | `head`      | writes the first sequences in a file                  |
 | `permute`   | randomly order sequence                               |
 | `reverse`   | reverse each sequence (or reverse complement)         |
 | `sample`    | randomly select entries from fasta file               |
 | `sniff`     | extract info about the sequence                       |
 | `sort`      | sort sequences                                        |
 | `split`     | split a fasta file into smaller files                 |
 | `stat`      | calculate sequence statistics                         |
 | `subseq`    | extract subsequence from each entry (revcomp if a\<b) |
 | `tail`      | writes the last sequences in a file                   |
 | `translate` | translate a DNA sequence into a protein sequence      |
 | `uniq`      | count, omit, or merge repeated entries                |
 | `wc`        | roughly emulates the UNIX wc command                  |


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

Sample data
===========

The FASTA files used in the examples below are available in the
`sample-data/anncaliia_algerae` folder in the `smof` github repo
([here](https://github.com/incertae-sedis/smof)).

UNIX-like commands
==================

This group of subcommands include commands based off UNIX builtins.

## `smof head` and `tail`

These functions mimic their GNU counterparts but on the entry, rather than
line, level. For example `smof head` prints the first entry in a file and `smof
-5` prints the first 5. Similarly for `smof tail`. 

```bash
smof head aa.faa
smof head -3 aa.faa
smof tail aa.faa
smof tail -3 aa.faa
smof tail +2 aa.faa | smof head
```

In addition to the GNU-like functionallity, `smof head` and `tail` can also
limit the sequence that is output. This can be useful for diagnostic purposes.

```bash
# print last 3 nucleotides (last codon) from the first 5 transcripts
smof head -l 3 -5 aa.transcripts.fna
# print the first codon
smof head -f 3 -5 aa.transcripts.fna
# print first and last
smof head -f 3 -l 3 -5 aa.transcripts.fna
```

This sort of diagnostics is easier done with `smof sniff`.

## `smof sort`

`smof sort` can be used to simply sort sequences alphabetically by header. It
can also sort by sequence length. One useful feature with no homolog in GNU
sort is the ability to sort by regex capture. For example, if the FASTA headers
are formated like 'locus|xxx|taxon|yyy|gi|zzz', you can sort them numerically
by taxon with the command `smof sort -nx 'taxon\|(\d+)'`.

```bash
# print the shortest sequence
smof sort -l aa.faa | smof head
# print the longest sequence
smof sort -l aa.faa | smof tail
# sort by the function in the header description
smof sort -x 'PRA339 (.*)' aa.faa | smof tail
```

## `smof sample`

`smof sample` allows extraction of a random sample of entries. With no
arguments, it reads the entire file into memory and outputs a random one.

```bash
# retrieve 1 sequence by default
smof sample aa.faa
smof sample -n 5 aa.faa
# set a random seed (useful for debugging and reproducible scripts)
smof sample --seed 42 aa.faa
```

## `smof split`

This command allows easily splitting of a large file into many smaller files.

You can split a large file several small files with equal numbers of sequences
```bash
smof split -n 5 -p zzz aa.faa
grep -c '>' aa.faa zzz*
rm zzz*
```

Of you can split a large file into many smaller files with a set maximum number
of sequences per file
```bash
smof split -qn 500 -p zzz aa.faa
grep -c '>' aa.faa zzz*
rm zzz*
```

## `smof uniq`

This command corresponds roughly to GNU uniq, but entries are considered
identical only if both header and sequence are exactly the same. As currently
implemented, I don't find much use for this command.

## `smof wc`

Outputs the number of characters and entries in the fasta file. Generally `smof
stat` is better.

## `smof grep`

Whereas GNU grep searches lines for matches, `smof grep` searches either the
FASTA headers or the FASTA sequence.

Extract the entries by matches to the header (default)

``` bash
smof grep H312_03353 aa.faa
```

Extract entries by matches to a sequence 

```bash
smof grep --match-sequence SKSQ aa.faa
# or equivalently
smof grep -q SKSQ aa.faa
```

You can include flanking regions in the match
```bash
# match 3 residues downstream
smof grep -qA3 'SKSQ' aa.faa
# match 3 residues upstream 
smof grep -qB3 'SKSQ' aa.faa
# match 3 residues up- and downstream 
smof grep -qC3 'SKSQ' aa.faa
```

Inclusion of flanking regions is particularly useful in tandem with the -o
option, which extracts only the matching sequence
```bash
smof grep -qoA3 'SKSQ' aa.faa
```

Write the output in gff format
```bash
smof grep -q --gff SKSQ aa.faa
```

You can count the number of sequences with a match
```bash
smof grep -qc SKS aa.faa
```

Or the total number of matches
```bash
smof grep -qm SKSQ aa.faa
```

Or both
```bash
smof grep -qmc SKS aa.faa
```

Just like in GNU grep, you can invert a search. This search finds all genes
that are not annotated as being hypothetical genes.
```bash
smof grep -v hypothetical aa.faa
```

By default `smof grep` is case insensitive (unlike GNU grep), but it can be
made case sensitive
```bash
smof grep -I CoA aa.faa
```

You can search using patterns in a file
```bash
smof grep -f id-sample.txt aa.faa
```

This, however, can be a little slow, since it searchs each pattern in the file
against the entire header. A much faster approach is to extract a search
pattern from the headers (or sequence) and then lookup the header pattern in
the set of search patterns.

```bash
smof grep -f id-sample.txt -w '\| (\S+) \|' aa.faa
```

Count occurrences (on both strands) of a DNA pattern using IUPAC extended
nucleotide alphabet.
```bash
smof grep -qmbG YYNCTATAWAWASM aa.supercontigs.fna
```

You can search using a sequence query
```bash
# select 5 random sequences
smof sample -n 5 aa.faa | smof subseq -b 5 35 > rand.faa
smof grep -q --fastain rand.faa aa.faa
```

Or you can search for identical sequences shared between two fasta files
```bash
smof sample -n 5 aa.faa > rand.faa
smof grep -q --fastain rand.faa aa.faa 
```

Find non-overlapping open reading frames of length greater than 100 codons.
This is meant as an example of regex searching. This will NOT give you a great
answer. smof does not consider frames (nor will it ever). It will not find the
set of longest possible ORFs. If you want to identify ORFs, you should use a
specialized program. That said:

``` bash
smof grep -qPb --gff 'ATG(.{3}){99,}?(TAA|TGA|TAG)' aa.supercontigs.fna
```

## `smof md5sum`

This tool is useful if you want a checksum for a FASTA file that is independent
of format (e.g. column width or case).


String manipulation commands
============================

## `smof permute`

Permutes the letters of a sequence

## `smof reverse`

Reverses a sequence (does NOT take the reverse complement)

## `smof subseq`

``` bash
# extract a subsequence
smof grep H312_00003T0 aa.faa | smof subseq -b 10 20
# color the subsequences instead
smof grep H312_00003T0 aa.faa | smof subseq -b 10 20 -c red
```

If the start is higher than the end, and the sequence appears to be a DNA
sequence, then smof will take the reverse complement.

`smof subseq` can also read from a gff file. However, if you want to extract
many sequences from a fasta file using a gff file as a guide (or other gff/bed
manipulations), consider using a specialized tools such as `bedtools`.


Biological sequence tools
=========================

## `smof clean`

This command can be used to tidy a sequence. You can change the column width,
remove gaps and stops, convert all letters to one case and/or change irregular
characters to unknowns. By default, it removes whitespace in a sequence and
makes uniform, 80-character columns.

## `smof filter`

Output only sequence that meet a set of conditions.

If you want to only keep sequences that are longer than 100 letters

```bash
smof clean -x aa.faa | smof filter -l 100
```

Note that I call clean before filtering to remove the stop character, which
should not be included when calculating length.

Or shorter than 100 letters

```bash
smof clean -x aa.faa | smof filter -s 100 aa.faa
```

Or that have greater than 50% AFILMVW content (hydrophobic amino acids)

```bash
smof clean -x aa.faa | smof filter -c 'AFILMVW > .5' aa.faa
```

## `smof sniff`

This command runs a number of checks on a FASTA file and is useful in
diagnosing problems. For details, run `smof sniff -h`.

## `smof stat`

The default operation outputs number of sequences and summary statistics
concerning the sequence lengths.

```bash
smof stat aa.supercontigs.fna
 nseq:      431
 nchars:    12163397
 5sum:      445 3301 9555 30563 746881
 mean(sd):  28221 (58445)
 N50:       71704
```

'5sum' refers to the five number summary of the sequence lengths (minimum, 25%
quantile, median, 75% quantile, and maximum).

Statistics can also be calculated on a sequence-by-sequence level, which by
default outputs the sequence names (the first word of the header) and the
sequence length, e.g.

```bash
smof stat -q aa.supercontigs.fna | head
```

There are many other options. Run `smof stat -h` for descriptions.


Case study: exploring motifs in chloroplast genomes
===================================================

Alice is interested in the chloroplast *maturase* gene. Bob gives her a sample
dataset which includes 10 fasta files of proteins encoded by the chloroplast
genomes of 10 different plant species. These files are available in the
`sample-data/chloroplasts` directory.

You can find this dataset in the folder *doc/test-data/chloroplast-proteins*.

Her first step is to explore the data. She first counts the sequences in each
file with a simple grep command.

```
grep -c '>' *faa
```

Next she tests the sequences with `smof sniff`

```
smof sniff *faa
```

Producing the following output:

```
578 uniq sequences (757 total)
All prot
All uppercase
Protein Features:
  initial-Met:         755        99.7358%
  terminal-stop:       0          0.0000%
  internal-stop:       0          0.0000%
  selenocysteine:      0          0.0000%
Universal Features:
  unknown:             8          1.0568%
  ambiguous:           0          0.0000%
  gapped:              0          0.0000%
```

Everything looks pretty good. But two of the sequences don't start with a
methionine. Alice wants to find them. She does this using `smof grep` and a
Perl regular expressions.

```
smof grep -qP '^[^M]' *faa
```

She finds these genes are both from *Solanum lycopersicum* and are described in
the fasta headers as being *partial*.

Now Alice wants to find the *maturase* genes by pulling out every entry with
'maturase' in the fasta header.

```
smof grep maturase *faa
smof grep maturase *faa > maturase.faa
```

For a close look at the distribution of sequence lengths, Alice calls `smof
stat`

```
smof stat maturase.faa
```

Alice happens to be interested in the sequence WTQPQR from *Panicum virgatum*
and would like to know what the homologous regions are in the other species.

So Alice aligns the maturase genes with
[MUSCLE](http://nar.oxfordjournals.org/content/32/5/1792.short) and searches
for the motif using the GFF output option.

```
muscle -quiet < maturase.faa | tee maturase.aln | smof grep -q --gff WTQPQR
```

This is outputs the location of the match in standard GFF format, i.e. the
match is at position 329 to 334. Homologs to this sequence will be at the same
positions in the aligned fasta file output by MUSCLE.

```
smof subseq -b 329 334 maturase.aln
```

HMMER could then be used to analyze the by-site conservation of the sextuplet.
