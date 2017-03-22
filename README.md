[![Build Status](https://travis-ci.org/arendsee/smof.svg?branch=master)](https://travis-ci.org/arendsee/smof)

smof
====

UNIX-style FASTA tools

For instructions on installation and a detailed tutorial, see the
[project page](http://zbwrnz.github.io/smof).

Much to my chagrine, I have found another suite of tools that overlaps heavily
with my smof: [FAST](https://github.com/tlawrence3/FAST). Their work was
published in May of 2015. The paper is well written and cogently argued for the
need of bioinformatics tools that follow UNIX philosophy.

> Lawrence, Travis J., et al. "FAST: FAST Analysis of Sequences Toolbox." Frontiers in Genetics 6 (2015): 172.

FAST is elegant, well documented, and absolutely worth studying. I will write a
comparison of smof and FAST at some point in the future.

Install
=======

```
git clone https://github.com/arendsee/smof
cd smof
cp -s $PWD/smof.py ~/bin/smof
```

You should replace `~/bin/smof` with some folder that is in PATH.

Screenshots
===========


![main help](doc/images/main-help.png)

![subcommand help](doc/images/grep-help.png)

![header grep](doc/images/grep-maturase.png)

![muscle gapped](doc/images/muscle-gapped.png)

![muscle gapped gff](doc/images/muscle-gapped-gff.png)

![subseq](doc/images/subseq.png)

![motif counts](doc/images/motif-counts.png)
