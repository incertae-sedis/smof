#!/usr/bin/env bash

set -e
set -u

diff <(smof grep -P '(S16|S12|S13)' aa.faa) aa-sub.faa
diff <(smof grep -f pattern-list.txt aa.faa) aa-sub.faa
diff <(smof grep --fastain aa-sub.faa aa.faa) aa-sub.faa

diff <(smof head -3 aa.faa) <(head -17 aa.faa)
diff <(smof tail -3 aa.faa) <(tail -13 aa.faa)

diff <(smof subseq -b 1 3 aa.faa) <(smof subseq --gff <(smof grep --gff -qP '^...' aa.faa) aa.faa)

diff <(smof head -2 trna.fna | smof tail)    <(smof cut -f2 trna.fna)
diff <(smof tail trna.fna | smof head)       <(smof cut -f10 trna.fna)
diff <(smof tail -5 trna.fna | smof head -2) <(smof cut -f6-7 trna.fna)

smof consensus trna.aln > /dev/null
smof consensus -t trna.aln > /dev/null

    # clean               cleans fasta files
    # filter              extracts sequences meeting the given conditions
    # md5sum              calculate an md5 checksum for the input sequences
    # permute             randomly order sequence
    # reverse             reverse each sequence (or reverse complement)
    # sample              randomly select entries from fasta file
    # sniff               extract info about the sequence
    # sort                sort sequences
    # split               split a fasta file into smaller files
    # stat                calculate sequence statistics
    # translate           translate a DNA sequence into a protein sequence
