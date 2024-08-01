#!/usr/bin/env bash

set -e
set -u

diff <(smof grep -P '(YP_008964013.1|YP_008964017.1|YP_008964016.1)' aa.faa) aa-sub.faa
diff <(smof grep -f ids.txt aa.faa) aa-sub.faa
diff <(smof grep --fastain aa-sub.faa aa.faa) aa-sub.faa

diff <(smof head -3 aa.faa) <(head -17 aa.faa)
diff <(smof tail -3 aa.faa) <(tail -13 aa.faa)

diff <(smof head -2 trna.fna | smof tail)    <(smof cut -f2 trna.fna)
diff <(smof tail trna.fna | smof head)       <(smof cut -f10 trna.fna)
diff <(smof tail -5 trna.fna | smof head -2) <(smof cut -f6-7 trna.fna)

smof consensus trna.aln > /dev/null
smof consensus -t trna.aln > /dev/null
