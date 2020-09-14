#!/usr/bin/env bash

set -e
set -u

diff <(smof grep -f ids.txt aa.faa) aa-sub.faa
diff <(smof grep --fastain aa-sub.faa aa.faa) aa-sub.faa

diff <(smof head -3 aa.faa) <(head -17 aa.faa)
diff <(smof tail -3 aa.faa) <(tail -13 aa.faa)
