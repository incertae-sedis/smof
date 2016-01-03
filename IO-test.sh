#!/bin/bash

# This script tests the rough correctness of IO
# Exhaustive unit testing is restricted to runtest.py

cd test-data

echo -n "Testing clean ... "
diff clean-test_-u_-w30_-r_-t_n_-x.fa \
     <(../smof.py clean -u -w 30 -rx -t n clean-test.fa) > /dev/null
echo $?

echo -n "Testing filter ... "
diff filter-test_-l5_GC50.fa \
     <(../smof.py filter -l5 -c 'GC > .5' filter-test.fa) > /dev/null
echo $?

echo -n "Testing grep ... "
diff grep-test_-q_DF.fa \
     <(../smof.py grep -q DF grep-test.fa) > /dev/null
echo $?

[[ md5sum ]] && {
    echo -n 'Testing md5sum ... '
    diff <(tr -d '>\n' < md5sum.fa | md5sum | awk '{print $1}') \
         <(../smof.py md5sum md5sum.fa)
    echo $?
}

echo -n 'Testing head ... '
diff <(head -2 grep-test.fa) <(../smof.py head grep-test.fa) > /dev/null
echo $?

echo -n 'Testing tail ... '
diff <(tail -2 grep-test.fa) <(../smof.py tail grep-test.fa) > /dev/null
echo $?

echo -n 'Testing multiple file input ... '
diff <( ../smof.py grep -q 'a' a.fa b.fa c.fa) <(cat a.fa b.fa c.fa | ../smof.py grep -q 'a') > /dev/null
echo $?

echo -n 'Testing --fastain option ... '
diff <(head -2 grep-test.fa) <(../smof.py grep --fastain b.fa grep-test.fa) > /dev/null
echo $?

echo -n 'Testing subseq coloring ... '
diff subseq-colored.fa \
    <(../smof.py subseq -b 1 5 -c blue -Y subseq.fa |
      ../smof.py subseq -b 7 7 -Yc green |
      ../smof.py subseq -b 2 3 -c red -Y |
      ../smof.py subseq -b 15 100 -c magenta -Y |
      ../smof.py subseq -b 10 20 -c green -Y |
      ../smof.py subseq -b 14 16 -c red) > /dev/null
echo $?

echo -n 'Testing uniq pack/unpack ... '
diff pack-unpack.fa \
    <(../smof.py uniq -p pack-unpack.fa | ../smof.py uniq -P)
echo $?
