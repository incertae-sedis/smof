#!/usr/bin/env bash

# This script tests the rough correctness of IO
# Exhaustive unit testing is restricted to runtest.py

total_failed=0

runtest () {
    msg=$1
    exp=$2
    obs=$3
    echo -n "Testing $msg ... "
    diff $2 $3 > /dev/null
    if [[ $? == 0 ]]; then
        echo "OK"
    else
        echo "FAIL"
        total_failed=$(( $total_failed + 1 ))
    fi
}

cd test-data

runtest "clean" \
        clean-test_-u_-w30_-r_-t_n_-x.fa \
        <(../smof.py clean -u -w 30 -rx -t n clean-test.fa)

runtest "filter" \
        filter-test_-l5_GC50.fa \
        <(../smof.py filter -l5 -c 'GC > .5' filter-test.fa) > /dev/null

runtest "grep" \
        grep-test_-q_DF.fa \
        <(../smof.py grep -q DF grep-test.fa)

[[ md5sum ]] && {
    runtest 'md5sum' \
            <(tr -d '>\n' < md5sum.fa | md5sum | awk '{print $1}') \
            <(../smof.py md5sum md5sum.fa)
}

runtest 'head' \
        <(head -2 grep-test.fa) \
        <(../smof.py head grep-test.fa)

runtest "tail" \
        <(tail -2 grep-test.fa) \
        <(../smof.py tail grep-test.fa)

runtest 'multiple file input' \
        <( ../smof.py grep -q 'a' a.fa b.fa c.fa) \
        <(cat a.fa b.fa c.fa | ../smof.py grep -q 'a')

runtest '--fastain option' \
        <(head -2 grep-test.fa) \
        <(../smof.py grep --fastain b.fa grep-test.fa)

runtest 'subseq coloring' \
        subseq-colored.fa \
        <(../smof.py subseq -b 1 5 -c blue -Y subseq.fa |
          ../smof.py subseq -b 7 7 -Yc green |
          ../smof.py subseq -b 2 3 -c red -Y |
          ../smof.py subseq -b 15 100 -c magenta -Y |
          ../smof.py subseq -b 10 20 -c green -Y |
          ../smof.py subseq -b 14 16 -c red -Y)

runtest 'uniq pack/unpack' \
        pack-unpack.fa \
        <(../smof.py uniq -p pack-unpack.fa | ../smof.py uniq -P)

# echo -n 'Testing grep file options -ql ... '
# diff grep_-ql_G.fa <(../smof.py grep -ql G [123].fa)
# echo $?
#
# echo -n 'Testing grep file options -qL ... '
# diff grep_-qL_G.fa <(../smof.py grep -qL G [123].fa)
# echo $?
#
# echo -n 'Testing grep file options -qc ... '
# diff grep_-qc_G.fa <(../smof.py grep -qc G [123].fa)
# echo $?
#
# echo -n 'Testing grep file options -qm ... '
# diff grep_-qm_G.fa <(../smof.py grep -qm G [123].fa)
# echo $?
#
# echo -n 'Testing grep file options -qcm ... '
# diff grep_-qcm_G.fa <(../smof.py grep -qcm G [123].fa)
# echo $?

if [[ $total_failed -eq 0 ]]
then
    exit 0
else
    exit 1
fi
