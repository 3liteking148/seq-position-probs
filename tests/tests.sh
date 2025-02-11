#! /bin/sh

cd $(dirname $0)

PATH=../bin:$PATH

{
    seq-position-probs dfam-test.hmm 10 200
    seq-position-probs dfam-test.hmm 50 500 dna-test.fa
} | diff -u tests.txt -
