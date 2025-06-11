#! /bin/sh

cd $(dirname $0)

PATH=../bin:$PATH

{
    seq-position-probs -n10 -l200 -b0 dfam-test.hmm
    seq-position-probs -n50 -b0 -s1 dfam-test.hmm dna-test.fa
    seq-position-probs -n50 -b0 dfam-test.hmm dna-test.fa
    seq-position-probs -n50 -b0 -s0 dfam-test.hmm dna-test.fa
    sed '2s/^/N/' dna-test.fa | seq-position-probs -n50 -b0 -s1 dfam-test.hmm -
    tr A N < dna-test.fa | seq-position-probs -n50 -b0 -s0 dfam-test.hmm -
    seq-position-probs -e0.001 -n100 -l400 dfam-test.hmm hg38-chr15-part.fa
} | diff -u tests.txt -
