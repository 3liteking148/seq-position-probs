#! /bin/sh

cd $(dirname $0)

PATH=../bin:$PATH

{
    seq-position-probs dfam-test.hmm 10 200
    seq-position-probs -s1 dfam-test.hmm 50 500 dna-test.fa
    seq-position-probs dfam-test.hmm 50 500 dna-test.fa
    seq-position-probs -s0 dfam-test.hmm 50 500 dna-test.fa
    sed '2s/^/N/' dna-test.fa | seq-position-probs -s1 dfam-test.hmm 50 500 -
    tr A N < dna-test.fa | seq-position-probs -s0 dfam-test.hmm 50 500 -
    seq-position-probs -e0.001 -b100 dfam-test.hmm 100 400 hg38-chr15-part.fa
} | diff -u tests.txt -
