#! /bin/sh

cd $(dirname $0)

PATH=../bin:$PATH

{
    seq-position-probs -t10 -l200 -b0 dfam-test.hmm
    seq-position-probs -t50 -b0 -e0 -s1 dfam-test.hmm dna-test.fa
    seq-position-probs -t50 -b0 -e0 dfam-test.hmm dna-test.fa
    seq-position-probs -t50 -b0 -e0 -s0 dfam-test.hmm dna-test.fa
    sed '2s/^/N/' dna-test.fa |
	seq-position-probs -t50 -b0 -e0 -s1 dfam-test.hmm -
    tr A N < dna-test.fa | seq-position-probs -t50 -b0 -e0 -s0 dfam-test.hmm -
    seq-position-probs -e0.001 -t100 -l400 dfam-test.hmm hg38-chr15-part.fa
} |
    grep -v seq-position-probs | diff -u tests.txt -
