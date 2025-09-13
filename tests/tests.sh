#! /bin/sh

cd $(dirname $0)

PATH=../bin:$PATH

{
    dummer -t10 -l200 -b0 dfam-test.hmm
    dummer -t50 -b0 -e0 -s1 dfam-test.hmm dna-test.fa
    dummer -t50 -b0 -e0 dfam-test.hmm dna-test.fa
    dummer -t50 -b0 -e0 -s0 dfam-test.hmm dna-test.fa
    sed '2s/^/N/' dna-test.fa |
	dummer -t50 -b0 -e0 -s1 dfam-test.hmm -
    dummer -e0.001 -t100 -l400 dfam-test.hmm hg38-chr15-part.fa
    dummer PF05369.hmm mtmb1.fa

    cat hakoLTR.stk Notch.stk | dummer-build --countonly -
    dummer-build --countonly --symfrac=0.15 Notch.stk
    dummer-build --countonly --enone hakoLTR.stk
    cat hakoLTR.stk Notch.stk | dummer-build --maxiter=5 -
} |
    grep -v DUMMER | diff -u tests.txt -
