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
    sed 's/aaaaaaaaaa/nnnnnnnnnn/' hg38-chr15-part.fa |
	dummerl -e0.01 -t100 -l400 dfam-test.hmm -
    dummer -t100 -l200 -e0.01 DF000001253.hmm chr22problem.fa
    dummer PF05369.hmm mtmb1.fa

    cat hakoLTR.stk Notch.stk | dummer-build --countonly -
    dummer-build --countonly --symfrac=0.15 Notch.stk
    dummer-build --countonly --enone hakoLTR.stk
    cat hakoLTR.stk Notch.stk | dummer-build --maxiter=5 -
    dummer-build --countonly PF05369.alignment.seed
    dummer-build --maxdiff 0.02 ABC_tran_2.stk
    dummer-build --countonly --enone --dmix hakoLTR.mixdchlet --gapprior my.gapPrior hakoLTR.stk
    dummer-build --counts --maxiter=5 hakoLTR.stk
} |
    grep -v DUMMER | diff -u tests.txt -
