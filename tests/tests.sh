#! /bin/sh

cd $(dirname $0)

PATH=../bin:$PATH

seq-position-probs dfam-test.hmm 10 200 | diff -u tests.txt -
