# seq-position-probs

This software finds similarities between genetic sequences and
"profiles".  A profile is a set of position-specific letter, deletion,
and insertion probabilities.

This isn't very useful: it's a proof-of-principle for the paper [A
simple way to find related sequences with position-specific
probabilities](https://doi.org/10.1101/2025.03.14.643233).

To compile it, just do `make`.

## Similarity scores between profiles and random sequences

    seq-position-probs profile.hmm numOfSequences sequenceLength

where `profile.hmm` has profile(s) in HMMER3/f format.  For each
profile versus each sequence, it shows three kinds of score:

* The maximum end-anchored score
* The maximum start-anchored score
* The maximum mid-anchored score

Here, an "anchor" means a pair of coordinates (*i*,*j*) in the profile
and the sequence.  It gets each kind of score for all possible
anchors, and shows the maximum.

It then estimates Gumbel parameters &lambda; and *K* for each kind of
score, assuming the scores are Gumbel-distributed.

* `lamMM` is &lambda; estimated by method-of-moments.
* `kMM` is *K* estimated by method-of-moments.
* `kMM1` is *K* estimated by method-of-moments assuming &lambda; = 1.
* `lamML` is &lambda; estimated by maximum-likelihood.
* `kML` is *K* estimated by maximum-likelihood.
* `kML1` is *K* estimated by maximum-likelihood assuming &lambda; = 1.

### Border

The `-b` option adds a "border" to each random sequence.  For example,
`-b100` adds a border of length 100.  If the sequence length is L and
border length is B, it first generates a random sequence of length L,
then appends B/L copies of this sequence to itself.  This aims to
avoid [edge effects][] on the distribution of scores.

## Searching profiles against real sequences

    seq-position-probs profile.hmm randomSeqNum randomSeqLength realSeqs.fasta

This prints the maximum end-, start-, and mid-anchored score, for each
profile versus each sequence.  It also shows an *E*-value for each
score.  This *E*-value means: the expected number of distinct
sequence regions with equal or higher score, if we compared *all* the
profiles to random sequences with the same length as *all* of
`realSeqs.fasta`.

For nucleotide profiles, it will search both strands of
`realSeqs.fasta`.  Use option `-s1` to search forward strands only, or
`-s0` for reverse strands only.

## Details

* The random sequences have letter frequencies equal to the profile's
  "background" letter frequencies.

* A profile's background letter frequencies are set to the geometric
  mean of its position-specific letter probabilities ([Barrett et
  al. 1997](https://doi.org/10.1093/bioinformatics/13.2.191)).

* A score is: log<sub>2</sub>[probability ratio].

* An *E*-value is: *K*<sub>tot</sub> *N* / 2^score, where
  *K*<sub>tot</sub> is the sum over profiles of `kML1`, and *N* is
  the sum of sequence lengths.

[edge effects]: https://doi.org/10.1093/nar/29.2.351
