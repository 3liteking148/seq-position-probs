# seq-position-probs

This software finds similarities between genetic sequences and
"profiles".  A profile is a set of position-specific letter, deletion,
and insertion probabilities.

This isn't very useful: it's a proof-of-principle for the paper "A
simple way to find related sequences with position-specific
probabilities".

To compile it, just do `make`.

## Similarity scores between profiles and random sequences

    seq-position-probs profile.hmm numOfSequences sequenceLength

where `profile.hmm` has profile(s) in HMMER3/f format.  It calculates
three kinds of score:

* The maximum end-anchored score
* The maximum start-anchored score
* The maximum mid-anchored score

It then estimates a Gumbel *K* parameter for each kind of score.  This
becomes more accurate if you increase numOfSequences and
sequenceLength.

* `Kharm` is the maximum-likelihood estimate, which uses the harmonic
  mean of exponentiated scores ([Frith
  2024](https://pubmed.ncbi.nlm.nih.gov/39152037/)).

* `Kgeom` is the method-of-moments estimate, which uses the geometric
  mean of exponentiated scores ([Yu et
  al. 2002](https://doi.org/10.1093/bioinformatics/18.6.864)).

## Searching profiles against real sequences

    seq-position-probs profile.hmm randomSeqNum randomSeqLength realSeqs.fasta

This prints the maximum end-, start-, and mid-anchored score, for each
profile versus each sequence.  It also shows an *E*-value for each
score.  This *E*-value means: the expected number of distinct
similarities with equal or higher score, if we compared *all* the
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
  *K*<sub>tot</sub> is the sum over profiles of `Kharm`, and *N* is
  the sum of sequence lengths.
