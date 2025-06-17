# seq-position-probs

This software finds similarities between genetic sequences and
"profiles".  A profile is a set of position-specific letter, deletion,
and insertion probabilities.

This is a proof-of-principle for the paper [A simple way to find
related sequences with position-specific probabilities][frith2025].

To compile it, just do `make`.

It can compare sequences in FASTA format to profiles in HMMER3/f format:

    seq-position-probs profiles.hmm sequences.fasta

The output shows similar regions in [MAF][] format:

    a score=44.6 E=8.15e-05 anchor=227,9133
    s myProfile   195 61 +   262 atacatgtaaagcgcttagaacagtgcctggcacatagtaagcgctcaataaatgttagct
    s mySequence 9102 60 + 10000 ctacaccttcaggacttagg-ctgtgcctggcacatagtaggtgctcagtagacactggtt

* The `score` reflects the likelihood that these regions are related
  rather than random.

* The *E*-value (`E=`) is the expected number of distinct sequence
  regions with equal or higher score, if we compared all profiles in
  `profiles.hmm` to random sequences with the same length as all of
  `sequences.fasta`.

* The `anchor` shows profile,sequence coordinates.  It means there are
  similar regions around these coordinates.  (These are
  [mid-anchored][frith2025] similarities.)

The `s` lines show a representative alignment.  This aligns letters
whose probability of being aligned is > 0.5, among all possible
alignments with that anchor.

## Options

Show all options and default values:

    seq-position-probs --help

Compare each profile to just the forward strand of each DNA sequence:

    seq-position-probs -s1 profiles.hmm sequences.fasta

`-s0` means reverse strands only, `-s1` means forward strands only,
and `-s2` means both strands (the default).  This is ignored for
proteins.

Get similarities with *E*-value at most (say) 0.01:

    seq-position-probs -e0.01 profiles.hmm sequences.fasta

`-e0` has a special meaning: for each profile versus each strand, it
shows the maximum [end-anchored][frith2025],
[start-anchored][frith2025], and [mid-anchored][frith2025] scores (in
that order).  It gets each kind of score for all possible anchors, and
shows the maximum.

## Random sequences

To calculate *E*-values, it needs to estimate a *K* parameter for each
profile.  To do that, it compares the profile to random sequences.  To
see details of this, give it a profile file only:

    seq-position-probs profiles.hmm

For each profile versus each sequence, it shows the maximum
[end-anchored][frith2025], [start-anchored][frith2025], and
[mid-anchored][frith2025] scores.

For each kind of score, it then shows various estimates of *K* (and
another &lambda; parameter that isn't used currently).

* `lamMM` is &lambda; estimated by the method of moments.
* `kMM` is *K* estimated by the method of moments.
* `kMM1` is *K* estimated by the method of moments assuming &lambda; = 1.
* `lamML` is &lambda; estimated by maximum-likelihood.
* `kML` is *K* estimated by maximum-likelihood.
* `kML1` is *K* estimated by maximum-likelihood assuming &lambda; = 1.
* `lamLM` is &lambda; estimated by the method of L-moments.
* `kLM` is *K* estimated by the the method of L-moments.

(`kLM1` isn't shown, because it's identical to `kMM1`.)

Only `kMM1` is used to calculate *E*-values.

These options affect the random sequences.

- `-t T`, `--trials T`: generate this many random sequences.

- `-l L`, `--length L`: length of each random sequence.

- `-b B`, `--border B`: add this size border to each random sequence.
  It first generates a random sequence of length L, then appends B/L
  copies of this sequence to itself.  This aims to avoid [edge
  effects][] on the distribution of scores.

## Details

* The random sequences have letter frequencies equal to the profile's
  "background" letter frequencies.

* A profile's background letter frequencies are set to the geometric
  mean of its position-specific letter probabilities ([Barrett et
  al. 1997](https://doi.org/10.1093/bioinformatics/13.2.191)).

* A score is: log<sub>2</sub>[probability ratio].

* An *E*-value is: *K*<sub>tot</sub> *N* / 2^score, where
  *K*<sub>tot</sub> is the sum over profiles of *K*, and *N* is the
  sum of sequence lengths.

[frith2025]: https://doi.org/10.1101/2025.03.14.643233
[edge effects]: https://doi.org/10.1093/nar/29.2.351
[MAF]: https://genome.ucsc.edu/FAQ/FAQformat.html#format5
