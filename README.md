# DUMMER

DUMMER (Dumb Uncomplicated Match ModelER) aims to find distant
relationships between genetic sequences (nucleotide or protein).  It's
similar to [HMMER][], but much simpler, and aspires to be better.

In more detail, it finds similar regions between sequences and
"profiles".  A profile is a set of position-specific letter, deletion,
and insertion probabilities: typically made from a family of related
sequences.

This is a proof-of-principle for the paper [Simple and thorough
detection of related sequences with position-varying probabilities of
substitutions, insertions, and deletions][frith2025].

## Motivation: HMMER &rarr; DUMMER

Good features of HMMER (also in DUMMER):

* It uses position-varying probabilities of not only letters, but also
  starting and extending insertions and deletions.

* It integrates evidence from alternative ways of aligning two
  regions.

In other words, it uses simple probabilities as thoroughly as
possible, to find subtly related regions sensitively.

But, if this is the best way, why isn't it just the (standard) way,
used by all sequence search tools?

[HMMER's theory][] has excessive complexity and minor biases, it's
optimized to answer the wrong question, and its *E*-value conjectures
seem not entirely right ([Frith 2025][frith2025]).  HMMER is optimized
to judge whether a whole sequence contains a match: it would rather
find two matches in two short sequences than three matches in one long
sequence.

DUMMER is optimized to find matches: it doesn't care whether they lie
in long/short same/different sequences.  It's vastly simpler
(comparable to widely-used classic alignment): the hope is to help
thorough probability calculation become "the way" used by other tools
too.

## Current status

* DUMMER uses no heuristic shortcuts: so it's "as sensitive as
  possible", but slow and memory-consuming.

* DUMMER doesn't avoid similarities of unrelated simple sequences,
  like atatatatatatatat, which evolve frequently and independently.

* It can fail due to overflow (numbers getting too big).  This only
  happens when there are very strong similarities.

* The *E*-values are over-estimated when the profile or sequence is
  short.

## Setup

You can get the highest version number from
https://gitlab.com/mcfrith/seq-position-probs/-/tags (or `git clone`
it).  Using the command line, go into the downloaded directory and do
`make`.  That puts the programs in a `bin` subdirectory.

## Usage

DUMMER can compare sequences in FASTA format to profiles in HMMER3/f
format.  You can make profiles with `dummer-build` (new and
lightly-tested).  Or you can (ab)use [HMMER][] profiles, whose
definition doesn't quite fit: you can get DNA profiles from [Dfam][],
or protein profiles from [Pfam][].  Run it like this:

    dummer profiles.hmm sequences.fasta

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
  similar regions around these coordinates.

In more detail, the score comes from adding the probabilities of all
possible alignments passing through the anchor:

score = log<sub>2</sub>[ (sum of alignment probabilities) / (probability of length-0 alignment) ]

DUMMER tries all possible anchors, and outputs all whose *E*-values
are &le; a threshold and are local optima.

The `s` lines show a representative alignment.  This aligns letters
whose probability of being aligned is > 0.5, among all possible
alignments with that anchor.

## Fast, low-memory version

`dummerl` uses half as much memory, and is faster, but is more
prone to numeric overflow.  (It uses single-precision instead of
double-precision floating-point numbers.)

## Options

Show all options and default values:

    dummer --help

Compare each profile to just the forward strand of each DNA sequence:

    dummer -s1 profiles.hmm sequences.fasta

`-s0` means reverse strands only, `-s1` means forward strands only,
and `-s2` means both strands (the default).  This is ignored for
proteins.

Get similarities with *E*-value at most (say) 0.01:

    dummer -e0.01 profiles.hmm sequences.fasta

## Unusual symbols in sequences

Capital/lowercase letters are treated the same.  For nucleotide
sequences, `U` (uracil) is converted to `T` (thymine).  For proteins,
`U` (selenocysteine) is treated the same as `C` (cysteine), and `O`
(pyrrolysine) the same as `K` (lysine).

`dummer` treats other unusual symbols as barriers that break the
sequence into contigs (contiguous sequence).

## Rarely useful features

### Strongest similarities

`-e0` has a special meaning: for each profile versus each contig strand, it
shows the maximum end-anchored, start-anchored, and mid-anchored
scores (in that order).  These scores sum over all alignments ending
at, starting at, or passing through the anchor.  DUMMER gets each kind
of score for all possible anchors, and shows the maximum.

### Random sequences

To calculate *E*-values, it needs to estimate a *K* parameter for each
profile.  To do that, it compares the profile to random sequences.  To
see details of this, give it a profile file only:

    dummer profiles.hmm

For each profile versus each sequence, it shows the maximum
end-anchored, start-anchored, and mid-anchored scores.

For each kind of score, it then shows various estimates of *K* (and
another &lambda; parameter that isn't used currently):

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

These options affect the random sequences:

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

* An *E*-value is: *K*<sub>tot</sub> *N* / 2^score, where
  *K*<sub>tot</sub> is the sum over profiles of *K*, and *N* is the
  sum of sequence lengths.

* &lambda; comes from a not-yet-successful attempt to get better
  *E*-values for short profiles.  Here, the *E*-value for one profile
  is: *KN* / (2^score)^&lambda;.  That would surely work and give us
  better *E*-values, if only we could determine a good value for &lambda;.

## dummer-build

New and experimental.  It makes a profile from an aligned family of
related sequences (in [Stockholm][] format):

    dummer-build alignments.stk > profiles.hmm

It's basically copied from HMMER's hmmbuild, but here are some differences:

* It doesn't assume the alignment is perfect: it optimizes the profile
  with an iterative algorithm (Baum-Welch, aka Expectation
  Maximization).

* End gaps (gaps beyond the start or end of a sequence) are treated as
  missing data (like hmmbuild --fragthresh=1).

* It always does Henikoff position-based sequence weighting: no other
  options.  This down-weights sequences that are similar to each other
  (e.g. human and chimp versions of a sequence).

* Option `--enone` scales the absolute sequence weights so that the
  maximum weight is 1.  That seems good for probabilities fitted to
  the input sequences.

* Often, however, we want probabilities fitted to more-distantly
  related sequences than those in the input.  It does that by
  something similar to HMMER's "entropy weighting", but slightly
  different.  At each position, if the total letter count at that
  position exceeds an "effective sequence number (EFFN)" threshold,
  those counts are downscaled so their total equals the threshold.
  Insertion and deletion counts are treated similarly.

[Dfam]: https://dfam.org/home
[Pfam]: https://www.ebi.ac.uk/interpro/entry/pfam/#table
[frith2025]: https://doi.org/10.1101/2025.03.14.643233
[dmmbuild]: https://github.com/Padraig20/dmmbuild
[edge effects]: https://doi.org/10.1093/nar/29.2.351
[HMMER]: http://hmmer.org
[HMMER's theory]: https://doi.org/10.1371/journal.pcbi.1000069
[MAF]: https://genome.ucsc.edu/FAQ/FAQformat.html#format5
[Stockholm]: https://en.wikipedia.org/wiki/Stockholm_format
