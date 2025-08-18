// These gap priors & descriptions are copied from HMMER 3.4, and
// probably not optimal for DUMMER.

struct GapPriors {
  double match, insStart, delStart;
  double insEnd, insExtend;
  double delEnd, delExtend;
};

// For protein sequences, originally trained by Graeme Mitchison in
// the mid-1990's. Notes have been lost, but we believe they were
// trained on an early version of Pfam.
const GapPriors mitchisonAaGapPriors = {0.7939, 0.0278, 0.0135,
					0.1551, 0.1331,
					0.9002, 0.5630};

// For nucleotide sequences, trained on a portion of the rmark dataset.
// Roughly, learned from rmark benchmark - hand-beautified (trimming
// overspecified significant digits).
const GapPriors wheelerNtGapPriors = {2, 0.1, 0.1,
				      0.12, 0.4,
				      0.5, 1};
