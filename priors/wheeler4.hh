// This Dirichlet mixture & description were copied from HMMER 3.4
// It doesn't satisfy Chargaff parity (A=T, C=G): undesirable for DNA?

// Match emission priors are trained on Rmark3 database
// Xref: ~wheelert/notebook/2011/0325_nhmmer_new_parameters

const double wheeler4[] = {
  0.24,
  0.16,  0.45,  0.12,   0.39,

  0.26,
  0.09,  0.03,  0.09,   0.04,

  0.08,
  1.29,  0.40,  6.58,   0.51,

  0.42,
  1.74,  1.49,  1.57,   1.95,
};
