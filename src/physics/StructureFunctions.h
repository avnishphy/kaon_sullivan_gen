#pragma once
#include "../models/PdfModel.h"

// LO structure functions calculator for spin-1/2 targets (mesons/hadrons).
// This class computes F2 (and optionally FL) in LO from PDFs.
// LO formulas used: F2(x,Q2) = x * sum_f (e_f^2) [ q_f(x,Q2) + qbar_f(x,Q2) ]
// FL at LO is 0 (massless); we provide a simple approximate LO expression (0).
class StructureFunctions {
public:
  // Compute F2 at LO using pdf_model.
  // pdf_model->pdf_flavor(flavor,x,Q2) returns f(x,Q2)
  static double F2_LO(const PdfModel &pdf_model, double x, double Q2);

  // Compute FL (leading-power) — LO massless result is 0; we return 0 here.
  static double FL_LO(const PdfModel &pdf_model, double x, double Q2) { (void)pdf_model; (void)x; (void)Q2; return 0.0; }
};
