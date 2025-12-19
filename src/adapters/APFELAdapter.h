#pragma once
#include "../models/PdfModel.h"
#include <string>
#include <vector>

// Adapter skeleton for APFEL++ evolution & coefficient functions.
// This is a thin wrapper: given an input "PDF" (either analytic or LHAPDF grids),
// APFEL++ will evolve to target Q2 and compute structure functions at NLO/NNLO.
//
// NOTE: APFEL++ has its own data types and build system. You must link the apfelxx library
// and add #include <apfel/apfelxx.h> or the appropriate header. Below is a *design skeleton*.
class APFELAdapter {
public:
  APFELAdapter() = default;
  ~APFELAdapter() = default;

  // Initialize APFEL engine; set order (LO/NLO/NNLO), alpha_s etc.
  void initialize(int order = 0 /* 0=LO,1=NLO,2=NNLO */);

  // Evolve PDFs from their initial scale to target Q2 and return a PdfModel-like object.
  // Option 1: implement an object that wraps apfel grids and implements PdfModel interface.
  // Option 2: fill an interpolation grid and let your code query it.
  // This function should return a new PdfModel (heap-allocated) implementing PdfModel,
  // or nullptr on failure. Caller responsible for delete.
  PdfModel* evolve_PDF_to_Q2(const PdfModel &inputPDF, double Q2_target);

  // Compute structure functions (F2, FL) at target Q2, using APFEL coefficient convolutions.
  // For performance, this should be implemented inside APFELAdapter using apfel convolution routines.
  double compute_F2_with_APFEL(const PdfModel &pdf, double x, double Q2);
};
