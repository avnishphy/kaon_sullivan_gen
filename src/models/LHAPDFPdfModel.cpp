#include "LHAPDFPdfModel.h"
#include <LHAPDF/LHAPDF.h>   // LHAPDF6 header
#include <cmath>
#include <iostream>

LHAPDFPdfModel::LHAPDFPdfModel(const std::string &setname, int member)
  : setname_(setname), member_(member), pdf_(nullptr), q20_(0.0)
{
  // Two common ways to instantiate:
  // 1) mkPDF("SETNAME", member)  e.g. mkPDF("CT18NNLO/0") or mkPDF("CT18NNLO",0)
  // 2) Make a PDFSet and mkPDFs
  // We'll try mkPDF(setname, member) which is flexible.
  try {
    pdf_ = LHAPDF::mkPDF(setname_, member_);
  } catch (const std::exception &e) {
    std::cerr << "LHAPDFPdfModel: error creating PDF '" << setname_
              << "' member " << member_ << ": " << e.what() << "\n";
    pdf_ = nullptr;
    return;
  }

  // LHAPDF works in Q (GeV), but grids are defined in Q^2 usually. We can't always extract Q0.
  // If you need the set's nominal Q0, inspect pdf_->q2min()/q2max() if provided (some versions provide).
  // For portability, we leave q20_ = 0.0 (unknown).
}

LHAPDFPdfModel::~LHAPDFPdfModel(){
  if (pdf_) delete pdf_;
}

double LHAPDFPdfModel::model_initial_scale() const {
  return q20_;
}

double LHAPDFPdfModel::pdf_flavor(int flavor, double x, double Q2) const {
  if (!pdf_) return 0.0;
  if (x <= 0.0 || x >= 1.0) return 0.0;
  if (Q2 <= 0.0) {
    // If user passed Q2 but it's zero -> use Q = 1.0 GeV (conservative)
    Q2 = 1.0;
  }
  double Q = std::sqrt(Q2); // LHAPDF::PDF::xfxQ takes Q, not Q^2
  // LHAPDF returns x * f(flavor, x, Q), so we divide by x to return f(x,Q2).
  // Some implementations also provide xfxQ2; the portable approach is xfxQ.
  double xfx = pdf_->xfxQ(flavor, x, Q);
  double f = xfx / x;
  return f;
}
