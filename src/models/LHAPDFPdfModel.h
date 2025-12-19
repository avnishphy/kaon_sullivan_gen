#pragma once
#include "PdfModel.h"
#include <string>
#include <memory>

namespace LHAPDF { class PDF; } // forward declaration from LHAPDF

// Concrete PdfModel that wraps an LHAPDF6 set.
// Example usage:
//   LHAPDFPdfModel model("NNPDF31_nnlo_as_0118", 0);    // member 0
class LHAPDFPdfModel : public PdfModel {
public:
  // setname like "NNPDF31_nnlo_as_0118" or "CT18NNLO/0" (member inlined)
  // member = -1 : use default member 0 created through mkPDF(setname, 0)
  LHAPDFPdfModel(const std::string &setname, int member = 0);
  virtual ~LHAPDFPdfModel();

  double pdf_flavor(int flavor, double x, double Q2) const override;
  double model_initial_scale() const override;
  std::string model_name() const override { return "LHAPDF:" + setname_; }

private:
  std::string setname_;
  int member_;
  LHAPDF::PDF *pdf_; // raw pointer per LHAPDF API; destructor must delete
  double q20_; // initial Q^2 of the grid (if known); 0 if unknown
};
