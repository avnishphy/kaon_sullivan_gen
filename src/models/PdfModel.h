#pragma once
#include <string>

// Abstract interface for PDF-based models.
// Implementations must provide pdf_flavor(flavor, x, Q2) where 'flavor' follows
// PDG/LHAPDF integer codes: d=1,u=2,s=3,c=4,b=5,t=6, g=21, anti-quarks negative (e.g. -2 = anti-u)
// NOTE: Many users prefer using  1..6  for d,u,s,c,b,t (as in LHAPDF). We keep LHAPDF convention.
class PdfModel {
public:
  virtual ~PdfModel() = default;

  // Return f_flavor(x, Q2) (NOT x*f). Use Q2 in GeV^2.
  // Flavor integers should match LHAPDF / PDG conventions.
  virtual double pdf_flavor(int flavor, double x, double Q2) const = 0;

  // Convenience: return total PDF for given flavor-scheme/combination (e.g. valence)
  // Default: call pdf_flavor directly
  virtual double pdf_xf(int flavor, double x, double Q2) const {
    return pdf_flavor(flavor, x, Q2) * x;
  }

  // Some models have a natural initial scale (GeV^2); 0.0 means not specified.
  virtual double model_initial_scale() const { return 0.0; }

  // human-readable name
  virtual std::string model_name() const = 0;
};
