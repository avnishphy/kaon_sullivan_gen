#include "StructureFunctions.h"
#include <array>

double StructureFunctions::F2_LO(const PdfModel &pdf_model, double x, double Q2) {
  // quark charges squared (u,d,s,c,b,t) in PDG indexed order flavor->charge2
  // LHAPDF flavor codes: d=1, u=2, s=3, c=4, b=5, t=6
  static const std::array<double, 7> qcharge2 = {0.0,  // 0 unused
    1.0/9.0,   // d (1)
    4.0/9.0,   // u (2)
    1.0/9.0,   // s (3)
    4.0/9.0,   // c (4) (same as u)
    1.0/9.0,   // b (5)
    4.0/9.0    // t (6) (same as u)
  };

  double sum = 0.0;
  // For each quark flavor include q + qbar
  for (int flav = 1; flav <= 6; ++flav) {
    double q = pdf_model.pdf_flavor(flav, x, Q2);
    double qbar = pdf_model.pdf_flavor(-flav, x, Q2);
    double e2 = qcharge2[flav];
    sum += e2 * (q + qbar);
  }
  double F2 = x * sum;
  return F2;
}
