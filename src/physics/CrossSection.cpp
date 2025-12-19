// implementing a very basic xsec just as the template workflow for now

#include "CrossSection.h"
#include <cmath>

static constexpr double alpha_em = 1.0 / 137.035999; // fine-structure constant ~1/137

double CrossSection::d2sigma_dx_dQ2_LO(const PdfModel &pdf_model, double x, double Q2, double s) {
  if (x <= 0.0 || x >= 1.0) return 0.0;
  if (Q2 <= 0.0) return 0.0;
  if (s <= 0.0) return 0.0;

  double y = Q2 / (x * s);
  if (y <= 0.0 || y >= 1.0) return 0.0;

  double F2 = StructureFunctions::F2_LO(pdf_model, x, Q2);
  double FL = StructureFunctions::FL_LO(pdf_model, x, Q2);

  // One-photon exchange formula (unpolarized):
  // d^2σ/(dx dQ^2) = (4π α^2) / (x Q^4) * [ (1 - y + y*y/2) * F2 - (y*y/2) * FL ]
  double prefactor = 4.0 * M_PI * alpha_em * alpha_em / (x * std::pow(Q2, 2));
  double kinematic = (1.0 - y + 0.5 * y * y) * F2 - 0.5 * y * y * FL;

  return prefactor * kinematic;
}
