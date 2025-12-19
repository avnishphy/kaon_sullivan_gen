#pragma once
#include "StructureFunctions.h"
#include "../models/PdfModel.h"

// Compute double-differential DIS cross section (unpolarized) d^2σ/dx dQ2
// For electron scattering off a hadron (lepton mass neglected).
// Uses one-photon exchange LO formula:
// d^2σ/dx dQ2 = (4πα^2 / (x Q^4)) * [ 1 - y + y^2/2 ] * F2(x,Q2)  (neglecting FL term for LO)
//
// You must provide s = 2 E_e * E_h (in lab frame) or the hadron+electron invariant s.
// y is computed as Q^2 / (x s).
class CrossSection {
public:
  // Compute d^2σ/dx dQ2 using LO F2 (barn units via alpha constant); returns value in pb/GeV^2? We keep natural units: GeV^-2.
  // It's up to the caller to convert units. Use alpha=1/137.035999.
  static double d2sigma_dx_dQ2_LO(const PdfModel &pdf_model, double x, double Q2, double s);
};
