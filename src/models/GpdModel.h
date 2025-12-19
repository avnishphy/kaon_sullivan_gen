#pragma once
#include <string>
#include <stdexcept>
#include "PdfModel.h"

// Abstract GPD interface: returns GPD H/E/... depending on `name`.
// gpd(name, x, xi, t, Q2)
// - x in [0,1]
// - xi skewness
// - t in GeV^2 (negative for spacelike: e.g. -0.2)
// - Q2 in GeV^2
class GpdModel {
public:
  virtual ~GpdModel() = default;

  // Return GPD value (no builtin prefactor). Forward limit for PDFs often is H(x,0,0,Q2).
  virtual double gpd(const std::string &name, double x, double xi, double t, double Q2) const = 0;

  // Default: forward-limit PDF from GPD = H(x, xi->0, t->0, Q2).
  // Some GPDs need a Jacobian or 2x factor; check your model. We expose this helper to get a PDF.
  virtual double pdf_from_gpd(double x, double Q2) const {
    // by default use H forward limit:
    return gpd("H", x, 0.0, 0.0, Q2);
  }

  virtual double model_initial_scale() const { return 0.0; }
  virtual std::string model_name() const = 0;
};
