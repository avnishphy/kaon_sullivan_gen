#include "PARTONSAdapter.h"
#include <stdexcept>

// PARTONS is a complete C++ framework. Integration steps:
//
// - Build PARTONS and link your code to it
// - Instantiate a partons::Partons object, configure DB and modules
// - Register your GPD model in PARTONS (or implement a wrapper that queries GpdModel)
// - Use PARTONS services to evolve and compute the forward limit or DVMP amplitudes
//
// See PARTONS docs: https://partons.cea.fr/partons/doc/html/usage.html

void PARTONSAdapter::initialize(const std::string &config_file) {
  (void)config_file;
  throw std::runtime_error("PARTONSAdapter::initialize: implement Partons initialization and DB configuration (see PARTONS docs).");
}

PdfModel* PARTONSAdapter::pdf_from_gpd_after_evolution(const GpdModel &gpd, double Q2_target) {
  (void)gpd; (void)Q2_target;
  throw std::runtime_error("PARTONSAdapter::pdf_from_gpd_after_evolution: must be implemented using PARTONS API.");
}
