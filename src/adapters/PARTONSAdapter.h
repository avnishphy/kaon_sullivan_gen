#pragma once
#include "../models/GpdModel.h"
#include "../models/PdfModel.h"
#include <string>

// PARTONS adapter skeleton: use PARTONS to compute GPD evolution and Compton/meson production amplitudes.
// Typical usage:
//  - configure PARTONS (database, config xml)
//  - load/create GPD model into PARTONS
//  - call PARTONS service to evolve GPDs and compute observables
class PARTONSAdapter {
public:
  PARTONSAdapter() = default;
  ~PARTONSAdapter() = default;

  // Initialize PARTONS runtime (load modules, DB path)
  void initialize(const std::string &config_file);

  // Given a GPD model, compute forward-limit PDF (H(x,0,0,Q2)) optionally after PARTONS evolution.
  // Return a PdfModel pointer that wraps the output.
  PdfModel* pdf_from_gpd_after_evolution(const GpdModel &gpd, double Q2_target);

  // Compute meson electroproduction amplitude or structure functions directly via PARTONS
  // (more advanced; refer to PARTONS documentation).
};
