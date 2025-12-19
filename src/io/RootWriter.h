// #pragma once

// #include <string>
// #include "TFile.h"
// #include "TTree.h"

// class RootWriter {
// public:
//     RootWriter(const std::string &fname);
//     ~RootWriter();

//     void initialize();
//     void fill_event(int evtid, double x, double Q2, double t,
//                     double weight, double pdf_k, double gpd);
//     void write();

// private:
//     std::string fname_;
//     TFile *tfile_;
//     TTree *tree_;

//     int evtid_;
//     double x_;
//     double Q2_;
//     double t_;
//     double weight_;
//     double pdf_k_;
//     double gpd_;
// };

#pragma once
// RootWriter.h
// Writes full GenEvent information into a ROOT TTree.

#include <string>

class TFile;
class TTree;

// Forward-declare GenEvent (defined in core/Generator.h)
struct GenEvent;

class RootWriter {
public:
  explicit RootWriter(const std::string &fname);
  ~RootWriter();

  // Initialize ROOT file and TTree + branches
  void initialize();

  // Fill a minimal event (backwards-compatible helper)
  void fill_event(int evtid, double x, double Q2, double t, double weight, double pdf, double gpd);

  // Fill the full event record (preferred)
  void fill_event_full(const GenEvent &ev);

  // Write and close file
  void finalize();

private:
  std::string fname_;
  TFile *tfile_ = nullptr;
  TTree *tree_ = nullptr;

  // Branch variables (mirror GenEvent contents)
  int   evtid_;
  double x_, Q2_, t_, z_lightcone_;
  // four-vectors components
  double k_in_px_, k_in_py_, k_in_pz_, k_in_E_;
  double P_in_px_, P_in_py_, P_in_pz_, P_in_E_;
  double k_out_px_, k_out_py_, k_out_pz_, k_out_E_;
  double L_px_, L_py_, L_pz_, L_E_;
  double kaon_px_, kaon_py_, kaon_pz_, kaon_E_;
  double photon_px_, photon_py_, photon_pz_, photon_E_;

  // physics summary
  double xB_, y_, W2_, s_;
  double xL_, pT_L_;
  double theta_e_, phi_e_, theta_L_, phi_L_, rapidity_L_;
  double t0_, t_calc_;

  // weights
  double weight_;
  int    unweighted_;
  int    is_physics_based_;
  // small diagnostics
  double pdf_diag_, gpd_diag_;
};
