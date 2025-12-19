// #include "RootWriter.h"
// #include <iostream>

// RootWriter::RootWriter(const std::string &fname)
//     : fname_(fname), tfile_(nullptr), tree_(nullptr)
// {}

// RootWriter::~RootWriter()
// {
//     if(tfile_) tfile_->Close();
// }

// void RootWriter::initialize()
// {
//     tfile_ = TFile::Open(fname_.c_str(), "RECREATE");
//     if (!tfile_ || tfile_->IsZombie()) {
//         std::cerr << "ERROR: Could not open ROOT file " << fname_ << "\n";
//         return;
//     }

//     tree_ = new TTree("events", "kaon sullivan events");

//     tree_->Branch("evtid",    &evtid_,   "evtid/I");
//     tree_->Branch("x",        &x_,       "x/D");
//     tree_->Branch("Q2",       &Q2_,      "Q2/D");
//     tree_->Branch("t",        &t_,       "t/D");
//     tree_->Branch("weight",   &weight_,  "weight/D");
//     tree_->Branch("pdf_k",    &pdf_k_,   "pdf_k/D");
//     tree_->Branch("gpd",      &gpd_,     "gpd/D");
// }

// void RootWriter::fill_event(int evtid, double x, double Q2, double t,
//                             double weight, double pdf_k, double gpd)
// {
//     evtid_ = evtid;
//     x_ = x;
//     Q2_ = Q2;
//     t_ = t;
//     weight_ = weight;
//     pdf_k_ = pdf_k;
//     gpd_ = gpd;

//     tree_->Fill();
// }

// void RootWriter::write()
// {
//     tfile_->Write();
// }

// RootWriter.cpp
#include "RootWriter.h"
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <iomanip>

// Include GenEvent so we can access fields for fill_event_full
#include "../core/Generator.h"

RootWriter::RootWriter(const std::string &fname) : fname_(fname) {}
RootWriter::~RootWriter() {
  // ensure finalize called (safe if already finalized)
  finalize();
}

void RootWriter::initialize() {
  tfile_ = TFile::Open(fname_.c_str(), "RECREATE");
  if (!tfile_ || tfile_->IsZombie()) {
    std::cerr << "[RootWriter] ERROR: Could not open ROOT file: " << fname_ << std::endl;
    tfile_ = nullptr;
    return;
  }
  tree_ = new TTree("events", "Generated events");

  // Minimal/diagnostic branches
  tree_->Branch("evtid", &evtid_, "evtid/I");
  tree_->Branch("x", &x_, "x/D");
  tree_->Branch("Q2", &Q2_, "Q2/D");
  tree_->Branch("t", &t_, "t/D");
  tree_->Branch("z_lightcone", &z_lightcone_, "z_lightcone/D");

  // four-vectors: incoming electron
  tree_->Branch("k_in_px", &k_in_px_, "k_in_px/D");
  tree_->Branch("k_in_py", &k_in_py_, "k_in_py/D");
  tree_->Branch("k_in_pz", &k_in_pz_, "k_in_pz/D");
  tree_->Branch("k_in_E",  &k_in_E_,  "k_in_E/D");

  // proton
  tree_->Branch("P_in_px", &P_in_px_, "P_in_px/D");
  tree_->Branch("P_in_py", &P_in_py_, "P_in_py/D");
  tree_->Branch("P_in_pz", &P_in_pz_, "P_in_pz/D");
  tree_->Branch("P_in_E",  &P_in_E_,  "P_in_E/D");

  // scattered electron
  tree_->Branch("k_out_px", &k_out_px_, "k_out_px/D");
  tree_->Branch("k_out_py", &k_out_py_, "k_out_py/D");
  tree_->Branch("k_out_pz", &k_out_pz_, "k_out_pz/D");
  tree_->Branch("k_out_E",  &k_out_E_,  "k_out_E/D");

  // Lambda
  tree_->Branch("Lambda_px", &L_px_, "Lambda_px/D");
  tree_->Branch("Lambda_py", &L_py_, "Lambda_py/D");
  tree_->Branch("Lambda_pz", &L_pz_, "Lambda_pz/D");
  tree_->Branch("Lambda_E",  &L_E_,  "Lambda_E/D");

  // Kaon (effective)
  tree_->Branch("kaon_px", &kaon_px_, "kaon_px/D");
  tree_->Branch("kaon_py", &kaon_py_, "kaon_py/D");
  tree_->Branch("kaon_pz", &kaon_pz_, "kaon_pz/D");
  tree_->Branch("kaon_E",  &kaon_E_,  "kaon_E/D");

  // Photon (exclusive)
  tree_->Branch("photon_px", &photon_px_, "photon_px/D");
  tree_->Branch("photon_py", &photon_py_, "photon_py/D");
  tree_->Branch("photon_pz", &photon_pz_, "photon_pz/D");
  tree_->Branch("photon_E",  &photon_E_,  "photon_E/D");

  // physics vars
  tree_->Branch("xB", &xB_, "xB/D");
  tree_->Branch("y",  &y_,  "y/D");
  tree_->Branch("W2", &W2_, "W2/D");
  tree_->Branch("s",  &s_,  "s/D");
  tree_->Branch("xL", &xL_, "xL/D");
  tree_->Branch("pT_L", &pT_L_, "pT_L/D");

  tree_->Branch("theta_e", &theta_e_, "theta_e/D");
  tree_->Branch("phi_e", &phi_e_, "phi_e/D");
  tree_->Branch("theta_L", &theta_L_, "theta_L/D");
  tree_->Branch("phi_L", &phi_L_, "phi_L/D");
  tree_->Branch("rapidity_L", &rapidity_L_, "rapidity_L/D");

  tree_->Branch("t0", &t0_, "t0/D");
  tree_->Branch("t_calc", &t_calc_, "t_calc/D");

  // weights & flags
  tree_->Branch("weight", &weight_, "weight/D");
  tree_->Branch("unweighted", &unweighted_, "unweighted/I");
  tree_->Branch("is_physics_based", &is_physics_based_, "is_physics_based/I");

  // diagnostics
  tree_->Branch("pdf_diag", &pdf_diag_, "pdf_diag/D");
  tree_->Branch("gpd_diag", &gpd_diag_, "gpd_diag/D");
}

void RootWriter::fill_event(int evtid, double x, double Q2, double t, double weight, double pdf, double gpd) {
  if (!tree_) return;
  evtid_ = evtid;
  x_ = x; Q2_ = Q2; t_ = t;
  weight_ = weight;
  pdf_diag_ = pdf;
  gpd_diag_ = gpd;
  unweighted_ = 0;
  is_physics_based_ = 0;
  // clear four-vectors
  k_in_px_ = k_in_py_ = k_in_pz_ = k_in_E_ = 0.0;
  P_in_px_ = P_in_py_ = P_in_pz_ = P_in_E_ = 0.0;
  k_out_px_ = k_out_py_ = k_out_pz_ = k_out_E_ = 0.0;
  L_px_ = L_py_ = L_pz_ = L_E_ = 0.0;
  kaon_px_ = kaon_py_ = kaon_pz_ = kaon_E_ = 0.0;
  photon_px_ = photon_py_ = photon_pz_ = photon_E_ = 0.0;
  xB_ = y_ = W2_ = s_ = xL_ = pT_L_ = 0.0;
  theta_e_ = phi_e_ = theta_L_ = phi_L_ = rapidity_L_ = 0.0;
  t0_ = t_calc_ = 0.0;
  tree_->Fill();
}

void RootWriter::fill_event_full(const GenEvent &ev) {
  if (!tree_) return;
  evtid_ = static_cast<int>(ev.evtid);
  x_ = ev.x; Q2_ = ev.Q2; t_ = ev.t; z_lightcone_ = ev.z_lightcone;

  // incoming electron
  k_in_px_ = ev.k_in.Px();
  k_in_py_ = ev.k_in.Py();
  k_in_pz_ = ev.k_in.Pz();
  k_in_E_  = ev.k_in.E();

  // proton
  P_in_px_ = ev.P_in.Px();
  P_in_py_ = ev.P_in.Py();
  P_in_pz_ = ev.P_in.Pz();
  P_in_E_  = ev.P_in.E();

  // scattered electron
  k_out_px_ = ev.k_out.Px();
  k_out_py_ = ev.k_out.Py();
  k_out_pz_ = ev.k_out.Pz();
  k_out_E_  = ev.k_out.E();

  // Lambda
  L_px_ = ev.Lambda.Px();
  L_py_ = ev.Lambda.Py();
  L_pz_ = ev.Lambda.Pz();
  L_E_  = ev.Lambda.E();

  // kaon
  kaon_px_ = ev.kaon.Px();
  kaon_py_ = ev.kaon.Py();
  kaon_pz_ = ev.kaon.Pz();
  kaon_E_  = ev.kaon.E();

  // photon
  photon_px_ = ev.photon.Px();
  photon_py_ = ev.photon.Py();
  photon_pz_ = ev.photon.Pz();
  photon_E_  = ev.photon.E();

  // physics
  xB_ = ev.xB; y_ = ev.y; W2_ = ev.W2; s_ = ev.s;
  xL_ = ev.xL; pT_L_ = ev.pT_L;
  theta_e_ = ev.theta_e; phi_e_ = ev.phi_e;
  theta_L_ = ev.theta_L; phi_L_ = ev.phi_L;
  rapidity_L_ = ev.rapidity_L;
  t0_ = ev.t0; t_calc_ = ev.t_calc;

  weight_ = ev.weight;
  unweighted_ = ev.unweighted;
  is_physics_based_ = ev.is_physics_based ? 1 : 0;

  // diagnostics: attempt to store model diagnostics if any
  pdf_diag_ = 0.0;
  gpd_diag_ = 0.0;

  tree_->Fill();
}

void RootWriter::finalize() {
  if (!tfile_) return;
  tfile_->cd();
  if (tree_) tree_->Write();
  tfile_->Close();
  delete tfile_;
  tfile_ = nullptr;
  tree_ = nullptr;
  std::cout << "[RootWriter] Wrote ROOT file: " << fname_ << std::endl;
}
