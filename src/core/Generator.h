#pragma once
// Generator.h -- extended generator with EIC beams and many kinematic variables.
//
// Requires:
//  - yaml-cpp
//  - ROOT (TLorentzVector, TVector3)
//  - PdfModel / GpdModel interfaces (from models/)
//  - RootWriter extended with fill_event_full(const GenEvent&)

#include <string>
#include <memory>
#include <random>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <yaml-cpp/yaml.h>

#include "../models/PdfModel.h"
#include "../models/GpdModel.h"
#include "../io/RootWriter.h"

namespace PDG {
    constexpr int electron = 11;
    constexpr int positron = -11;
    constexpr int proton = 2212;
    constexpr int neutron = 2112;
    constexpr int kaon_plus = 321;
    constexpr int kaon_minus = -321;
    constexpr int kaon_zero = 310;
    constexpr int lambda0 = 3122;
    constexpr int photon = 22;
}

enum class ReactionTopology {
    Inclusive_Sullivan,
    Exclusive_DVCS_onKaon,
    Exclusive_DVMP_onKaon
};

struct GenEvent {
    unsigned int evtid = 0;

    // core kinematics (inputs)
    double x = 0.0;
    double Q2 = 0.0;
    double t = 0.0;
    double z_lightcone = 0.0; // meson fraction, if used

    // four-vectors
    TLorentzVector k_in;   // incident electron
    TLorentzVector P_in;   // incident proton
    TLorentzVector k_out;  // scattered electron
    TLorentzVector Lambda; // outgoing Lambda (recoil baryon)
    TLorentzVector kaon;   // effective/virtual kaon
    TLorentzVector photon; // for exclusive DVCS-like

    // physics variables
    double xB = 0.0;     // Bjorken x (x_B)
    double y  = 0.0;     // inelasticity
    double W2 = 0.0;     // invariant mass squared of gamma*-N
    double s  = 0.0;     // cms energy squared
    double xL = 0.0;     // (L.k) / (P.k)
    double pT_L = 0.0;   // transverse momentum of Lambda (GeV)
    double theta_e = 0.0;// electron scattering angle (rad)
    double phi_e = 0.0;  // electron azimuth (rad)
    double theta_L = 0.0;// lambda polar angle (rad)
    double phi_L = 0.0;  // lambda azimuth
    double rapidity_L = 0.0;
    double t0 = 0.0;     // kinematic minimum of |t|
    double t_calc = 0.0; // (P - L)^2 computed

    // weights / flags
    double weight = 0.0;         // physics weight (cross section)
    int unweighted = 0;          // 1 if accepted as unweighted event, else 0
    bool is_physics_based = true; // false if no PDF/GPD provided -> pure MC
};

struct GeneratorConfig {
    unsigned int n_events = 1000;
    unsigned int seed = 12345;
    double E_beam_e = 10.0;   // electron beam energy (GeV)
    double E_beam_p = 100.0;  // proton beam energy (GeV)
    std::string output_dir = "./";
    std::string output_root = "events.root";
    ReactionTopology reaction = ReactionTopology::Inclusive_Sullivan;

    // kinematic ranges
    double x_min = 0.01, x_max = 0.8;
    double Q2_min = 1.0, Q2_max = 10.0;
    double t_min = -1.0, t_max = 0.0;

    // Sullivan specifics
    double z_default = 0.5;

    // unweighting
    double max_weight = 0.0;
    unsigned int prescan_events = 1000;

    bool verbose = true;
};

class Generator {
public:
    Generator(const std::string &cfg_file, PdfModel* pdf_model = nullptr, GpdModel* gpd_model = nullptr);
    ~Generator();

    // initialize: parse config and prepare writers and RNG
    void initialize();

    // run generator: produce unweighted events
    void run();

    // finalize: close writers
    void finalize();

    // set/get models
    void set_pdf_model(PdfModel* pdf){ pdf_model_ = pdf; }
    void set_gpd_model(GpdModel* gpd){ gpd_model_ = gpd; }

private:
    // helpers
    void parse_config();
    GenEvent sample_weighted_event(unsigned int id);
    void compute_full_kinematics(GenEvent &ev, double z_lightcone);
    double compute_weight_from_models(const GenEvent &ev);
    void prescan_estimate_max_weight();
    bool unweight_and_write_event(GenEvent &gev);

    // output helpers placeholders
    void write_hepmc(const GenEvent &ev);
    void write_lund(const GenEvent &ev);
    void hand_to_pythia(const GenEvent &ev);

    std::string cfg_file_;
    GeneratorConfig cfg_;

    PdfModel* pdf_model_ = nullptr; // not owned
    GpdModel* gpd_model_ = nullptr; // not owned

    std::unique_ptr<RootWriter> root_writer_;

    // RNG
    std::mt19937 rng_;
    std::uniform_real_distribution<double> dist_x_;
    std::uniform_real_distribution<double> dist_Q2_;
    std::uniform_real_distribution<double> dist_t_;
    std::uniform_real_distribution<double> dist_u01_;

    // internal
    double estimated_max_weight_ = 0.0;
    unsigned long n_generated_ = 0;
    unsigned long n_accepted_ = 0;
};
