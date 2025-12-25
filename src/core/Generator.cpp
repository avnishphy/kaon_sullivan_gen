// Generator.cpp -- implementation of extended generator
#include "Generator.h"

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <random>
#include <limits>
#include <algorithm>

#include "../physics/StructureFunctions.h"
#include "../physics/CrossSection.h"

// portable PI (compatible with C++17)
namespace {
    constexpr double PI = 3.14159265358979323846;
}

// ----------------- Constructor / Destructor -----------------
Generator::Generator(const std::string &cfg_file, PdfModel* pdf_model, GpdModel* gpd_model)
 : cfg_file_(cfg_file), pdf_model_(pdf_model), gpd_model_(gpd_model)
{
    // default seed; will be overridden in initialize()
    rng_.seed(12345);
    dist_u01_ = std::uniform_real_distribution<double>(0.0, 1.0);
}

Generator::~Generator() {
    finalize();
}

// ----------------- parse_config -----------------
void Generator::parse_config() {
    YAML::Node root = YAML::LoadFile(cfg_file_);

    if (root["run"]) {
        const YAML::Node &r = root["run"];
        if (r["n_events"]) cfg_.n_events = r["n_events"].as<unsigned>();
        if (r["seed"]) cfg_.seed = r["seed"].as<unsigned>();
        if (r["reaction"]) {
            std::string s = r["reaction"].as<std::string>();
            if (s == "inclusive_sullivan") cfg_.reaction = ReactionTopology::Inclusive_Sullivan;
            else if (s == "exclusive_dvcs") cfg_.reaction = ReactionTopology::Exclusive_DVCS_onKaon;
            else if (s == "exclusive_dvmp") cfg_.reaction = ReactionTopology::Exclusive_DVMP_onKaon;
        }
    }

    if (root["beam"]) {
        const YAML::Node &b = root["beam"];
        if (b["E_e"]) cfg_.E_beam_e = b["E_e"].as<double>();
        if (b["E_p"]) cfg_.E_beam_p = b["E_p"].as<double>();
    }

    if (root["kinematics"]) {
        const YAML::Node &k = root["kinematics"];
        if (k["x"]) {
            if (k["x"]["min"]) cfg_.x_min = k["x"]["min"].as<double>();
            if (k["x"]["max"]) cfg_.x_max = k["x"]["max"].as<double>();
        }
        if (k["Q2"]) {
            if (k["Q2"]["min"]) cfg_.Q2_min = k["Q2"]["min"].as<double>();
            if (k["Q2"]["max"]) cfg_.Q2_max = k["Q2"]["max"].as<double>();
        }
        // t not used for sampling any more but keep values if present in YAML
        if (k["t"]) {
            if (k["t"]["min"]) cfg_.t_min = k["t"]["min"].as<double>();
            if (k["t"]["max"]) cfg_.t_max = k["t"]["max"].as<double>();
        }
    }

    if (root["output"]) {
        const YAML::Node &o = root["output"];
        if (o["directory"]) cfg_.output_dir = o["directory"].as<std::string>();
        if (o["filename"]) cfg_.output_root = o["filename"].as<std::string>();
    }

    if (root["sullivan"]) {
        const YAML::Node &s = root["sullivan"];
        if (s["z_default"]) cfg_.z_default = s["z_default"].as<double>();
    }

    if (root["unweighting"]) {
        const YAML::Node &u = root["unweighting"];
        if (u["max_weight"]) cfg_.max_weight = u["max_weight"].as<double>();
        if (u["prescan_events"]) cfg_.prescan_events = u["prescan_events"].as<unsigned>();
    }

    // New optional parameters controlling Lambda pT sampling (kept optional)
    if (root["Lambda"]) {
        const YAML::Node &L = root["Lambda"];
        if (L["pT0"]) cfg_.Lambda_pT0 = L["pT0"].as<double>();            // default 0.2
        if (L["pT_max"]) cfg_.Lambda_pT_max = L["pT_max"].as<double>();    // default 1.5
    }

    if (root["run"] && root["run"]["verbose"]) cfg_.verbose = root["run"]["verbose"].as<bool>();
}

// ----------------- initialize -----------------
void Generator::initialize() {
    parse_config();

    // set defaults for new Lambda pT params if not provided
    if (!(cfg_.Lambda_pT0 > 0.0)) cfg_.Lambda_pT0 = 0.20;     // GeV
    if (!(cfg_.Lambda_pT_max > 0.0)) cfg_.Lambda_pT_max = 1.50; // GeV

    // RNG and sampling distributions
    rng_.seed(cfg_.seed);
    dist_x_ = std::uniform_real_distribution<double>(cfg_.x_min, cfg_.x_max);
    dist_Q2_ = std::uniform_real_distribution<double>(cfg_.Q2_min, cfg_.Q2_max);
    // do not sample t up-front; t is computed
    dist_u01_ = std::uniform_real_distribution<double>(0.0, 1.0);

    // prepare ROOT writer
    std::string fullpath = cfg_.output_dir;
    if (!fullpath.empty() && fullpath.back() != '/') fullpath += '/';
    fullpath += cfg_.output_root;
    root_writer_ = std::make_unique<RootWriter>(fullpath);
    root_writer_->initialize();

    // warn if no physics model provided
    if (!pdf_model_ && !gpd_model_) {
        std::cerr << "[Generator] WARNING: No PDF or GPD model provided. Events will be PURE MC (no physics weights).\n";
    } else {
        if (cfg_.verbose) {
            if (pdf_model_) std::cout << "[Generator] Using PDF model: " << pdf_model_->model_name() << "\n";
            if (gpd_model_) std::cout << "[Generator] Using GPD model: " << gpd_model_->model_name() << "\n";
        }
    }

    // estimate max weight unless given
    if (cfg_.max_weight <= 0.0) {
        prescan_estimate_max_weight();
        if (cfg_.verbose) std::cout << "[Generator] Estimated max weight = " << estimated_max_weight_ << "\n";
    } else {
        estimated_max_weight_ = cfg_.max_weight;
    }

    n_generated_ = 0;
    n_accepted_ = 0;

    if (cfg_.verbose) {
        std::cout << "[Generator] Initialization complete. Output: " << fullpath << "\n";
        std::cout << "  E_e = " << cfg_.E_beam_e << " GeV, E_p = " << cfg_.E_beam_p << " GeV\n";
        std::cout << "  Lambda pT sampling: pT0 = " << cfg_.Lambda_pT0 << " GeV, pT_max = " << cfg_.Lambda_pT_max << " GeV\n";
    }
}

// ----------------- prescan -----------------
void Generator::prescan_estimate_max_weight() {
    double maxw = 0.0;
    unsigned int nscan = std::min<unsigned int>(cfg_.prescan_events, 20000);
    for (unsigned int i=0; i<nscan; ++i) {
        GenEvent ev;
        ev.x = dist_x_(rng_);
        ev.Q2 = dist_Q2_(rng_);

        // sample electron direction isotropically (uniform cos theta, uniform phi)
        double u = dist_u01_(rng_);
        double cos_th = 2.0*u - 1.0;
        ev.theta_e = std::acos(std::clamp(cos_th, -1.0, 1.0));
        ev.phi_e = 2.0 * PI * dist_u01_(rng_);

        // compute kinematics (this will sample Lambda with pT-based proposal)
        ev.z_lightcone = cfg_.z_default;
        compute_full_kinematics(ev, ev.z_lightcone);

        double w = compute_weight_from_models(ev);
        if (w > maxw) maxw = w;
    }
    estimated_max_weight_ = maxw * 1.2 + 1e-20;
    if (estimated_max_weight_ <= 0.0) estimated_max_weight_ = 1e-12;
}

// ----------------- sample one weighted event -----------------
GenEvent Generator::sample_weighted_event(unsigned int id) {
    GenEvent ev;
    ev.evtid = id;

    const unsigned int max_attempts = 10000;
    unsigned int attempt = 0;
    bool success = false;

    while (attempt < max_attempts && !success) {
        ++attempt;

        // sample primary kinematics
        ev.x = dist_x_(rng_);
        ev.Q2 = dist_Q2_(rng_);

        // ----------------------------
        // 1) scattered electron direction: isotropic sampling (spherical symmetry)
        //    sample cos(theta) uniformly in [-1,1] and phi uniformly in [0,2pi)
        // ----------------------------
        double u_th = dist_u01_(rng_);
        double cos_th_e = 2.0 * u_th - 1.0;            // uniform in [-1,1]
        ev.theta_e = std::acos(std::clamp(cos_th_e, -1.0, 1.0));
        ev.phi_e = 2.0 * PI * dist_u01_(rng_); // uniform phi

        // compute kinematics (this will set k_in, P_in, sample Lambda with pT-based proposal and derive K = P - L)
        ev.z_lightcone = cfg_.z_default; // placeholder; we'll derive xL from L later
        compute_full_kinematics(ev, ev.z_lightcone);

        // check that kaon has positive energy and finite values (otherwise resample)
        if (!std::isfinite(ev.kaon.E()) || ev.kaon.E() <= 1e-9) {
            continue; // resample
        }

        // enforce |t| < 0.9 already in compute_full_kinematics? if not, check here
        if (!std::isfinite(ev.t) || std::abs(ev.t) >= 0.9) {
            continue; // resample
        }

        // enforce that the inclusive X (ep->eXL) has a positive mass
        if (!std::isfinite(ev.X_out.E()) || ev.X_out.E() <= 1e-9 || ev.X_out.M2() <= 1e-9){
            continue; // resample
        }

        // All good
        success = true;
    }

    if (!success) {
        // failed to find a configuration
        ev.is_physics_based = (pdf_model_ || gpd_model_);
        ev.weight = 0.0;
        ev.unweighted = 0;
        return ev;
    }

    // compute the physics weight (or placeholder if no model)
    if (!pdf_model_ && !gpd_model_) {
        ev.is_physics_based = false;
        ev.weight = 1.0; // pure MC weight (placeholder)
    } else {
        ev.is_physics_based = true;
        ev.weight = compute_weight_from_models(ev);
    }
    ev.unweighted = 0;
    return ev;
}

// ----------------- compute_full_kinematics -----------------
// Proton-side sub-process is P -> K + L with L generated forward-peaked in pT (exponential).
// Kaon is derived as K = P - L and may be off-shell. The electron scatters off that kaon (eK->eX) for weights.
void Generator::compute_full_kinematics(GenEvent &ev, double /*z_lightcone_unused*/) {
    // masses (GeV)
    const double m_e = 0.000511;
    const double m_p = 0.938272081;
    const double m_L = 1.115683; // Lambda mass (on-shell)
    const double m_K = 0.493677; // Kaon nominal mass (we allow off-shell K)

    // Incoming beams: head-on collider convention:
    // Electron moves +z, Proton moves -z (opposite directions).
    double E_e = cfg_.E_beam_e;
    double E_p = cfg_.E_beam_p;

    // electron 4-vector k_in
    double pz_e = std::sqrt(std::max(0.0, E_e*E_e - m_e*m_e));
    ev.k_in.SetPxPyPzE(0.0, 0.0, +pz_e, E_e);

    // proton 4-vector P_in (opposite direction)
    double pz_p = -std::sqrt(std::max(0.0, E_p*E_p - m_p*m_p));
    ev.P_in.SetPxPyPzE(0.0, 0.0, pz_p, E_p);

    // compute s = (k + P)^2 (legacy)
    TLorentzVector total = ev.k_in + ev.P_in;
    ev.s = total.M2();

    // ----------------------------
    // Outgoing electron kinematics from sampled angles and Q2
    // Using invariant relation (massless approximation where appropriate):
    //   Q^2 = 2 k·k' ≈ 2 E_e E_out (1 - cos theta)
    // => E_out = Q^2 / (2 E_e (1 - cos theta))
    // Keep safe clamps when cos theta ~ 1.
    // ----------------------------
    double cos_th = std::cos(ev.theta_e);
    double denom = 2.0 * E_e * (1.0 - cos_th);
    double E_out = 0.0;
    if (denom > 1e-12) {
        E_out = ev.Q2 / denom;
    } else {
        // extremely forward scattering: clamp to near-beam energy
        E_out = std::max(m_e + 1e-6, 0.99 * E_e);
    }
    if (!std::isfinite(E_out) || E_out <= m_e) E_out = m_e + 1e-6;

    // set outgoing electron momentum vector from sampled direction
    double p_out = std::sqrt(std::max(0.0, E_out*E_out - m_e*m_e));
    double px = p_out * std::sin(ev.theta_e) * std::cos(ev.phi_e);
    double py = p_out * std::sin(ev.theta_e) * std::sin(ev.phi_e);
    double pz = p_out * std::cos(ev.theta_e);
    ev.k_out.SetPxPyPzE(px, py, pz, E_out);

    // virtual photon q = k - k'
    TLorentzVector q = ev.k_in - ev.k_out;

    // store Q2-check (optional)
    double Q2_check = -q.M2();
    (void)Q2_check;

    // compute P·q and y if possible
    double Pdotq = ev.P_in.Dot(q);
    double kdotP_val = ev.k_in.Dot(ev.P_in);
    if (kdotP_val > 0.0) {
        ev.y = Pdotq / kdotP_val;
    } else {
        if (ev.x > 0.0 && ev.s > 0.0) ev.y = ev.Q2 / (ev.x * ev.s);
        else ev.y = 0.0;
    }

    // ----------------------------
    // Proton-side: sample Lambda (L) with forward-peaked pT distribution (truncated exponential)
    // Then derive kaon K = P - L (kaon may be off-shell).
    //
    // Strategy:
    //  - sample pT from truncated exponential with scale pT0 up to pT_max
    //  - sample phi uniform in [0,2pi)
    //  - for given pT compute pz_max allowed so that E_L <= E_p - eps (so E_K > 0)
    //  - sample pz uniformly in [-pz_max, pz_max]
    // ----------------------------
    double pT0 = cfg_.Lambda_pT0;
    double pT_max = cfg_.Lambda_pT_max;
    if (!(pT0 > 0.0)) pT0 = 0.20;
    if (!(pT_max > 0.0)) pT_max = 1.50;

    // truncated exponential inverse CDF:
    // CDF(p) = (1 - exp(-p/p0)) / (1 - exp(-pTmax/p0))
    double exp_cut = std::exp(-pT_max / pT0);
    double u_pT = dist_u01_(rng_);
    double one_minus_exp_cut = 1.0 - exp_cut;
    if (one_minus_exp_cut <= 0.0) one_minus_exp_cut = 1e-12;
    double pT = -pT0 * std::log(1.0 - u_pT * one_minus_exp_cut);

    // phi for Lambda
    ev.phi_L = 2.0 * PI * dist_u01_(rng_);

    // compute maximum allowed |pz| given E_L <= E_p - eps
    double eps = 1e-6;
    double E_L_allowed_max = std::max(m_L, (E_p > m_L + eps ? E_p - eps : m_L));
    double max_pz2 = std::max(0.0, E_L_allowed_max*E_L_allowed_max - m_L*m_L - pT*pT);
    double pz_max = std::sqrt(max_pz2);

    // if pz_max is zero, then pT is too large; clamp and set pz=0
    double pz_L = 0.0;
    if (pz_max > 0.0) {
        // sample pz uniformly in [-pz_max, pz_max]
        double u_pz = dist_u01_(rng_);
        pz_L = (2.0 * u_pz - 1.0) * pz_max;
    } else {
        pz_L = 0.0;
    }

    // build Lambda 4-vector
    double pxL = pT * std::cos(ev.phi_L);
    double pyL = pT * std::sin(ev.phi_L);
    double pL2 = pT*pT + pz_L*pz_L;
    double E_L = std::sqrt(std::max(0.0, m_L*m_L + pL2));
    ev.Lambda.SetPxPyPzE(pxL, pyL, pz_L, E_L);

    // Now derive Kaon as K = P - L (may be off-shell)
    ev.kaon = ev.P_in - ev.Lambda;

    // compute xL = (L·k)/(P·k)
    if (kdotP_val != 0.0) {
        ev.xL = ev.Lambda.Dot(ev.k_in) / kdotP_val;
    } else {
        ev.xL = 0.0;
    }

    // compute transverse momentum of Lambda relative to beam (pT)
    ev.pT_L = std::sqrt(ev.Lambda.Px()*ev.Lambda.Px() + ev.Lambda.Py()*ev.Lambda.Py());

    // compute t = (P - L)^2
    ev.t_calc = (ev.P_in - ev.Lambda).M2();
    ev.t = ev.t_calc; // computed t

    // if |t| too large, we'll rely on caller to resample; still compute angles
    TVector3 pL_vec = ev.Lambda.Vect();
    ev.theta_L = pL_vec.Theta();
    ev.rapidity_L = 0.5 * std::log((ev.Lambda.E() + pL_vec.Z()) / (ev.Lambda.E() - pL_vec.Z() + 1e-12));

    // compute xB from invariants: xB = Q2 / (2 P·q) if P·q > 0
    if (Pdotq > 0.0) {
        ev.xB = ev.Q2 / (2.0 * Pdotq);
        ev.x = ev.xB;
    } else {
        ev.xB = 0.0;
    }

    // recompute y from invariants if possible
    if (kdotP_val > 0.0) ev.y = Pdotq / kdotP_val;

    // safety: clamp NaNs / infs
    if (!std::isfinite(ev.xB)) ev.xB = 0.0;
    if (!std::isfinite(ev.y)) ev.y = 0.0;
    if (!std::isfinite(ev.pT_L)) ev.pT_L = 0.0;

    // compute t0: t0 = ((1-xL)/xL) * (m_L^2 - xL * m_p^2)
    if (ev.xL > 0.0) {
        ev.t0 = ((1.0 - ev.xL) / ev.xL) * (m_L*m_L - ev.xL * m_p*m_p);
    } else {
        ev.t0 = 0.0;
    }

    // X = k_in + P_in - k_out - Lambda (ep->eXL)
    ev.X_out = ev.k_in + ev.P_in - ev.k_out - ev.Lambda;

    // W2 = (P_in + k_in - k_out)**2 = (Lambda + X_out)^2
    ev.W2 = (ev.Lambda + ev.X_out).M2();

}

// ----------------- compute_weight from models -----------------
double Generator::compute_weight_from_models(const GenEvent &ev) {
    // If no models, 1.0 (pure MC)
    if (!pdf_model_ && !gpd_model_) return 1.0;

    // Wrap GPD as PDF if necessary
    PdfModel *pdf_view = nullptr;
    std::unique_ptr<PdfModel> temp_pdf;

    if (pdf_model_) {
        pdf_view = pdf_model_;
    } else if (gpd_model_) {
        // minimal wrapper: forward limit H->pdf (toy)
        struct GpdToPdfWrapper : public PdfModel {
            const GpdModel* g_;
            GpdToPdfWrapper(const GpdModel* g): g_(g) {}
            double pdf_flavor(int flavor, double x, double Q2) const override {
                (void)flavor;
                return g_->pdf_from_gpd(x, Q2);
            }
            std::string model_name() const override { return std::string("gpd_forward_") + g_->model_name(); }
        };
        temp_pdf = std::make_unique<GpdToPdfWrapper>(gpd_model_);
        pdf_view = temp_pdf.get();
    }

    // Use electron-kaon invariant s_eK = (k + K)^2 for eK subprocess weight
    TLorentzVector kplusK = ev.k_in + ev.kaon;
    double s_eK = kplusK.M2();

    // optional: APFEL evolution call could go here to evolve pdf_view to ev.Q2
    // For now compute LO F2 and then LO cross section using s_eK (so eK -> eX)
    double F2 = StructureFunctions::F2_LO(*pdf_view, ev.x, ev.Q2);
    double FL = StructureFunctions::FL_LO(*pdf_view, ev.x, ev.Q2);

    // If s_eK is not physically reasonable, guard it
    if (!std::isfinite(s_eK) || s_eK <= 0.0) {
        // fallback to proton-level s if needed
        s_eK = ev.s > 0.0 ? ev.s : 1.0;
    }

    double d2 = CrossSection::d2sigma_dx_dQ2_LO(*pdf_view, ev.x, ev.Q2, s_eK);
    // Note: must include Jacobian if sampling not uniform in true measure.

    return d2;
}

// ----------------- unweight & write -----------------
// Only store/write events that pass |t| < 0.9
bool Generator::unweight_and_write_event(GenEvent &gev) {
    // Apply t cut before writing: if event fails cut, treat as rejected
    if (!std::isfinite(gev.t)) return false;
    if (std::abs(gev.t) >= 0.9) {
        // do not store this event
        return false;
    }

    double u = dist_u01_(rng_);
    double thresh = gev.weight / estimated_max_weight_;
    if (thresh > 1.0) {
        // update estimate conservatively
        estimated_max_weight_ = gev.weight * 1.2;
        thresh = gev.weight / estimated_max_weight_;
    }
    if (u <= thresh) {
        gev.unweighted = 1;
        // write ROOT extended event -> requires RootWriter::fill_event_full(const GenEvent&)
        root_writer_->fill_event_full(gev);

        // write other outputs (placeholders)
        write_hepmc(gev);
        write_lund(gev);
        // hand_to_pythia(gev); // optional if Pythia linked

        ++n_accepted_;
        return true;
    } else {
        return false;
    }
}

// ----------------- run -----------------
void Generator::run() {
    if (cfg_.verbose) {
        std::cout << "[Generator] Running " << cfg_.n_events << " target unweighted events...\n";
    }

    unsigned int produced_id = 0;
    while (n_accepted_ < cfg_.n_events) {
        GenEvent gev = sample_weighted_event(produced_id);
        ++n_generated_;
        ++produced_id;

        if (gev.weight <= 0.0 && gev.is_physics_based) {
            if (cfg_.verbose) std::cout << "[Generator] Event " << gev.evtid << " has zero weight (physics) -> skipping\n";
            continue;
        }

        // If pure MC (no physics) weight is 1.0; accept only if it passes t cut
        if (!gev.is_physics_based) {
            if (std::isfinite(gev.t) && std::abs(gev.t) < 0.9) {
                gev.unweighted = 1;
                root_writer_->fill_event_full(gev);
                ++n_accepted_;
            }
            continue;
        }

        // physics-based: unweight and write only if |t| < 0.9 (handled inside unweight_and_write_event)
        unweight_and_write_event(gev);

        if (cfg_.verbose && (n_generated_ % 1000 == 0)) {
            std::cout << "[Generator] attempted: " << n_generated_ << ", accepted: " << n_accepted_ << "\n";
        }

        // safe guard to avoid infinite loop if estimated_max_weight very small
        if (n_generated_ > cfg_.n_events * 1000 && cfg_.verbose) {
            std::cerr << "[Generator] WARNING: many generation attempts without enough accepted events. Check max_weight and sampling.\n";
        }
    }

    if (cfg_.verbose) {
        std::cout << "[Generator] Finished. Attempts: " << n_generated_ << ", accepted: " << n_accepted_ << "\n";
    }
}

// ----------------- placeholders for other outputs -----------------
void Generator::write_hepmc(const GenEvent &ev) {
    (void)ev;
    // TODO: implement HepMC writing
}
void Generator::write_lund(const GenEvent &ev) {
    (void)ev;
    // TODO: implement LUND writer
}
void Generator::hand_to_pythia(const GenEvent &ev) {
    (void)ev;
    // TODO: implement Pythia handing
}

// ----------------- finalize -----------------
void Generator::finalize() {
    if (root_writer_) {
        // I expect RootWriter::finalize() to write and close the file
        root_writer_->finalize();
        root_writer_.reset();
    }
}
