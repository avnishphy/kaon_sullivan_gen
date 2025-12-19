// Generator.cpp -- implementation of extended generator
#include "Generator.h"

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>

#include "../physics/StructureFunctions.h"
#include "../physics/CrossSection.h"

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

    if (root["run"] && root["run"]["verbose"]) cfg_.verbose = root["run"]["verbose"].as<bool>();
}

// ----------------- initialize -----------------
void Generator::initialize() {
    parse_config();

    // RNG and sampling distributions
    rng_.seed(cfg_.seed);
    dist_x_ = std::uniform_real_distribution<double>(cfg_.x_min, cfg_.x_max);
    dist_Q2_ = std::uniform_real_distribution<double>(cfg_.Q2_min, cfg_.Q2_max);
    dist_t_ = std::uniform_real_distribution<double>(cfg_.t_min, cfg_.t_max);
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
        ev.t = dist_t_(rng_);
        compute_full_kinematics(ev, cfg_.z_default);
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
    ev.x = dist_x_(rng_);
    ev.Q2 = dist_Q2_(rng_);
    ev.t = dist_t_(rng_);
    // default z (should be sampled from Sullivan flux in realistic version)
    ev.z_lightcone = cfg_.z_default;
    compute_full_kinematics(ev, ev.z_lightcone);

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
void Generator::compute_full_kinematics(GenEvent &ev, double z_lightcone) {
    // masses
    const double m_e = 0.000511;
    const double m_p = 0.938272081;
    const double m_L = 1.115683; // Lambda mass
    const double m_K = 0.493677;

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

    // compute s = (k + P)^2:
    TLorentzVector total = ev.k_in + ev.P_in;
    ev.s = total.M2();

    // compute y using invariant definition y = (P·q)/(P·k)
    // but first build k_out via Q2, x and approximations:
    // We have chosen to sample x and Q2. Use y = Q2 / (x * s)
    if (ev.x > 0 && ev.s > 0) {
        ev.y = ev.Q2 / (ev.x * ev.s);
    } else {
        ev.y = 0.0;
    }

    // scattered electron energy in collider approximate (lab frame):
    // In collider, relate Q2 = -q^2 = 4 E_e E_e' sin^2(theta/2) is more complicated.
    // Use invariant relation: k'·P = (1 - y) k·P
    double kdotP = ev.k_in.Dot(ev.P_in);
    double kprime_dot_P = (1.0 - ev.y) * kdotP;
    // Solve for k' energy in lab (E') approximately by assuming direction same as incoming electron (small angle)
    // But better: compute E' from invariants: Q2 = - (k - k')^2 = 2 k·k' (neglect m_e)
    // So, k·k' = Q2/2. We know k·k' = E_e*E_out - |p_e||p_out| cos theta ~ E_e*E_out (for small angles).
    // We'll approximate E_out using y: E_out ≈ (1 - y) * E_e (collinear approx)
    double E_out = (1.0 - ev.y) * E_e;
    if (E_out <= m_e) E_out = m_e + 1e-6;

    // scattering angle theta from Q2 relation: Q2 = 4 E_e E_out sin^2(theta/2)
    double arg = ev.Q2 / (4.0 * E_e * E_out);
    if (arg < 0.0) arg = 0.0;
    if (arg > 1.0) arg = 1.0;
    ev.theta_e = 2.0 * std::asin(std::sqrt(arg));
    ev.phi_e = 0.0;

    // k' 4-vector
    double p_out = std::sqrt(std::max(0.0, E_out*E_out - m_e*m_e));
    double px = p_out * std::sin(ev.theta_e) * std::cos(ev.phi_e);
    double py = p_out * std::sin(ev.theta_e) * std::sin(ev.phi_e);
    double pz = p_out * std::cos(ev.theta_e);
    ev.k_out.SetPxPyPzE(px, py, pz, E_out);

    // virtual photon q = k - k'
    TLorentzVector q = ev.k_in - ev.k_out;
    // Q2 from kinematics check
    double Q2_check = -q.M2();
    (void)Q2_check; // optional check

    // W2 = (P + q)^2
    TLorentzVector Pplusq = ev.P_in + q;
    ev.W2 = Pplusq.M2();

    // Kaon and Lambda kinematics (placeholder)
    // We put the kaon collinear with q direction with energy z * E_p (very approximate)
    ev.z_lightcone = z_lightcone;
    double E_kaon = std::max(m_K + 1e-6, z_lightcone * E_p);
    TVector3 qdir = q.Vect().Unit();
    double pka = std::sqrt(std::max(0.0, E_kaon*E_kaon - m_K*m_K));
    ev.kaon.SetPxPyPzE(qdir.X()*pka, qdir.Y()*pka, qdir.Z()*pka, E_kaon);

    // Lambda recoiling: approximate from momentum conservation: Lambda = P + q - kaon
    TLorentzVector proton4(0.0, 0.0, pz_p, E_p);
    ev.Lambda = proton4 + q - ev.kaon;

    // compute x_L = (L·k)/(P·k)
    double kdotP_val = ev.k_in.Dot(ev.P_in);
    if (kdotP_val != 0.0) {
        ev.xL = ev.Lambda.Dot(ev.k_in) / kdotP_val;
    } else {
        ev.xL = 0.0;
    }

    // compute transverse momentum of Lambda relative to beam (pT)
    ev.pT_L = std::sqrt(ev.Lambda.Px()*ev.Lambda.Px() + ev.Lambda.Py()*ev.Lambda.Py());

    // compute t = (P - L)^2
    ev.t_calc = (ev.P_in - ev.Lambda).M2();
    ev.t = ev.t_calc; // store computed t (overwrites sampled t as consistency)

    // compute theta, phi, rapidity of Lambda
    TVector3 pL = ev.Lambda.Vect();
    ev.theta_L = pL.Theta();
    ev.phi_L = pL.Phi();
    ev.rapidity_L = 0.5 * std::log((ev.Lambda.E() + pL.Z()) / (ev.Lambda.E() - pL.Z() + 1e-12));

    // compute t0: t0 = ((1-xL)/xL) * (m_L^2 - xL * m_p^2)
    if (ev.xL > 0.0) {
        ev.t0 = ((1.0 - ev.xL) / ev.xL) * (m_L*m_L - ev.xL * m_p*m_p);
    } else {
        ev.t0 = 0.0;
    }

    // store xB from invariants: xB = Q2 / (2 P·q)
    double Pdotq = ev.P_in.Dot(q);
    if (Pdotq > 0.0) {
        ev.xB = ev.Q2 / (2.0 * Pdotq);
    } else {
        ev.xB = 0.0;
    }

    // y alternative definition
    if (kdotP_val > 0.0) ev.y = Pdotq / kdotP_val;

    // safety: if some values are NaN or inf, clamp them
    if (!std::isfinite(ev.xB)) ev.xB = 0.0;
    if (!std::isfinite(ev.y)) ev.y = 0.0;
    if (!std::isfinite(ev.W2)) ev.W2 = 0.0;
    if (!std::isfinite(ev.pT_L)) ev.pT_L = 0.0;
}

// ----------------- compute weight from models -----------------
double Generator::compute_weight_from_models(const GenEvent &ev) {
    // If no models, 0.0 (should not be called if checked earlier)
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

    // optional: APFEL evolution call could go here to evolve pdf_view to ev.Q2
    // For now compute LO F2 and then LO cross section
    double F2 = StructureFunctions::F2_LO(*pdf_view, ev.x, ev.Q2);
    double FL = StructureFunctions::FL_LO(*pdf_view, ev.x, ev.Q2);
    double d2 = CrossSection::d2sigma_dx_dQ2_LO(*pdf_view, ev.x, ev.Q2, ev.s);
    // Note: must include Jacobian if sampling not uniform in true measure.

    return d2;
}

// ----------------- unweight & write -----------------
bool Generator::unweight_and_write_event(GenEvent &gev) {
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
        // You must implement that method in RootWriter so that the entire GenEvent is stored.
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

        // If pure MC (no physics) weight is 1.0; do simple acceptance (always accept)
        if (!gev.is_physics_based) {
            gev.unweighted = 1;
            root_writer_->fill_event_full(gev);
            ++n_accepted_;
            continue;
        }

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
