#include "Generator.h"
#include "SullivanFlux.h"
#include "io/RootWriter.h"
#include <random>
#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>

Generator::Generator(const std::string &config_file, ModelAPI* model)
: cfg_file_(config_file), model_(model) {}

Generator::~Generator(){}

void Generator::initialize()
{
    YAML::Node cfg = YAML::LoadFile(cfg_file_);

    if (!cfg["run"] || !cfg["kinematics"] || !cfg["output"]) {
    throw std::runtime_error("Invalid YAML: missing top-level sections");
    }

    // --- Run control ---
    n_events_ = cfg["run"]["n_events"].as<unsigned>();
    seed_     = cfg["run"]["seed"].as<unsigned>();

    // --- Kinematics ---
    x_min_  = cfg["kinematics"]["x"]["min"].as<double>();
    x_max_  = cfg["kinematics"]["x"]["max"].as<double>();

    if (x_min_ >= x_max_)
        throw std::runtime_error("Invalid x range");

    Q2_min_ = cfg["kinematics"]["Q2"]["min"].as<double>();
    Q2_max_ = cfg["kinematics"]["Q2"]["max"].as<double>();

    if (Q2_min_ <= 0)
        throw std::runtime_error("Q2 must be positive");

    t_min_  = cfg["kinematics"]["t"]["min"].as<double>();
    t_max_  = cfg["kinematics"]["t"]["max"].as<double>();

    // --- Output ---
    std::string out_dir  = cfg["output"]["directory"].as<std::string>();
    std::string out_name = cfg["output"]["filename"].as<std::string>();
    out_file_ = out_dir + "/" + out_name;

    writer_ = new RootWriter(out_file_);
    writer_->initialize();

    std::cout << "[Generator] Initialized\n"
              << "  Events : " << n_events_ << "\n"
              << "  Output : " << out_file_ << "\n";
}


void Generator::run(){
    std::mt19937 rng(seed_);
    std::uniform_real_distribution<double> ux(x_min_, x_max_);
    std::uniform_real_distribution<double> uQ2(Q2_min_, Q2_max_);
    std::uniform_real_distribution<double> ut(t_min_, t_max_);

    for(unsigned int i=0;i<n_events_;++i){
        double x = ux(rng);
        double Q2 = uQ2(rng);
        double t = ut(rng);

        // simple toy flow: compute kaon PDF at x, Q2
        double pdf_k = model_->pdf_flavor(1, x, Q2); // '1' = toy kaon valence
        double gpd = model_->gpd("H", x, 0.0, t, Q2);

        double flux = SullivanFlux::meson_flux(t, 0.5); // toy flux depends on t and z

        double weight  = flux * pdf_k * gpd;

        // Build a minimal event record: x, Q2, t, weight, pdf, gpd
        writer_->fill_event(i, x, Q2, t, weight, pdf_k, gpd);
    }
}

void Generator::finalize(){
    if(writer_) writer_->write();
    std::cout << "Generation complete.\n";
}