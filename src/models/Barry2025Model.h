#pragma once
#include "ModelAPI.h"

// https://arxiv.org/abs/2510.11979v1
class Barry2025Model : public ModelAPI{
    public:
        Barry2025Model();
        virtual ~Barry2025Model();
        double pdf_flavor(int pid, double x, double Q2) override;
        double gpd(const std::string &name, double x, double xi, double t, double Q2) override;
        double model_initial_scale() const override {return 0.5;} //GeV^2 (toy)
        std::string model_name() const override {return "Barry2025";}

    private:
        double beta_a_, beta_b_, norm_;
};