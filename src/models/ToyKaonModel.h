#pragma once
#include "ModelAPI.h"

class ToyKaonModel : public ModelAPI{
    public:
        ToyKaonModel();
        virtual ~ToyKaonModel();
        double pdf_flavor(int pid, double x, double Q2) override;
        double gpd(const std::string &name, double x, double xi, double t, double Q2) override;
        double model_initial_scale() const override {return 0.5;} //GeV^2 (toy)
        std::string model_name() const override {return "toy_kaon";}

    private:
        double beta_a_, beta_b_, norm_;
};