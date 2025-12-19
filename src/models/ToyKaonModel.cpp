#include "ToyKaonModel.h"
#include <cmath>


ToyKaonModel::ToyKaonModel() : beta_a_(0.5), beta_b_(1.5) {
    // normalization for x^a (1-x)^b Beta function
    // norm = 1/B(a+1, b+1)
    auto Bfunc = [](double a, double b){
    // Use simple approximation via tgamma
    return std::tgamma(a) * std::tgamma(b) / std::tgamma(a + b);
    };
    norm_ = 1.0 / Bfunc(beta_a_ + 1.0, beta_b_ + 1.0);
}


ToyKaonModel::~ToyKaonModel(){}


double ToyKaonModel::pdf_flavor(int pid, double x, double Q2){
    // pid ignored in this toy implementation; return simple Beta distribution on x
    if(x <= 0.0 || x >= 1.0) return 0.0;
    double val = norm_ * std::pow(x, beta_a_) * std::pow(1.0 - x, beta_b_);
    // no Q2 evolution in toy—just return val
    return val;
}


// Very simple GPD: factorized form H(x,xi,t) = pdf(x) * exp(B t) * (1/(1+xi^2))
double ToyKaonModel::gpd(const std::string &name, double x, double xi, double t, double Q2){
    double pdfx = pdf_flavor(1, x, Q2);
    double B = 1.5; // GeV^-2
    double t_factor = std::exp(B * t);
    double xi_factor = 1.0 / (1.0 + xi*xi);
    return pdfx * t_factor * xi_factor;
}