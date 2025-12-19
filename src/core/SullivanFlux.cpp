#include "SullivanFlux.h"
#include <cmath>

namespace SullivanFlux{
    double meson_flux(double t, double z){
        // Very simple toy flux: exponential t-dependence times a z-shape.
        double B = 3.0; // GeV^2 slope (toy)
        double z0 = 0.5;
        double fz = std::exp(-((z - z0)*(z - z0))/0.02);
        double ft = std::exp(B*t);
        return fz * ft;
    }
}