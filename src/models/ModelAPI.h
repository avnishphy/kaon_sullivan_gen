#pragma once
#include <string>

class ModelAPI{
    public:
        virtual ~ModelAPI(){}
        // flavor id: user-defined (eg., 1 -> kaon valence). Return pdf(x, Q2)
        virtual double pdf_flavor(int pid, double x, double Q2) = 0;
        // basic GPD: name like "H"; xi, t, Q2
        virtual double gpd(const std::string &name, double x, double xi, double t, double Q2) = 0;
        virtual double model_initial_scale() const = 0;
        virtual std::string model_name() const = 0;
};