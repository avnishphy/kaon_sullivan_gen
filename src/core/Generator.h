#pragma once
#include <string>
#include "models/ModelAPI.h"

class Generator{
    public:
        Generator(const std::string &config_file, ModelAPI* model);
        ~Generator();
        void initialize();
        void run();
        void finalize();
    private:
        std::string cfg_file_;
        ModelAPI* model_;
        unsigned int n_events_;
        unsigned int seed_;
        std::string out_file;
        double x_min_, x_max_, Q2_min_, Q2_max_, t_min_, t_max_;
        // writer (ROOT)
        class RootWriter* writer_ = nullptr;
};