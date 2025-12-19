#include <iostream>
#include "core/Generator.h"
#include "models/ToyKaonModel.h"

int main(int argc, char** argv){
    if(argc < 2){
        std::cerr << "Usage: " << argv[0] << " conf.yaml\n";

        return 1;
    }

    std::string cfg = argv[1];
    // simple config parsing (very small) — in real code use yaml-cpp
    // For this starter, we just hardcode reading example values in Generator::from_config

    ToyKaonModel model; //toy model
    Generator gen(cfg, &model);
    gen.initialize();
    gen.run();
    gen.finalize();
    return 0;
}