#include <iostream>
#include "core/Generator.h"
#include "models/ToyKaonModel.h"

int main(int argc, char** argv){
    // allow optional config path on the command line; otherwise use example
    std::string cfg_path = "examples/kaon_default.yaml";
    if (argc > 1) cfg_path = argv[1];

    // Instantiate toy model
    ToyKaonModel toy_model;

    // create generator with config and toy model
    Generator gen(cfg_path, &toy_model);

    // run workflow
    std::cout << "Initializing generator (cfg: " << cfg_path << ")...\n";
    gen.initialize();

    std::cout << "Running generator...\n";
    gen.run();

    std::cout << "Finalizing generator...\n";
    gen.finalize();

    std::cout << "Test run complete. Output file should be 'events.root' (or the name in you config).\n";
    return 0;
}