#include "core/Generator.h"
#include <iostream>

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: kaon_sullivan_gen <config.yaml>\n";
        return 1;
    }

    std::string cfg = argv[1];

    Generator gen(cfg);  // no PDF/GPD → pure MC mode
    gen.initialize();
    gen.run();
    gen.finalize();

    return 0;
}
