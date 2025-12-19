# kaon-sullivan-gen (starter)


Minimal skeleton for a Kaon-Sullivan event generator.


## What you get
- A tiny C++ generator that samples kinematics, evaluates a toy kaon PDF and GPD, computes a toy Sullivan flux, and writes events into a ROOT TTree.
- Clean modular layout to add APFEL++, PARTONS, LHAPDF, HepMC3, Pythia8 etc.


## Requirements
- CMake (>= 3.12)
- A C++17 compiler
- ROOT (for the minimal writer) — optional for building if you adapt the IO
- (Optional) GoogleTest if you want to run tests


## Build
```bash
mkdir build && cd build
cmake ..
make -j