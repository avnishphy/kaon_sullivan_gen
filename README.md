# Kaon Sullivan Event Generator (`kaon-sullivan-gen`)

> A modular, extensible Monte Carlo event generator focused on **Sullivan (meson-cloud)** processes with emphasis on **kaon exchange** (p → K + Λ).  
> Designed for EIC-style collisions and fixed-target studies. Produces ROOT event files and is prepared to connect to LHAPDF, APFEL++, PARTONS, HepMC, and PYTHIA.

---

## Table of contents

- [Overview](#overview)  
- [Goals & Scope](#goals--scope)  
- [Repository layout](#repository-layout)  
- [Design principles & architecture](#design-principles--architecture)  
- [Kinematics and physics quantities stored](#kinematics-and-physics-quantities-stored)  
- [Models & adapters](#models--adapters)  
- [YAML steering file (example)](#yaml-steering-file-example)  
- [Build instructions](#build-instructions)  
- [Quick run / examples](#quick-run--examples)  
- [ROOT output: tree & branches](#root-output-tree--branches)  
- [CMake / LHAPDF linking notes (important)](#cmake--lhapdf-linking-notes-important)  
- [Running without physics models](#running-without-physics-models)  
- [Development notes & future work](#development-notes--future-work)  
- [Troubleshooting](#troubleshooting)  
- [Contributing / License](#contributing--license)

---

## Overview

`kaon-sullivan-gen` produces simulated scattering events where the target proton fluctuates into a kaon and a hyperon (Λ), and the hard interaction occurs on the virtual kaon. The generator:

- builds consistent relativistic kinematics for incoming beams and outgoing particles,
- supports electron + proton beams (EIC head-on convention),
- optionally uses PDF/GPD models to compute structure functions and cross sections,
- writes full event records to ROOT files,
- is extensible to LHAPDF / APFEL++ / PARTONS / HepMC / PYTHIA backends.

At this stage the code is a production-grade skeleton: kinematics, IO, and model abstractions are implemented; physics backends are pluggable.

---

## Goals & Scope

- Provide a modular generator capable of producing:
  - **Inclusive Sullivan**: `e p -> e X Λ`
  - **Exclusive (example)**: `e p -> e γ K Λ` (DVCS-like placeholder)
- Store full four-vectors and detailed diagnostics per event
- Allow plugging-in:
  - LHAPDF PDF grids
  - APFEL++ evolution and coefficient functions
  - PARTONS GPD modules (for skewness and GPD evolution)
- Allow running **without** physics models for pure MC / acceptance tests

---

## Repository layout


---

## Design principles & architecture

- **Separation of concerns**
  - `Generator` handles sampling, kinematics, event-loop, IO, unweighting orchestration.
  - `PdfModel` and `GpdModel` are abstract interfaces — concrete implementations (LHAPDF-based or analytic toy) are separate.
  - `StructureFunctions` and `CrossSection` compute physics quantities from `PdfModel`.
  - `APFELAdapter` and `PARTONSAdapter` are where evolution / amplitude code should be called.

- **PDG-based identifiers and four-vectors**
  - All particles are identified with PDG codes and their 4-momenta stored explicitly.

- **Config-driven**
  - YAML steering files determine beam energies, kinematic ranges, output names, model choices, and run options.

- **Unweighted event generation**
  - The generator pre-scans to estimate a maximum weight and uses acceptance-rejection to output unweighted events. The raw `weight` is stored in the ROOT tree for reproducibility and reweighting.

---

## Kinematics & physics quantities stored (per event)

- **Four-vectors** (TLorentzVector): `k_in`, `P_in`, `k_out`, `Lambda`, `kaon`, `photon` (if applicable).
- **DIS variables**: `Q2`, `xB` (Bjorken x), `y`, `W2`, `s`.
- **Sullivan variables**: `xL = (L·k)/(P·k)`, `t = (P − L)^2`, `t0` (kinematic minimum for |t|), `pT_L` (Λ transverse momentum), `theta/phi` angles, `rapidity`.
- **Models & weights**: `weight` (physics weight), `unweighted` (0/1), `is_physics_based` (bool), `pdf_diag`, `gpd_diag`.
- **Diagnostics**: checks for physicality (y in [0,1], positive energies), fallback behaviors.

**Note:** Currently the kaon/Lambda construction uses a simple placeholder (kaon collinear with virtual photon, etc.). Replace with a physics-accurate Sullivan flux and off-shell kinematics when implementing the model.

---

## Models & adapters (current & planned)

### Implemented
- **ToyKaonPdfModel** — analytic beta-function PDF for kaon-like shape (for validation).
- **ToyKaonGpdModel** — factorized GPD H(x,ξ,t) = q(x) e^{B t} / (1 + ξ²).
- **LHAPDFPdfModel** — adapter wrapping LHAPDF6 (returns f(x,Q²) with LHAPDF conventions handled).

### Stubs (skeletons)
- **APFELAdapter** — where a real APFEL++ integration should live (evolution & coefficient convolution).
- **PARTONSAdapter** — where GPD evolution & amplitude computations should be delegated.

---

## YAML steering file (example)

Create or edit a config file like `examples/kaon_no_models.yaml` below:

```yaml
# examples/kaon_no_models.yaml
run:
  n_events: 2000
  seed: 424242
  reaction: inclusive_sullivan
  verbose: true

beam:
  E_e: 10.6     # electron energy (GeV)
  E_p: 100.0    # proton energy (GeV) - collider

kinematics:
  x:
    min: 0.01
    max: 0.5
  Q2:
    min: 1.0
    max: 5.0
  t:
    min: -1.0
    max: -0.01

sullivan:
  z_default: 0.5

unweighting:
  max_weight: 0.0        # 0 -> prescan will estimate
  prescan_events: 500

output:
  directory: output
  filename: events_no_models.root
