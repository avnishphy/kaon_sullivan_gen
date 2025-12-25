# Generator README

> Detailed description of the extended generator implementation (file: `Generator.cpp`) — how sampling is done, what is sampled, what physics cuts are applied, and everything a reader needs to reproduce/understand the generator's behaviour.

---

## Table of contents

* [Overview](#overview)
* [High-level flow](#high-level-flow)
* [Configuration (YAML) and defaults](#configuration-yaml-and-defaults)
* [What variables are sampled and how](#what-variables-are-sampled-and-how)
* [Kinematics construction (details inside `compute_full_kinematics`)](#kinematics-construction-details-inside-compute_full_kinematics)
* [Physics weight calculation](#physics-weight-calculation)
* [Unweighting procedure & prescan](#unweighting-procedure--prescan)
* [Physics / analysis cuts currently applied](#physics--analysis-cuts-currently-applied)
* [Output formats and hooks](#output-formats-and-hooks)
* [Safety checks, numerical clamps and resampling rules](#safety-checks-numerical-clamps-and-resampling-rules)
* [Limitations, caveats, and recommended tuning](#limitations-caveats-and-recommended-tuning)
* [Minimal usage example (C++) and example `config.yml`][#minimal-usage-example-c-and-example-configyml)
* [Author / contact / TODOs](#author--contact--todos)

---

## Overview

This generator produces events for processes where the proton side is treated as a two-body system `P -> K + L` (proton → kaon + Lambda) and the electron scatters off the resulting `K` (kaon) — e.g. inclusive Sullivan-like configurations or related exclusive topologies. The generator supports *physics-weighted* events using a `PdfModel` or a `GpdModel` if provided, else it falls back to pure MC events (weight `= 1`).

All major logic is in:

* sampling: `sample_weighted_event()`
* kinematics: `compute_full_kinematics()`
* weights: `compute_weight_from_models()`
* unweighting + writing: `unweight_and_write_event()` and `run()`

---

## High-level flow

1. Parse configuration file (YAML) — see config keys below.
2. Initialize RNG and sampling distributions for `x` and `Q²`.
3. If `max_weight` not provided by configuration, run a *prescan* to estimate `max_weight`.
4. Loop: sample weighted events (`sample_weighted_event`) until `n_events` unweighted events have been accepted.
5. For each sampled event:

   * build electron and proton-side kinematics (`compute_full_kinematics`)
   * compute physics weight (if a model is available)
   * apply physics cuts (e.g. `|t| < 0.9`) and unweight (accept with probability `weight/estimated_max_weight`)
   * on acceptance, write event to output (ROOT writer) and optional hooks.

---

## Configuration (YAML) and defaults

The generator reads options from a YAML file. Important keys (and default behaviors) are:

```yaml
run:
  n_events:      # required target number of UNWEIGHTED events to produce
  seed:          # RNG seed (unsigned)
  reaction:      # "inclusive_sullivan" | "exclusive_dvcs" | "exclusive_dvmp"
  verbose: true/false

beam:
  E_e:  # electron beam energy (GeV)
  E_p:  # proton beam energy (GeV)

kinematics:
  x:
    min:  # lower bound for uniform sampling of x
    max:  # upper bound for uniform sampling of x
  Q2:
    min:  # lower bound for uniform sampling of Q2 (GeV^2)
    max:  # upper bound for uniform sampling of Q2 (GeV^2)
  t:     # (kept for compatibility but t is NOT sampled up-front anymore)

output:
  directory:   # output directory
  filename:    # output root filename

sullivan:
  z_default:  # placeholder light-cone fraction used in some parts (double)

unweighting:
  max_weight:     # user-provided max weight (optional). If <=0, prescan is used.
  prescan_events: # how many events to sample for prescan (optional)

Lambda:
  pT0:    # scale parameter for Lambda pT exponential (default: 0.20 GeV)
  pT_max: # truncation for Lambda pT (default: 1.50 GeV)
```

**Defaults**

* `Lambda.pT0` default = `0.20 GeV` if not provided.
* `Lambda.pT_max` default = `1.50 GeV`.
* If `unweighting.max_weight` missing or `<= 0`, `prescan_estimate_max_weight()` runs (up to 20,000 events or `prescan_events`, whichever is smaller) and sets `estimated_max_weight = max_found * 1.2 + 1e-20`.
* RNG seed default in code is `12345` but overwritten by `run.seed`.

---

## What variables are sampled and how

The generator samples the following **primary** variables:

1. **`x`** (Bjorken-x / partonic x)

   * Distribution: **uniform** in `[x_min, x_max]` using `std::uniform_real_distribution`.

2. **`Q²`**

   * Distribution: **uniform** in `[Q2_min, Q2_max]` using `std::uniform_real_distribution`.

3. **Outgoing electron direction**

   * **θ_e (theta)**: sampled *isotropically* by sampling `u ~ U(0,1)` then `cos(θ) = 2u - 1` (uniform in `[-1,1]`) and `θ = arccos(cosθ)`.

     * Thus the electron direction is uniform on the sphere (`d cosθ` uniform).
   * **φ_e (phi)**: uniform in `[0, 2π)` via `2π * U(0,1)`.

4. **Lambda (proton-side) transverse momentum pT**

   * Model: **truncated exponential (forward-peaked)** distribution for `pT`.
   * Parameters: `pT0` (scale) and `pT_max` (truncation).
   * Sampling method: inverse CDF of truncated exponential:

     * For `u ~ U(0,1)`, `pT = -pT0 * ln(1 - u * (1 - exp(-pT_max / pT0)))`.
   * `phi_L` (Lambda azimuth): uniform in `[0,2π)`.

5. **Lambda longitudinal momentum pz**

   * For a sampled `pT`, compute `pz_max` from the requirement `E_L <= E_p - eps` (so that kaon energy `E_K = E_p - E_L` remains positive).
   * `pz` is sampled **uniformly** in `[-pz_max, pz_max]` (if `pz_max > 0`).
   * If `pz_max` is zero (pT too large), `pz` is set to zero.

Derived quantities (not directly sampled, but computed):

* `E_out` (outgoing electron energy) derived from `(Q^2, θ_e)` using the approximation:

  ```
  Q^2 = 2 E_e E_out (1 - cosθ)  => E_out = Q^2 / (2 E_e (1 - cosθ))
  ```

  with care for near-forward (`cosθ ≈ 1`) cases (clamped to ≤ `0.99 * E_e` as fallback).
* Kaon four-vector `K = P_in - Lambda` (may be off-shell).
* `t` computed as `(P - L)^2`.
* `xL` computed as `(L·k_in) / (P·k_in)`.
* `xB` recomputed from invariants if `P·q > 0`: `xB = Q² / (2 P·q)` and then `ev.x` is set to `xB`.

**Random generators used**

* `std::mt19937` seeded with `cfg_.seed`.
* `std::uniform_real_distribution` for `x`, `Q²`, `u01`, etc.

**Resampling limit**
If the event generation for a single weighted event fails to produce acceptable kinematics, the attempt loop allows up to `max_attempts = 10000` before returning a dummy zero-weight event.

---

## Kinematics construction (details inside `compute_full_kinematics`)

Important steps and formulas implemented:

* **Masses** (hard-coded in the function):

  * `m_e = 0.000511 GeV` (electron)
  * `m_p = 0.938272081 GeV` (proton)
  * `m_L = 1.115683 GeV` (Lambda)
  * `m_K = 0.493677 GeV` (kaon nominal mass — kaon can be off-shell in this approach)

* **Beam 4-vectors**: head-on collider convention where electron momentum is `+z` and proton is `-z`.

  * electron input 4-vector `k_in` built from `E_e` and `pz_e = sqrt(E_e^2 - m_e^2)`.
  * proton input `P_in` built from `E_p` and `pz_p = -sqrt(E_p^2 - m_p^2)`.

* **Outgoing electron**:

  * `E_out` derived from `Q²` and `θ_e` using:

    ```
    denom = 2 * E_e * (1 - cosθ)
    E_out = Q² / denom   (if denom > 1e-12)
    ```

    * If denominator is too small (very forward scattering), `E_out` is clamped to `max(m_e + 1e-6, 0.99*E_e)`.

* **Virtual photon** `q = k_in - k_out` and checks on `Q²` via `-q.M2()`.

* **Lambda generation** (proton-side):

  * `pT` sampled from truncated exponential (see previous section).
  * `phi_L` uniform in `[0, 2π)`.
  * `pz_L` uniform in `[-pz_max, pz_max]` where `pz_max` ensures `E_L <= E_p - eps`.
  * Construct Lambda 4-vector `Lambda` with `E_L = sqrt(m_L^2 + pT^2 + pz_L^2)`.
  * Kaon 4-vector `kaon = P_in - Lambda` (can be off-shell).
  * `xL = (Lambda·k_in) / (P·k_in)`.
  * `t = (P - Lambda)^2` stored as `t_calc` and `t`.

* **Derived invariants**

  * `P·q` used to compute `y = P·q / (k·P)` and `xB = Q² / (2 P·q)` (when `P·q > 0`).
  * `X_out = k_in + P_in - k_out - Lambda` → used to check that the inclusive `X` has a positive mass `M²`.
  * `W² = (Lambda + X_out)²`.

* **t0**:

  ```
  t0 = ((1 - xL) / xL) * (m_L^2 - xL * m_p^2)   (if xL > 0)
  ```

---

## Physics weight calculation

* If neither `PdfModel* pdf_model_` nor `GpdModel* gpd_model_` is provided, the generator sets `weight = 1.0` (pure MC).
* If a `GpdModel` is provided but no `PdfModel`, a minimal wrapper `GpdToPdfWrapper` is constructed to provide a toy forward limit mapping `pdf(x,Q2) ← gpd.pdf_from_gpd(x,Q2)`. This is **not** a rigorous GPD → PDF reduction, it's a toy wrapper used to compute structure functions.
* The weight calculation uses:

  1. `StructureFunctions::F2_LO(*pdf_view, x, Q2)` and `StructureFunctions::FL_LO(*pdf_view, x, Q2)` to compute structure functions at LO.
  2. `CrossSection::d2sigma_dx_dQ2_LO(*pdf_view, x, Q2, s_eK)` is called to compute a LO double-differential cross section for the **electron–kaon** subprocess. Here:

     * `s_eK = (k_in + kaon)^2` — if `s_eK` is non-physical, it falls back to proton-level `s` (`ev.s`).
* **Important note**: The code contains a comment that a Jacobian might be needed depending on the sampling measure. *No explicit Jacobian is applied in the code.* That means the sampling is uniform in `(x, Q2, cosθ, φ)` with the weight `d2sigma_dx_dQ2_LO` used directly as the event weight — the user must be aware of the sampling measure when interpreting weights or differential cross-sections from the generated events.

---

## Unweighting procedure & prescan

**Prescan**

* If `cfg_.max_weight <= 0.0`, `prescan_estimate_max_weight()` runs with `nscan = min(prescan_events, 20000)`.
* It samples `nscan` weighted events and records the maximum found weight `maxw`.
* `estimated_max_weight_ = maxw * 1.2 + 1e-20` (safety margin).
* If `estimated_max_weight_ <= 0`, it's replaced by `1e-12` to avoid division by zero.

**Unweighting**

* For each physics-based sampled event with weight `w`, the acceptance probability is:

  ```
  P_accept = w / estimated_max_weight_
  ```
* A random `u ~ U(0,1)` is drawn and the event is accepted if `u <= P_accept`.
* If `P_accept > 1.0` (i.e. `w > estimated_max_weight_`), the code **conservatively updates**:

  ```
  estimated_max_weight_ = w * 1.2
  ```

  and recomputes `thresh = w / estimated_max_weight_` (now < 1.0).
* Only accepted events that pass the `|t| < 0.9` cut are written to the output.

**Target**

* The main `run()` method continues sampling until `n_accepted_ == cfg_.n_events` unweighted events have been written.

---

## Physics / analysis cuts currently applied

The generator enforces the following **hard-coded** checks / cuts:

1. **`|t| < 0.9`**

   * This is applied in two places:

     * `sample_weighted_event()` rejects (resamples) intermediate configurations where `|t| >= 0.9`.
     * `unweight_and_write_event()` also enforces `|t| < 0.9` and will not write events failing this criterion.
   * This is the primary physical acceptance cut present.

2. **Finite and positive energy checks**

   * After kinematics, the code verifies:

     * `ev.kaon.E()` must be finite and `> 1e-9`.
     * `ev.X_out.E()` must be finite and `> 1e-9`.
     * `ev.X_out.M2()` must be `> 1e-9` (so the inclusive `X` has positive mass).
   * If any of these fail, the event is resampled.

3. **Numerical checks**

   * Ensure `ev.t` is finite; otherwise resample.
   * `pz_max` is computed to ensure `E_L <= E_p - eps` so that kaon energy `E_K` is positive.

4. **Attempt limit**

   * `max_attempts = 10000` attempts per weighted event; if no acceptable kinematics found, the returned event has `weight=0` and will be skipped downstream.

---

## Output formats and hooks

* **Primary output**: A ROOT file written via `root_writer_` (constructed from `output.directory` and `output.filename`). The generator calls:

  ```cpp
  root_writer_->fill_event_full(gev);
  ```

  and at finalize: `root_writer_->finalize()`.

* **Other (placeholders)**

  * `write_hepmc(gev)` — placeholder (TODO) for HepMC writing.
  * `write_lund(gev)` — placeholder (TODO) for LUND format writing.
  * `hand_to_pythia(gev)` — placeholder (TODO) for handing over to Pythia event generator. Currently not implemented (functions are empty stubs).

---

## Safety checks, numerical clamps and resampling rules

The code includes several pragmatic protections:

* Clamp `cos(theta)` values to `[-1, 1]`.
* Fallback clamping for `E_out` for forward scattering.
* `one_minus_exp_cut` guard for truncated exponential (`>= 1e-12`).
* Fallback `estimated_max_weight` lower bound of `1e-12`.
* `eps = 1e-6` used when computing `E_L_allowed_max` to avoid `E_K <= 0`.
* All computed double fields are checked for `std::isfinite()` and for small positive thresholds (`> 1e-9`) before accepting.
* If too many attempts occur without enough accepted events, a warning is printed (helps detect poor sampling or badly estimated `max_weight`).

---

## Limitations, caveats, and recommended tuning

**Known limitations**

* Kaon may be **off-shell** (constructed as `K = P - Lambda`) — the physics consequences must be considered when analyzing results.
* `GpdModel` → `PdfModel` conversion is a **toy wrapper**; it is not a rigorous GPD→PDF mapping. If using GPDs you should provide a proper PDF interface or implement a more realistic forward limit.
* No Jacobian is applied for the transformation between the generator sampling measure and the true differential measure. If you want cross sections per physical phase-space measure, you must account for the sampling Jacobian.
* `t` is **not** sampled directly; it is derived from the sampled Lambda (pT/pz) distribution — so acceptance in `t` is determined implicitly and may be inefficient depending on the `pT0`/`pT_max` choices.
* Output formats (HepMC/LUND/Pythia) are placeholders — not implemented.

**Recommended tuning**

* If acceptance is low or unweighting is inefficient, set `unweighting.max_weight` to a known safe upper bound or increase `unweighting.prescan_events` to get a better `estimated_max_weight`.
* Tune `Lambda.pT0` and `Lambda.pT_max` to reflect realistic forward-peaked distributions (current defaults are `0.20 GeV` and `1.50 GeV`).
* Provide a realistic `PdfModel` (preferred) or implement a correct GPD→PDF mapping if using `GpdModel`.
* Consider sampling `Q²` or `x` from non-uniform distributions (e.g., 1/x or 1/Q²) if those reflect physics in your analysis and apply the appropriate Jacobian when computing weights.

---

## Minimal usage example (C++) and example `config.yml`

### Minimal C++ (pseudo) usage

```cpp
// create models (or pass nullptr to run pure-MC)
PdfModel* pdf = create_my_pdf_model(...);   // or nullptr
GpdModel* gpd = nullptr;

Generator gen("config.yml", pdf, gpd);
gen.initialize();
gen.run();
gen.finalize();
```

### Example `config.yml`

```yaml
run:
  n_events: 5000
  seed: 123456
  reaction: inclusive_sullivan
  verbose: true

beam:
  E_e: 11.0
  E_p: 100.0

kinematics:
  x:
    min: 0.01
    max: 0.5
  Q2:
    min: 1.0
    max: 10.0

output:
  directory: ./out
  filename: generator_out.root

sullivan:
  z_default: 0.5

unweighting:
  max_weight: 0.0        # 0 => prescan will compute it
  prescan_events: 10000

Lambda:
  pT0: 0.20
  pT_max: 1.50
```

---

## Internal `GenEvent` fields (fields you can expect to inspect in ROOT writer)

Typical fields set on `GenEvent` (as seen in the code and used by writers):

* `ev.evtid` — integer event id
* `ev.x`, `ev.Q2`, `ev.xB`, `ev.y`
* `ev.theta_e`, `ev.phi_e`, `ev.k_in`, `ev.k_out`
* `ev.P_in`, `ev.Lambda`, `ev.kaon`, `ev.X_out`
* `ev.pT_L`, `ev.theta_L`, `ev.phi_L`, `ev.rapidity_L`
* `ev.t`, `ev.t_calc`, `ev.t0`
* `ev.xL`
* `ev.W2`, `ev.s`
* `ev.weight` (physics weight)
* `ev.unweighted` (0/1 flag if accepted and written)
* `ev.is_physics_based` (bool)

(Exact structure depends on your `GenEvent` struct/class and `RootWriter::fill_event_full()` implementation.)

---

## Author / contact / TODOs

* **Author/maintainer**: (insert your name/email)
* **TODOs & possible improvements**:

  * Implement proper GPD→PDF mapping or expose interface to pass a PDF already evolved to `Q²`.
  * Add Jacobian handling and allow sampling in different measures (`x`, `Q²` importance-sampling).
  * Implement HepMC/LUND/Pythia output stubs.
  * Consider sampling `t` directly (or use importance sampling for `pT` / `pz` to improve acceptance).
  * Add concurrency (threaded sampling) with careful RNG handling if needed.
  * Improve logging around weight estimation and acceptance efficiency tracking.

---
