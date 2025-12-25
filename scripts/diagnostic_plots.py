# import uproot
# import awkward as ak
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages

# # -------------------------
# # Configuration
# # -------------------------
# ROOT_FILE = "output/events_no_models.root"
# TREE_NAME = "events"   # change if your tree name differs
# OUTPUT_PDF = "output/generator_diagnostics.pdf"

# # Physics constants
# MP = 0.938272
# ML = 1.115683

# # -------------------------
# # Helper plotting functions
# # -------------------------
# def plot_1d(ax, data, xlabel, title, bins=100, range=None, logy=False):
#     ax.hist(data, bins=bins, range=range, histtype="step", linewidth=1.5)
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel("Counts")
#     ax.set_title(title)
#     if logy:
#         ax.set_yscale("log")
#     ax.grid(alpha=0.3)

# def plot_2d(ax, x, y, xlabel, ylabel, title, bins=100, range=None):
#     h = ax.hist2d(x, y, bins=bins, range=range, cmap="viridis")
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#     ax.set_title(title)
#     plt.colorbar(h[3], ax=ax, label="Counts")

# # -------------------------
# # Load ROOT file
# # -------------------------
# with uproot.open(ROOT_FILE) as f:
#     tree = f[TREE_NAME]
#     data = tree.arrays(library="ak")

# # -------------------------
# # Extract variables
# # -------------------------
# x     = ak.to_numpy(data["x"])
# Q2    = ak.to_numpy(data["Q2"])
# t     = ak.to_numpy(data["t"])
# xB    = ak.to_numpy(data["xB"])
# y     = ak.to_numpy(data["y"])
# W2    = ak.to_numpy(data["W2"])
# s     = ak.to_numpy(data["s"])
# xL    = ak.to_numpy(data["xL"])
# pT_L  = ak.to_numpy(data["pT_L"])

# theta_e = ak.to_numpy(data["theta_e"])
# phi_e   = ak.to_numpy(data["phi_e"])
# theta_L = ak.to_numpy(data["theta_L"])
# phi_L   = ak.to_numpy(data["phi_L"])
# rap_L   = ak.to_numpy(data["rapidity_L"])

# # Lambda 4-vector diagnostics
# EL = ak.to_numpy(data["Lambda_E"])
# pL = np.sqrt(
#     ak.to_numpy(data["Lambda_px"])**2 +
#     ak.to_numpy(data["Lambda_py"])**2 +
#     ak.to_numpy(data["Lambda_pz"])**2
# )
# ML2 = EL**2 - pL**2

# # -------------------------
# # Derived diagnostics
# # -------------------------
# MX2 = (
#     ak.to_numpy(data["X_out_E"])**2
#     - ak.to_numpy(data["X_out_px"])**2
#     - ak.to_numpy(data["X_out_py"])**2
#     - ak.to_numpy(data["X_out_pz"])**2
# )


# # -------------------------
# # Create PDF
# # -------------------------
# with PdfPages(OUTPUT_PDF) as pdf:

#     # -------- 1D core DIS --------
#     fig, axs = plt.subplots(2, 2, figsize=(10, 8))
#     plot_1d(axs[0,0], x,  "x",  "Bjorken x", range=(0,1))
#     plot_1d(axs[0,1], Q2, "Q² [GeV²]", "Q²", logy=True)
#     plot_1d(axs[1,0], y,  "y",  "Inelasticity y", range=(0,1))
#     plot_1d(axs[1,1], W2, "W² [GeV²]", "Invariant mass W²")
#     plt.tight_layout()
#     pdf.savefig(fig)
#     plt.close()

#     # -------- t, xL, pT --------
#     fig, axs = plt.subplots(2, 2, figsize=(10, 8))
#     plot_1d(axs[0,0], t,    "t [GeV²]", "|t| distribution")
#     plot_1d(axs[0,1], xL,   "x_L", "Lambda momentum fraction", range=(0,1))
#     plot_1d(axs[1,0], pT_L, "pT_L [GeV]", "Lambda transverse momentum")
#     plot_1d(axs[1,1], MX2,  "M_X² [GeV²]", "Hadronic remnant mass²")
#     plt.tight_layout()
#     pdf.savefig(fig)
#     plt.close()

#     # -------- Angular distributions --------
#     fig, axs = plt.subplots(2, 2, figsize=(10, 8))
#     plot_1d(axs[0,0], theta_e, "θₑ [rad]", "Electron polar angle")
#     plot_1d(axs[0,1], phi_e,   "φₑ [rad]", "Electron azimuth")
#     plot_1d(axs[1,0], theta_L, "θ_Λ [rad]", "Lambda polar angle")
#     plot_1d(axs[1,1], phi_L,   "φ_Λ [rad]", "Lambda azimuth")
#     plt.tight_layout()
#     pdf.savefig(fig)
#     plt.close()

#     # -------- 2D correlations --------
#     fig, axs = plt.subplots(2, 2, figsize=(10, 8))
#     plot_2d(axs[0,0], x, Q2, "x", "Q²", "x vs Q²")
#     plot_2d(axs[0,1], xL, t, "x_L", "t", "t vs x_L")
#     plot_2d(axs[1,0], pT_L, xL, "pT_L", "x_L", "x_L vs pT_L")
#     plot_2d(axs[1,1], MX2, Q2, "M_X²", "Q²", "Q² vs M_X²")
#     plt.tight_layout()
#     pdf.savefig(fig)
#     plt.close()

#     # -------- Energy–momentum sanity --------
#     fig, axs = plt.subplots(1, 2, figsize=(10, 4))
#     plot_1d(axs[0], EL, "E_Λ [GeV]", "Lambda energy")
#     plot_1d(axs[1], pL, "|p_Λ| [GeV]", "Lambda momentum magnitude")
#     plt.tight_layout()
#     pdf.savefig(fig)
#     plt.close()

# print(f"[Diagnostics] Saved plots to {OUTPUT_PDF}")

#!/usr/bin/env python3
"""
diagnostic_plots.py

Reads ROOT events (uproot + awkward), computes diagnostics and
writes a multipage PDF with 1D and 2D histograms.

Auto-selects the highest-cycle tree if ROOT stored cycles like 'events;8'.
"""

import os
import sys
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from math import sqrt

# -------------------------
# Configuration
# -------------------------
ROOT_FILE = "output/events_no_models.root"
TREE_NAME = "events"   # base name (script will auto-select events;N if present)
OUTPUT_PDF = "output/generator_diagnostics.pdf"

# Physics constants (GeV)
MP = 0.938272
ML = 1.115683

# Plotting defaults
N_BINS = 100

# -------------------------
# Helper utilities
# -------------------------
def get_branch(tree, name):
    """Return a flattened numpy array for branch 'name' if present else None."""
    if name in tree.keys():
        arr = tree[name].array(library="ak")
        try:
            flat = ak.to_numpy(ak.ravel(arr))
        except Exception:
            flat = np.asarray(arr.tolist(), dtype=float)
        return flat
    else:
        print(f"[warning] branch '{name}' not found in tree; skipping related diagnostics.")
        return None

def safe_hist(ax, data, bins=N_BINS, rng=None, label=None, logy=False):
    if data is None or data.size == 0:
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return
    ax.hist(data, bins=bins, range=rng, histtype="step", linewidth=1.5, label=label)
    if label:
        ax.legend()
    if logy:
        ax.set_yscale("log")
    ax.grid(alpha=0.3)

def safe_hist2d(ax, x, y, bins=100, rng=None):
    if x is None or y is None or x.size == 0 or y.size == 0:
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return
    h = ax.hist2d(x, y, bins=bins, range=rng, cmap="viridis")
    plt.colorbar(h[3], ax=ax, label="Counts")
    ax.grid(alpha=0.3)

def percentile_range(arr, low_pct=0.5, high_pct=99.5, minwidth=1e-6):
    """Return (low, high) based on percentiles, with a minimum width."""
    if arr is None or arr.size == 0:
        return None
    lo = np.nanpercentile(arr, low_pct)
    hi = np.nanpercentile(arr, high_pct)
    if hi - lo < minwidth:
        hi = lo + minwidth
    return (lo, hi)

# -------------------------
# Open file and auto-select tree (handle cycles like 'events;8')
# -------------------------
if not os.path.exists(ROOT_FILE):
    print(f"[error] ROOT file not found: {ROOT_FILE}")
    sys.exit(1)

with uproot.open(ROOT_FILE) as f:
    available = list(f.keys())
    # find candidates whose base name equals TREE_NAME
    candidates = []
    for k in available:
        base = k.split(';')[0]
        if base == TREE_NAME:
            cycle = 0
            if ';' in k:
                try:
                    cycle = int(k.split(';')[1])
                except Exception:
                    cycle = 0
            candidates.append((cycle, k))
    if candidates:
        # pick the highest cycle (largest integer)
        candidates.sort()
        chosen_key = candidates[-1][1]
        tree = f[chosen_key]
        print(f"[info] Using tree key '{chosen_key}' (auto-selected from file).")
    else:
        # fallback: if exact match present use it, else show available keys and exit
        if TREE_NAME in available:
            tree = f[TREE_NAME]
            print(f"[info] Using tree key '{TREE_NAME}' (exact match).")
        else:
            print(f"[error] Tree '{TREE_NAME}' not found in file. Available keys:\n{available}")
            sys.exit(1)

    # show sample branches
    branches_available = list(tree.keys())
    print(f"[info] Found branches (sample): {branches_available[:40]} ...")

    # -------------------------
    # Extract branches (safe)
    # -------------------------
    x       = get_branch(tree, "x")
    Q2      = get_branch(tree, "Q2")
    t       = get_branch(tree, "t")
    xB      = get_branch(tree, "xB")
    y       = get_branch(tree, "y")
    W2      = get_branch(tree, "W2")
    s_arr   = get_branch(tree, "s")
    xL      = get_branch(tree, "xL")
    pT_L    = get_branch(tree, "pT_L")

    theta_e = get_branch(tree, "theta_e")
    phi_e   = get_branch(tree, "phi_e")
    theta_L = get_branch(tree, "theta_L")
    phi_L   = get_branch(tree, "phi_L")
    rap_L   = get_branch(tree, "rapidity_L")

    # Lambda four-vector components (your names)
    EL      = get_branch(tree, "Lambda_E")
    pxL     = get_branch(tree, "Lambda_px")
    pyL     = get_branch(tree, "Lambda_py")
    pzL     = get_branch(tree, "Lambda_pz")

    # hadronic remnant X (user-provided)
    X_E     = get_branch(tree, "X_out_E")
    X_px    = get_branch(tree, "X_out_px")
    X_py    = get_branch(tree, "X_out_py")
    X_pz    = get_branch(tree, "X_out_pz")

    # incoming/outgoing electron and proton components (if available)
    k_in_E  = get_branch(tree, "k_in_E")
    k_in_px = get_branch(tree, "k_in_px")
    k_in_py = get_branch(tree, "k_in_py")
    k_in_pz = get_branch(tree, "k_in_pz")

    k_out_E  = get_branch(tree, "k_out_E")
    k_out_px = get_branch(tree, "k_out_px")
    k_out_py = get_branch(tree, "k_out_py")
    k_out_pz = get_branch(tree, "k_out_pz")

    P_in_E  = get_branch(tree, "P_in_E")
    P_in_px = get_branch(tree, "P_in_px")
    P_in_py = get_branch(tree, "P_in_py")
    P_in_pz = get_branch(tree, "P_in_pz")

# -------------------------
# Derived arrays & safety masks
# -------------------------
def safe_vec_mag(px, py, pz):
    if px is None or py is None or pz is None:
        return None
    return np.sqrt(np.maximum(0.0, px*px + py*py + pz*pz))

# Lambda momentum magnitude & invariant (if components available)
pL = safe_vec_mag(pxL, pyL, pzL)
ML2 = None
if EL is not None and pL is not None:
    ML2 = EL**2 - pL**2

# MX2 from X_out branches (if present)
MX2 = None
MX_pos = None
if X_E is not None and X_px is not None and X_py is not None and X_pz is not None:
    MX2 = X_E**2 - (X_px**2 + X_py**2 + X_pz**2)
    MX_pos = np.where(MX2 > 0, np.sqrt(MX2), np.nan)
else:
    print("[info] Missing X_out_* branches; MX2 diagnostics will be skipped.")

# four-momentum residual: total_in - total_out - X_out  (should be ~0 if X_out is consistent)
residual4 = None
if (k_in_E is not None and P_in_E is not None and k_out_E is not None and EL is not None
    and k_in_px is not None and P_in_px is not None and k_out_px is not None and pxL is not None
    and X_E is not None and X_px is not None):
    diff_E = (k_in_E + P_in_E) - (k_out_E + EL) - X_E
    diff_px = (k_in_px + P_in_px) - (k_out_px + pxL) - X_px
    diff_py = (k_in_py + P_in_py) - (k_out_py + pyL) - X_py
    diff_pz = (k_in_pz + P_in_pz) - (k_out_pz + pzL) - X_pz
    residual4 = np.sqrt(diff_E**2 + diff_px**2 + diff_py**2 + diff_pz**2)
    print(f"[info] residual4: min={np.nanmin(residual4):.3e}, median={np.nanmedian(residual4):.3e}, max={np.nanmax(residual4):.3e}")
else:
    print("[info] Missing components for 4-vector residual check; skipping residual plot.")

# Build a global physical mask that we can apply to several plots:
n_events = None
for arr in (x, Q2, t, xB, y, W2, s_arr, EL):
    if arr is not None:
        n_events = arr.size
        break
if n_events is None:
    print("[error] No usable branches found; aborting.")
    sys.exit(1)
mask = np.ones(n_events, dtype=bool)

# apply basic physical constraints where available
if x is not None:
    mask &= np.isfinite(x) & (x > 0.0) & (x < 1.0)
if Q2 is not None:
    mask &= np.isfinite(Q2) & (Q2 > 0.0)
if y is not None:
    mask &= np.isfinite(y) & (y >= 0.0) & (y <= 1.0)
if EL is not None:
    mask &= np.isfinite(EL) & (EL >= ML)

# utility to apply mask safely
def apply_mask(arr):
    if arr is None:
        return None
    if arr.size != n_events:
        try:
            arr2 = np.asarray(arr)
            if arr2.size == n_events:
                arr = arr2
            else:
                return None
        except Exception:
            return None
    return arr[mask]

# masked versions for plotting
x_m   = apply_mask(x)
Q2_m  = apply_mask(Q2)
t_m   = apply_mask(t)
xB_m  = apply_mask(xB)
y_m   = apply_mask(y)
W2_m  = apply_mask(W2)
s_m   = apply_mask(s_arr)
xL_m  = apply_mask(xL)
pT_m  = apply_mask(pT_L)

theta_e_m = apply_mask(theta_e)
phi_e_m   = apply_mask(phi_e)
theta_L_m = apply_mask(theta_L)
phi_L_m   = apply_mask(phi_L)
rap_L_m   = apply_mask(rap_L)

EL_m  = apply_mask(EL)
pL_m  = apply_mask(pL)
ML2_m = apply_mask(ML2) if ML2 is not None else None

MX2_m = apply_mask(MX2) if MX2 is not None else None
MX_pos_m = apply_mask(MX_pos) if MX2 is not None else None

residual4_m = apply_mask(residual4) if residual4 is not None else None

# -------------------------
# Create PDF with diagnostics
# -------------------------
os.makedirs(os.path.dirname(OUTPUT_PDF), exist_ok=True)
with PdfPages(OUTPUT_PDF) as pdf:

    # 1) Core DIS 1D
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    safe_hist(axs[0,0], x_m if x_m is not None else np.array([]), bins=80, rng=(0,1))
    axs[0,0].set_title("Bjorken x")
    axs[0,0].set_xlabel("x")
    safe_hist(axs[0,1], Q2_m if Q2_m is not None else np.array([]), bins=80, rng=percentile_range(Q2_m) if Q2_m is not None else None, logy=True)
    axs[0,1].set_title("Q² [GeV²]")
    safe_hist(axs[1,0], y_m if y_m is not None else np.array([]), bins=80, rng=(0,1))
    axs[1,0].set_title("inelasticity y")
    safe_hist(axs[1,1], W2_m if W2_m is not None else np.array([]), bins=80, rng=percentile_range(W2_m))
    axs[1,1].set_title("W² [GeV²]")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()

    # 2) t, xL, pT, MX2
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    safe_hist(axs[0,0], t_m if t_m is not None else np.array([]), bins=80, rng=percentile_range(t_m))
    axs[0,0].set_title("t [GeV²]")
    safe_hist(axs[0,1], xL_m if xL_m is not None else np.array([]), bins=80, rng=(0,1))
    axs[0,1].set_title("x_L")
    safe_hist(axs[1,0], pT_m if pT_m is not None else np.array([]), bins=80, rng=percentile_range(pT_m))
    axs[1,0].set_title("pT_L [GeV]")
    if MX2_m is not None:
        safe_hist(axs[1,1], MX2_m, bins=100, rng=percentile_range(MX2_m, 0.5, 99.5))
        axs[1,1].set_title("M_X² [GeV²] (negative allowed)")
    else:
        axs[1,1].text(0.5, 0.5, "MX2 not available", ha="center", va="center")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()

    # 3) Angular distributions
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    safe_hist(axs[0,0], theta_e_m if theta_e_m is not None else np.array([]), bins=80, rng=percentile_range(theta_e_m))
    axs[0,0].set_title("theta_e [rad]")
    safe_hist(axs[0,1], phi_e_m if phi_e_m is not None else np.array([]), bins=80, rng=(0, 2*np.pi))
    axs[0,1].set_title("phi_e [rad]")
    safe_hist(axs[1,0], theta_L_m if theta_L_m is not None else np.array([]), bins=80, rng=percentile_range(theta_L_m))
    axs[1,0].set_title("theta_L [rad]")
    safe_hist(axs[1,1], phi_L_m if phi_L_m is not None else np.array([]), bins=80, rng=(0, 2*np.pi))
    axs[1,1].set_title("phi_L [rad]")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()

    # 4) 2D correlations (x vs Q2, xL vs t, pT vs xL, MX2 vs Q2)
    fig, axs = plt.subplots(2, 2, figsize=(11, 9))
    safe_hist2d(axs[0,0], x_m, Q2_m, bins=(80,80), rng=((0,1), percentile_range(Q2_m)))
    axs[0,0].set_title("x vs Q²")
    safe_hist2d(axs[0,1], xL_m, t_m, bins=(80,80), rng=((0,1), percentile_range(t_m)))
    axs[0,1].set_title("t vs x_L")
    safe_hist2d(axs[1,0], pT_m, xL_m, bins=(80,80), rng=(percentile_range(pT_m), (0,1)))
    axs[1,0].set_title("pT_L vs x_L")
    safe_hist2d(axs[1,1], MX2_m, Q2_m, bins=(80,80), rng=(percentile_range(MX2_m), percentile_range(Q2_m)))
    axs[1,1].set_title("M_X² vs Q²")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()

    # 5) Lambda energy & momentum
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    safe_hist(axs[0], EL_m if EL_m is not None else np.array([]), bins=80, rng=percentile_range(EL_m))
    axs[0].set_title("E_L [GeV]")
    safe_hist(axs[1], pL_m if pL_m is not None else np.array([]), bins=80, rng=percentile_range(pL_m))
    axs[1].set_title("|p_L| [GeV]")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()

    # 6) MX positive-only and MX negative diagnostics (if available)
    if MX2_m is not None:
        fig, axs = plt.subplots(1, 2, figsize=(11, 4))
        safe_hist(axs[0], MX_pos_m[~np.isnan(MX_pos_m)] if MX_pos_m is not None else np.array([]), bins=80, rng=percentile_range(MX_pos_m))
        axs[0].set_title("M_X (sqrt of positive M_X²)")
        neg_mask = (MX2_m < 0)
        axs[1].hist(MX2_m[neg_mask], bins=80)
        axs[1].set_title(f"Negative M_X² (count={np.sum(neg_mask)})")
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

    # 7) Four-momentum residual diagnostics (if available)
    if residual4_m is not None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        safe_hist(ax, residual4_m, bins=100, rng=percentile_range(residual4_m, 0.5, 99.5), logy=True)
        ax.set_title("4-vector residual magnitude (Euclidean) [GeV]")
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

    # 8) Quick pairwise checks (a few important 2D)
    fig, axs = plt.subplots(2, 2, figsize=(11, 9))
    safe_hist2d(axs[0,0], x_m, xL_m, bins=(80,80), rng=((0,1),(0,1)))
    axs[0,0].set_title("x vs x_L")
    safe_hist2d(axs[0,1], pT_m, t_m, bins=(80,80), rng=(percentile_range(pT_m), percentile_range(t_m)))
    axs[0,1].set_title("t vs pT_L")
    safe_hist2d(axs[1,0], x_m, MX2_m, bins=(80,80), rng=((0,1), percentile_range(MX2_m)))
    axs[1,0].set_title("x vs M_X²")
    safe_hist2d(axs[1,1], Q2_m, MX2_m, bins=(80,80), rng=(percentile_range(Q2_m), percentile_range(MX2_m)))
    axs[1,1].set_title("Q² vs M_X²")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()

print(f"[Diagnostics] Saved plots to {OUTPUT_PDF}")
