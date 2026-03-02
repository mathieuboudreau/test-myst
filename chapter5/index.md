---
title: "Chapter 5: Magnetization Transfer"
---

# Chapter 5: Magnetization Transfer

## Introduction

Magnetization transfer (MT) is a phenomenon arising from the interaction
between mobile water protons — the **free pool** — and protons bound within
macromolecular structures such as myelin, collagen, and proteins — the
**bound (semi-solid) pool** {cite}`wolff1989`. Although bound-pool protons
cannot be detected directly by conventional MRI (their transverse relaxation
time T2b is on the order of microseconds, producing no detectable FID), they
exert a measurable influence on the free-pool signal through magnetisation
exchange and dipolar cross-relaxation.

The bound pool has an extremely broad absorption lineshape (width ~10–100 kHz)
compared to the narrow Lorentzian of liquid water (~1 Hz in vivo). Applying an
off-resonance RF pulse — at a frequency offset $\Delta\omega$ well outside the
free-pool linewidth but within the bound-pool linewidth — saturates the bound
pool while leaving the free pool largely unaffected. This saturation
subsequently propagates to the free pool via exchange, reducing the observable
signal. The degree of signal reduction reflects the density and exchange
kinetics of the macromolecular pool.

### Physical Mechanism

The two key processes coupling the pools are:

1. **Chemical exchange** — direct exchange of individual protons between
   the two pools (e.g., hydroxyl or amide groups of macromolecules with
   water).
2. **Cross-relaxation** — transfer of longitudinal magnetisation via
   through-space dipolar interactions, without actual proton exchange.

In practice, both mechanisms are grouped into a single **exchange rate**
parameter in the two-pool model.

### Why MT Matters

MT imaging has a broad range of applications:

- **White matter integrity**: myelin sheaths are rich in bound-pool protons;
  demyelinating diseases reduce the MT effect {cite}`henkelman1993`.
- **Cartilage and musculoskeletal imaging**: dense collagen networks produce
  strong MT signals.
- **MR angiography**: MT pre-pulses suppress background tissue, enhancing
  blood vessel conspicuity.
- **Quantitative tissue characterisation**: MT parameters (bound pool
  fraction, exchange rate) are sensitive markers of microstructural
  integrity that complement T1 and T2.

### The Two-Pool Model

Quantitative MT methods describe the spin system as two coupled
compartments {cite}`henkelman1993`:

| Parameter | Free pool (a) | Bound pool (b) |
|-----------|---------------|----------------|
| T1 | T1a (~1–2 s) | T1b (~1 s) |
| T2 | T2a (~30–80 ms) | T2b (~8–15 µs) |
| Pool fraction | $1-f$ | $f$ (bound pool fraction) |
| Exchange | $k_f$ (free→bound) | $k_r = k_f/f$ (bound→free) |

The coupled **Bloch-McConnell equations** {cite}`mcconnell1958` govern
the evolution of both pools under RF irradiation and relaxation. In the
steady state under continuous-wave (CW) irradiation, the free-pool
Z-magnetisation $Z(\Delta f)$ — the MT spectrum — can be expressed in
closed form, providing the foundation for quantitative MT fitting.

### B1 Dependence

All MT methods depend on the RF transmit field B1⁺. The saturation
efficiency of any MT preparation pulse scales as $B_1^2$, meaning
that a 20 % B1 inhomogeneity (common at 3 T) introduces ~44 % error
in MT-derived parameters unless corrected. B1 maps (Chapter 2) are
therefore a prerequisite for accurate quantitative MT measurements.

## Chapter Overview

Three MT methods are presented, ordered from simplest to most
quantitative:

| Section | Method | Key feature |
|---------|--------|-------------|
| 5.1 | MT Ratio (MTR) | Semi-quantitative; single signal ratio |
| 5.2 | Quantitative MT (qMT) | Two-pool model; full Z-spectrum fitting |
| 5.3 | MT Saturation (MTsat) | Rapid approximation; three GRE contrasts |

Each section describes the pulse sequence, derives the signal model,
provides interactive simulations, and demonstrates a fitting workflow.
