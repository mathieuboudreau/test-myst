---
title: "Chapter 2: B1 Mapping"
---

# Chapter 2: B1 Mapping

## Introduction

In an ideal MRI experiment, the transmitted radiofrequency field (B1) would
be perfectly uniform throughout the imaging volume, producing the exact
nominal flip angle at every location. In practice, this is never achieved.
B1 inhomogeneity arises from a combination of factors that become
increasingly pronounced at higher field strengths:

- **RF coil characteristics** — the spatial sensitivity profile of any finite
  coil geometry is inherently non-uniform.
- **Electromagnetic wave effects** — at high field (≥ 3 T), the RF wavelength
  in tissue becomes comparable to the subject dimensions
  ($\lambda \approx 12$ cm at 3 T in tissue). Constructive and destructive
  interference of standing waves create large-scale B1 patterns.
- **Dielectric properties of tissue** — spatial variations in conductivity
  and permittivity locally perturb the RF field, creating hot and cold spots.

### Why B1 Inhomogeneity Matters

When the actual flip angle deviates from the nominal value, every signal
equation that assumes a specific flip angle is violated. The consequences
include:

1. **Contrast errors** in weighted images (e.g., T1-weighted GRE images
   show bright/dark bands from the B1 map imprinted on anatomy).
2. **Quantification errors** in qMRI methods that rely on the flip angle:
   - VFA T1 mapping: a 20% B1 error can produce > 50% T1 error {cite}`stollberger1996`.
   - MT saturation: the saturation pulse efficiency depends on B1².
3. **SAR miscalculation** — the deposited RF energy scales as B1², so
   accurate B1 knowledge is required for safety compliance.

### The B1+ Field

The relevant quantity is the **transmit B1 field**, denoted B1⁺ — the
circularly polarised component of the RF field that rotates in the same
sense as the nuclear spin precession. The flip angle experienced by
a spin is:

$$
\alpha_\text{actual}(\mathbf{r}) = \gamma \int_0^{\tau_p}
  |B_1^+(\mathbf{r}, t)|\,dt = F(\mathbf{r}) \cdot \alpha_\text{nominal}
$$

where $F(\mathbf{r})$ is the dimensionless **B1 scaling factor** (ideally
1 everywhere).

B1 mapping methods measure $F(\mathbf{r})$ so that downstream qMRI
analyses can correct for its spatial variation.

## Chapter Overview

Three complementary B1 mapping methods are covered in this chapter:

| Section | Method | Key feature |
|---------|--------|-------------|
| 2.1 | Double Angle Method (DAM) | Simple ratio; amplitude-based |
| 2.2 | Actual Flip-angle Imaging (AFI) | Efficient 3D; single-acquisition-pair |
| 2.3 | Bloch-Siegert Shift | Phase-based; T1-independent |

Each section describes the pulse sequence, derives the signal model,
provides an interactive simulation, and demonstrates a representative
data-fitting workflow.
