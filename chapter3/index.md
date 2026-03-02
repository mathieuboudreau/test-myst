---
title: "Chapter 3: T1 Mapping"
---

# Chapter 3: T1 Mapping

## Introduction

The longitudinal relaxation time T1 — also called the *spin-lattice* relaxation
time — describes how rapidly the net magnetisation recovers toward its thermal
equilibrium value $M_0$ after being displaced by an RF pulse:

$$
M_z(t) = M_0\left(1 - e^{-t/T_1}\right)
$$

T1 depends on the efficiency of energy exchange between excited spins and
their molecular environment (the lattice). It is sensitive to molecular
mobility and, in tissue, to field strength, macromolecular content, and the
local chemical environment. At 3 T, T1 spans roughly 400 ms in lipid-rich
white matter to 4500 ms in cerebrospinal fluid.

### Why Quantify T1?

Accurate T1 maps underpin a wide range of quantitative MRI applications:

- **Contrast-agent pharmacokinetics**: gadolinium shortens T1 locally;
  a pre-contrast T1 map is required for dynamic contrast-enhanced (DCE)
  MRI to quantify tracer concentration {cite}`haacke1999`.
- **B1 correction**: variable flip-angle T1 mapping requires a B1 map
  (Chapter 2) to correct flip-angle errors; T1 maps likewise feed back
  into B1-map quality checks.
- **Myelin water imaging and MT**: T1a is a free parameter in two-pool
  MT models (Chapter 5).
- **MR fingerprinting**: simultaneous T1/T2 maps from a single
  pseudo-random sequence.
- **Brain development and ageing**: T1 decreases as myelination
  progresses and increases in demyelinating lesions.

### T1 in the Steady State

Most clinical MRI uses sequences that do not allow full T1 recovery
between excitations. The steady-state longitudinal magnetisation for a
spoiled GRE with repetition time TR and flip angle $\alpha$ is given by
the **Ernst equation** {cite}`ernst1966`:

$$
M_z^{\text{SS}} = M_0\,\frac{1 - E_1}{1 - E_1\cos\alpha}, \qquad
E_1 = e^{-\text{TR}/T_1}
$$

This expression is the starting point for both the inversion recovery
and variable flip angle T1 mapping methods.

## Chapter Overview

Three T1 mapping methods are covered in this chapter, ordered from the
classical gold standard to modern rapid approaches:

| Section | Method | Key feature |
|---------|--------|-------------|
| 3.1 | Inversion Recovery (IR) | Gold standard; direct T1 measurement |
| 3.2 | Variable Flip Angle (VFA) | Fast 3D; requires B1 map |
| 3.3 | MP2RAGE | B1-robust; simultaneous T1 + uniform-contrast image |

Each section describes the pulse sequence, derives the signal model,
provides interactive simulations, and demonstrates a data-fitting workflow.
