---
title: "Chapter 4: T2 Mapping"
---

# Chapter 4: T2 Mapping

## Introduction

The transverse relaxation time T2 — also called the *spin-spin* relaxation
time — governs how quickly transverse magnetisation decays after excitation.
Two distinct mechanisms contribute to the loss of transverse coherence:

**True T2** (spin-spin relaxation) arises from random fluctuating magnetic
fields created by neighbouring spins exchanging energy. This is an irreversible
process. After an RF excitation:

$$
M_{xy}(t) = M_{xy}(0)\,e^{-t/T_2}
$$

**T2\*** (effective transverse relaxation) includes additional dephasing from
static field inhomogeneity $\Delta B_0(\mathbf{r})$:

$$
\frac{1}{T_2^*} = \frac{1}{T_2} + \frac{1}{T_2'}, \qquad
T_2' = \frac{1}{\gamma\,|\Delta B_0|}
$$

The dephasing from $\Delta B_0$ is reversible by a 180° refocusing pulse
(the basis of the spin echo). Gradient echoes are sensitive to T2*;
spin echoes refocus $\Delta B_0$ and measure true T2.

### Tissue T2 Values

At 3 T, T2 in the brain spans roughly:

| Tissue | T2 (ms) | T2* (ms) |
|--------|---------|---------|
| White matter | 70–80 | 30–40 |
| Gray matter | 90–110 | 40–50 |
| CSF | 1000–2000 | 100–200 |
| Myelin water | 10–20 | < 10 |

### Clinical and Scientific Relevance

- **Lesion detection**: multiple sclerosis plaques, oedema, and
  haemorrhage are distinguished from normal tissue by their T2/T2*.
- **Myelin water imaging**: the short-T2 water trapped between myelin
  bilayers (T2 ≈ 10–20 ms) can be resolved by multi-echo T2 fitting,
  providing a myelin water fraction (MWF) map.
- **BOLD fMRI**: the BOLD signal exploits T2* sensitivity to
  deoxyhaemoglobin concentration changes.
- **Quantitative susceptibility mapping (QSM)**: T2*-weighted multi-echo
  GRE phase data are used to compute local susceptibility maps.
- **Cartilage assessment**: T2 mapping is sensitive to cartilage
  degeneration in osteoarthritis.

## Chapter Overview

Three T2/T2* mapping methods are covered:

| Section | Method | Measures |
|---------|--------|----------|
| 4.1 | Multi-Echo Spin Echo (MESE) | True T2 |
| 4.2 | Multi-Echo Gradient Echo (ME-GRE) | T2* (and R2*) |
| 4.3 | T2 Preparation | T2-weighted contrast; T2 from LUT |
