---
title: "Chapter 4 Conclusion"
---

# Conclusion

This chapter presented three complementary approaches to mapping
transverse relaxation in tissue, covering both true T2 and the
effective T2*.

## Summary of Methods

**Multi-Echo Spin Echo (MESE/CPMG)** is the gold standard for T2
mapping. A 90° excitation followed by a train of 180° refocusing
pulses generates spin echoes at multiples of the inter-echo spacing TE.
The CPMG phase cycle ensures that imperfect refocusing accumulates
coherently. The monoexponential model

$$
S(n\cdot\text{TE}) = S_0\,e^{-n\cdot\text{TE}/T_2}
$$

is fitted by either log-linear regression or non-linear least squares.
Multi-component fitting of the full echo train resolves the fast myelin
water pool (T2 ≈ 10–20 ms) from intra/extra-cellular water, enabling
myelin water fraction mapping.

**Multi-Echo Gradient Echo (ME-GRE)** measures R2* = 1/T2* by
sampling the free-induction decay with multiple gradient echoes after
a single excitation. Because no refocusing pulse is applied, the signal
decays as

$$
|S(\text{TE})| = S_0\,e^{-\text{TE}/T_2^*}, \qquad
\frac{1}{T_2^*} = \frac{1}{T_2} + \gamma\,|\Delta B_0|
$$

R2* maps are sensitive to iron content, myelin, and oxygenation level
(BOLD). The complex signal phase carries $\Delta B_0$ information
exploited in quantitative susceptibility mapping (QSM).

**T2 Preparation** decouples T2 encoding from the imaging readout.
A 90°ₓ–CPMG–90°₋ₓ module converts T2 decay into initial longitudinal
magnetisation before a fast 3D readout. Adiabatic refocusing pulses make
the module robust to B1 inhomogeneity. Fitting

$$
S(\tau) = S_0\,e^{-\tau/T_2}
$$

to typically 4–6 T2prep durations yields T2 maps compatible with
breath-hold cardiac protocols.

## Comparison

| Property | MESE | ME-GRE | T2prep |
|----------|------|--------|--------|
| Measures | True T2 | T2* (and R2*) | True T2 |
| Refocuses ΔB0 | Yes | No | Yes (adiabatic) |
| B1 sensitivity | Moderate | Low | Low |
| Scan efficiency | Moderate | High | High |
| Multi-component | Yes (MWF) | No | No |
| Key application | qMRI, MWF | QSM, BOLD, R2* | Cardiac T2 |

## Connection to Other Chapters

- **B1 Mapping (Chapter 2)**: MESE T2 accuracy degrades with imperfect
  180° pulses; a B1 map and EPG fitting correct this bias.
- **T1 Mapping (Chapter 3)**: Combined T1/T2 maps support full
  relaxometry characterisation and feed into tissue segmentation.
- **Magnetization Transfer (Chapter 5)**: The two-pool qMT model
  requires T2 of the free pool (T2a) as a measured input. MESE provides
  this value directly.

## References

```{bibliography}
:filter: docname in docnames
```
