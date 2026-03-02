---
title: "Chapter 1 Conclusion"
---

# Conclusion

This chapter introduced the physical principles that underpin all MRI
experiments. Starting from the quantum mechanical concept of nuclear spin,
we derived the Larmor precession frequency and showed how a macroscopic net
magnetization arises at thermal equilibrium. The Bloch equations provided a
unified framework for describing RF excitation and the two relaxation
processes — longitudinal (T1) and transverse (T2, T2\*) — that govern
signal evolution.

We then explored three pulse sequences that exemplify different strategies
for generating image contrast:

- **Spin Echo (SE)** uses a 180° refocusing pulse to reverse the dephasing
  caused by B0 inhomogeneity. The echo amplitude decays with the true
  transverse relaxation time T2, making SE the gold standard for T2
  measurement. The CPMG variant extends this idea to a multi-echo train,
  enabling efficient T2 mapping from a single acquisition.

- **Gradient Echo (GRE)** replaces the 180° pulse with a gradient reversal.
  This dramatically shortens the minimum TR, enabling fast imaging, but
  leaves T2\* dephasing unrefocused. The Ernst equation describes the
  steady-state signal as a function of flip angle and TR/T1, and the Ernst
  angle identifies the flip angle that maximises SNR for a given TR.

- **Inversion Recovery (IR)** prepends a 180° inversion pulse to any
  readout sequence, creating powerful T1-weighted contrast. The zero-crossing
  of the longitudinal magnetisation — at $\text{TI}_\text{null} = T_1 \ln 2$
  — can be exploited to suppress specific tissues (fat in STIR, CSF in FLAIR).
  IR is also the basis for quantitative T1 mapping.

## Key Equations

| Technique | Signal equation |
|-----------|----------------|
| Spin Echo | $S = M_0 \bigl(1 - e^{-\text{TR}/T_1}\bigr)\,e^{-\text{TE}/T_2}$ |
| Gradient Echo | $S = M_0 \sin\alpha \dfrac{1-E_1}{1-E_1\cos\alpha}\,e^{-\text{TE}/T_2^*}$ |
| Inversion Recovery | $S = M_0 \bigl\lvert 1 - 2e^{-\text{TI}/T_1} + e^{-\text{TR}/T_1}\bigr\rvert$ |

## Looking Ahead

The pulse sequences and signal models introduced here are the building
blocks for the quantitative techniques covered in subsequent chapters.
In particular:

- **Chapter 2** uses GRE-based acquisitions to map the spatial variation of
  the transmitted RF field (B1).
- **Chapter 3** uses IR and VFA to produce quantitative T1 maps.
- **Chapter 4** uses multi-echo SE (CPMG) to produce quantitative T2 maps.

A thorough understanding of how TR, TE, TI, and flip angle control signal
contrast is essential for appreciating both the power and the limitations
of each quantitative mapping method.

## References

```{bibliography}
:filter: docname in docnames
```
