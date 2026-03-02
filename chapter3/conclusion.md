---
title: "Chapter 3 Conclusion"
---

# Conclusion

This chapter introduced three methods for mapping the longitudinal
relaxation time T1, spanning from the classical gold standard to modern
rapid whole-brain techniques.

## Summary of Methods

**Inversion Recovery (IR)** inverts the magnetisation and samples the
subsequent exponential recovery at multiple inversion times TI. Fitting
the three-parameter model

$$
S(\text{TI}) = A\,\bigl|1 - B\,e^{-\text{TI}/T_1}\bigr|
$$

yields a robust, direct T1 estimate. IR is the reference standard for
T1 accuracy but is slow (typically 20–40 min for whole-brain 3D) and
must be corrected for imperfect inversion efficiency when B1 is
inhomogeneous.

**Variable Flip Angle (VFA)** exploits the flip-angle dependence of
the spoiled-GRE steady-state signal. The Deoni linearisation

$$
\frac{S}{\sin\alpha} = E_1\,\frac{S}{\tan\alpha} + M_0(1 - E_1)
$$

extracts T1 from the slope of a straight line, enabling fast 3D
acquisitions (< 5 min). The critical weakness is extreme sensitivity to
B1 inhomogeneity: VFA without a co-registered B1 map is unreliable at
3 T and above.

**MP2RAGE** acquires two GRE images within a single inversion recovery
and forms the ratio

$$
\text{MP2RAGE} = \frac{\text{GRE}_1 \cdot \text{GRE}_2}
                      {|\text{GRE}_1|^2 + |\text{GRE}_2|^2}
$$

This ratio cancels M0 and the receive-field profile and is substantially
less sensitive to B1 transmit inhomogeneity than VFA, making it the
preferred 3D T1 mapping method at 7 T. T1 is extracted by inverting a
pre-computed look-up table.

## Comparison

| Property | IR | VFA | MP2RAGE |
|----------|-----|-----|---------|
| Scan time (3D whole-brain) | Very long | Fast | Moderate |
| B1 transmit sensitivity | Moderate | Very high | Low |
| Receive-field bias | Present | Present | Eliminated |
| Requires B1 map | No (but helpful) | Yes (mandatory) | No |
| T1 accuracy | Highest | High (with B1 corr.) | High |
| Implementation | Simple | Simple | Moderate |

## Connection to Other Chapters

T1 mapping connects tightly to the rest of the course:

- **B1 Mapping (Chapter 2)**: VFA requires a B1 map for correction;
  B1 maps in turn use T1 estimates to correct finite-TR biases in DAM.
- **Magnetization Transfer (Chapter 5)**: T1a is a key free parameter
  in two-pool qMT fitting. Providing a T1 map as a prior significantly
  stabilises the Z-spectrum fit.
- **T2 Mapping (Chapter 4)**: Multi-echo T2 sequences are often combined
  with a T1 map to provide a full relaxometry profile.

## References

```{bibliography}
:filter: docname in docnames
```
