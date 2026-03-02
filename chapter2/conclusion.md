---
title: "Chapter 2 Conclusion"
---

# Conclusion

This chapter introduced the problem of B1 field inhomogeneity in MRI and
presented three widely-used methods for mapping the transmitted RF field (B1⁺)
across the imaging volume.

## Summary of Methods

**Double Angle Method (DAM)** is conceptually simple: two spoiled GRE
acquisitions with nominal flip angles $\alpha$ and $2\alpha$ yield a signal
ratio that, in the fully-relaxed case, depends only on the actual flip angle:

$$
\alpha_\text{actual} = \arccos\!\left(\frac{S(2\alpha)}{2\,S(\alpha)}\right)
$$

DAM is robust and easy to implement, but requires long TR (typically $\geq 5T_1$)
to minimise T1 bias, making it slow for whole-brain 3D acquisitions.

**Actual Flip-angle Imaging (AFI)** achieves far better scan efficiency by
acquiring two signals within a single TR period using two interleaved time
intervals (TR1 and TR2). The signal ratio is largely T1-independent when
$n = \text{TR2}/\text{TR1} \geq 5$ and TR1 $\leq 0.5\,T_1$, enabling fast
3D B1 maps with sub-minute acquisition times.

**Bloch-Siegert Shift** encodes B1 into the image phase rather than magnitude.
An off-resonance Fermi pulse induces a phase shift:

$$
\Delta\phi = 2\,K_\text{BS}\,(B_1^+)^2
$$

This approach is inherently T1-independent and insensitive to magnitude
variations, but requires accurate knowledge of the pulse shape constant
$K_\text{BS}$ and careful phase unwrapping.

## Comparison

| Property | DAM | AFI | Bloch-Siegert |
|----------|-----|-----|---------------|
| Signal domain | Magnitude | Magnitude | Phase |
| T1 dependence | High | Low | None |
| Acquisition time | Slow | Fast | Moderate |
| SAR | Low | Low | Higher (BS pulse) |
| Implementation complexity | Low | Moderate | Moderate |

## Why B1 Mapping Matters for qMRI

Every quantitative MRI method that relies on a specific flip angle is
vulnerable to B1 inhomogeneity. In particular:

- **VFA T1 mapping** (Chapter 3): the Ernst equation assumes a known flip
  angle; B1 errors propagate to T1 errors super-linearly. B1 correction
  is essential for accurate VFA T1 maps, especially at 3 T and above.
- **Magnetization Transfer** (Chapter 5): the MT saturation efficiency
  scales as $B_1^2$, so an uncharacterised 20% B1 deviation causes ~44%
  error in the MT saturation parameter.

Acquiring a B1 map as part of every quantitative protocol is increasingly
recognised as a best practice in the qMRI community.

## References

```{bibliography}
:filter: docname in docnames
```
