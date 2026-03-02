---
title: "Chapter 1: A Quick Primer of MRI Physics"
---

# Chapter 1: A Quick Primer of MRI Physics

## Introduction

Magnetic Resonance Imaging (MRI) is a non-ionizing medical imaging technique
that exploits the quantum mechanical property of nuclear spin to generate images
with exquisite soft-tissue contrast. Unlike X-ray computed tomography, MRI does
not use ionizing radiation and can produce contrast sensitive to a wide variety
of tissue properties simply by altering the pulse sequence timing parameters.

This chapter provides the physical foundation needed to understand the
quantitative MRI methods covered in the rest of the course. We begin with
nuclear spin and magnetization, then describe how RF pulses tip the
magnetization into the transverse plane to produce a detectable signal, and
finally discuss how two complementary relaxation processes — longitudinal (T1)
and transverse (T2, T2\*) — determine the time course of that signal.

### Nuclear Spin and the Larmor Frequency

The hydrogen nucleus (a single proton) possesses an intrinsic quantum
mechanical angular momentum called **spin**, with a spin quantum number
$I = \tfrac{1}{2}$. Associated with this spin is a magnetic dipole moment
$\boldsymbol{\mu} = \gamma \mathbf{I}$, where $\gamma$ is the
**gyromagnetic ratio** ($\gamma / 2\pi = 42.577\ \text{MHz T}^{-1}$ for
$^1\text{H}$).

When placed in an external static magnetic field $B_0$ (oriented along
$\hat{z}$), individual proton spins precess about $\hat{z}$ at the **Larmor
frequency**:

$$
\omega_0 = \gamma B_0
$$

For a clinical 3 T scanner, $f_0 = \omega_0 / 2\pi \approx 128\ \text{MHz}$,
which lies in the RF band.

### Net Magnetization

A macroscopic tissue voxel contains on the order of $10^{22}$ protons. At
thermal equilibrium, the slight Boltzmann preference for the low-energy
spin-up state produces a small but measurable net magnetization along $B_0$:

$$
M_0 = \frac{\rho \hbar^2 \gamma^2 B_0}{4 k_\text{B} T}
$$

where $\rho$ is the proton density, $k_\text{B}$ is Boltzmann's constant, and
$T$ is absolute temperature. This net magnetization $\mathbf{M}_0 = M_0
\hat{z}$ is the source of all MRI signal.

### The Bloch Equations

The time evolution of $\mathbf{M}$ in the presence of $B_0$, an RF field
$\mathbf{B}_1$, and relaxation is described by the **Bloch equations**
{cite}`bloch1946`:

$$
\frac{d\mathbf{M}}{dt} = \gamma\,\mathbf{M} \times \mathbf{B}
  - \frac{M_x \hat{x} + M_y \hat{y}}{T_2}
  - \frac{(M_z - M_0)\hat{z}}{T_1}
$$

In the rotating frame at $\omega_0$, the RF field appears stationary, and the
Bloch equations simplify to a tractable set of linear ODEs.

### Relaxation

Following RF excitation two independent relaxation mechanisms return
$\mathbf{M}$ to equilibrium:

**Longitudinal (T1) relaxation** — also called *spin-lattice* relaxation —
describes the recovery of $M_z$ toward $M_0$:

$$
M_z(t) = M_0 \left(1 - e^{-t/T_1}\right)
$$

T1 reflects energy exchange between spins and their molecular environment
(the lattice) and depends strongly on field strength and tissue type.
Typical brain values at 3 T: WM $\approx 840\ \text{ms}$,
GM $\approx 1300\ \text{ms}$, CSF $\approx 4500\ \text{ms}$.

**Transverse (T2) relaxation** — also called *spin-spin* relaxation —
describes the dephasing of the transverse magnetization $M_{xy}$:

$$
M_{xy}(t) = M_{xy}(0)\,e^{-t/T_2}
$$

T2 arises from irreversible, random spin-spin interactions. An additional
contribution from B0 inhomogeneity yields the faster **effective transverse
relaxation time** T2\*:

$$
\frac{1}{T_2^*} = \frac{1}{T_2} + \frac{1}{T_2'}
$$

where $T_2'$ captures the (reversible) dephasing due to local field
inhomogeneity. The spin-echo technique can refocus $T_2'$ effects and
measure true T2.

### From Signal to Image

After excitation, the precessing transverse magnetization induces a voltage
in a nearby receiver coil — the free induction decay (FID). Spatial encoding
is achieved by superimposing time-varying linear gradient fields
$\mathbf{G} = (G_x, G_y, G_z)$ on $B_0$, which cause spins at different
locations to precess at different frequencies. The received signal
$s(t)$ maps out **k-space**, and a Fourier transform yields the image.

## Chapter Overview

This chapter covers three fundamental pulse sequences that form the basis
of most MRI contrast mechanisms:

| Section | Technique | Key feature |
|---------|-----------|-------------|
| 1.1 | Spin Echo (SE) | Refocuses B0 inhomogeneity; measures T2 |
| 1.2 | Gradient Echo (GRE) | Fast acquisition; sensitive to T2\* |
| 1.3 | Inversion Recovery (IR) | T1 contrast; tissue nulling |

Each section walks through the pulse sequence design, signal equation,
interactive simulation, and a worked data-fitting example.
