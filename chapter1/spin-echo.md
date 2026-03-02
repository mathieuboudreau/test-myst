---
title: "1.1 Spin Echo"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 1.1 Spin Echo

The spin echo (SE) is one of the oldest and most important pulse sequences
in MRI, first described by Hahn in 1950 {cite}`hahn1950`. Its defining
feature — a 180° refocusing pulse — eliminates the signal attenuation
caused by static B0 inhomogeneity, allowing the true transverse relaxation
time T2 to be measured independently of hardware imperfections.

## Pulse Sequence

A spin-echo sequence consists of two RF pulses separated by a time
$\tau = \text{TE}/2$:

1. **90° excitation pulse** — tips $M_z \to M_{xy}$, starting transverse
   dephasing.
2. **180° refocusing pulse** at time $\tau$ — inverts the phase of all
   isochromats, so that faster-precessing spins end up *behind* slower ones.
3. **Echo** at time TE = 2τ — all isochromats rephase, producing a maximum
   in the received signal.

Because static B0 offsets are refocused, the echo amplitude decays only with
the *irreversible* T2, not the faster T2\*.

```{code-cell} python
---
tags: [hide-input]
mystnb:
  figure:
    name: fig-se-sequence
    caption: |
      Spin-echo pulse sequence timing diagram. The 90° sinc excitation pulse
      (RF channel) is followed at time $\tau = \text{TE}/2$ by the 180° refocusing
      pulse. Slice-select gradients (Gz), a phase-encode blip (Gy), a pre-phaser
      and readout gradient (Gx), and the ADC acquisition window are shown.
    align: center
---

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

fig, axes = plt.subplots(4, 1, figsize=(10, 6), sharex=True)
fig.subplots_adjust(hspace=0.05)

t = np.linspace(0, 1, 1000)
TR = 1.0
TE = 0.4
tau = TE / 2

def rect(t, t0, width, height=1.0):
    return height * ((t >= t0) & (t <= t0 + width)).astype(float)

def sinc_pulse(t, t0, width, height=1.0):
    x = (t - (t0 + width / 2)) / (width / 6)
    v = height * np.sinc(x)
    v[(t < t0) | (t > t0 + width)] = 0
    return v

# RF channel
rf = sinc_pulse(t, 0.02, 0.06, 0.6) + sinc_pulse(t, tau - 0.04, 0.08, 1.0)
axes[0].fill_between(t, 0, rf, color='royalblue', alpha=0.8)
axes[0].axhline(0, color='k', linewidth=0.8)
axes[0].annotate('90°', (0.05, 0.65), ha='center', fontsize=9)
axes[0].annotate('180°', (tau, 1.05), ha='center', fontsize=9)
axes[0].set_ylabel('RF', fontsize=9)
axes[0].set_ylim(-0.2, 1.3)
axes[0].yaxis.set_tick_params(labelleft=False)

# Gz (slice select)
gz = rect(t, 0.01, 0.08, 0.7) + rect(t, 0.09, 0.04, -0.4) + \
     rect(t, tau - 0.05, 0.1, 0.7) + rect(t, tau + 0.05, 0.03, -0.4)
axes[1].fill_between(t, 0, gz, where=gz > 0, color='forestgreen', alpha=0.7)
axes[1].fill_between(t, 0, gz, where=gz < 0, color='tomato', alpha=0.7)
axes[1].axhline(0, color='k', linewidth=0.8)
axes[1].set_ylabel('Gz', fontsize=9)
axes[1].yaxis.set_tick_params(labelleft=False)

# Gy (phase encode)
gy = rect(t, 0.16, 0.06, 0.6)
axes[2].fill_between(t, 0, gy, color='darkorange', alpha=0.7)
axes[2].axhline(0, color='k', linewidth=0.8)
axes[2].set_ylabel('Gy', fontsize=9)
axes[2].yaxis.set_tick_params(labelleft=False)

# Gx (readout) + ADC
gx = rect(t, 0.16, 0.06, -0.5) + rect(t, TE - 0.08, 0.16, 0.5)
axes[3].fill_between(t, 0, gx, where=gx > 0, color='purple', alpha=0.7)
axes[3].fill_between(t, 0, gx, where=gx < 0, color='plum', alpha=0.7)
# ADC window
axes[3].fill_between(t, -0.15, -0.05,
                     where=((t >= TE - 0.08) & (t <= TE + 0.08)),
                     color='gold', alpha=0.9, label='ADC')
axes[3].axhline(0, color='k', linewidth=0.8)
axes[3].set_ylabel('Gx / ADC', fontsize=9)
axes[3].yaxis.set_tick_params(labelleft=False)

# Annotations
for ax in axes:
    ax.axvline(tau, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.axvline(TE, color='gray', linestyle=':', linewidth=0.8, alpha=0.6)

axes[3].annotate('', xy=(TE, -0.25), xytext=(0, -0.25),
                 arrowprops=dict(arrowstyle='<->', color='k'))
axes[3].text(TE / 2, -0.35, 'TE', ha='center', fontsize=9)
axes[3].annotate('', xy=(TR, -0.25), xytext=(0, -0.25),
                 arrowprops=dict(arrowstyle='<->', color='dimgray'))
axes[3].text(TR / 2, -0.45, 'TR', ha='center', fontsize=9, color='dimgray')

axes[3].set_xlim(0, TR)
axes[3].set_xlabel('Time', fontsize=9)
axes[3].set_ylim(-0.55, 0.8)

axes[0].set_title('Spin-Echo Pulse Sequence Diagram', fontsize=11, pad=8)
plt.tight_layout()
plt.show()
```

The multi-echo variant — the **CPMG** sequence {cite}`carr1954,meiboom1958` —
applies a train of 180° pulses at intervals of TE to produce multiple echoes
in a single TR, enabling efficient T2 mapping.

## Mathematical Model

### Signal Equation

At the time of the echo, the transverse magnetisation is:

$$
S(\text{TE}, \text{TR}) = M_0\,\bigl(1 - e^{-\text{TR}/T_1}\bigr)\,e^{-\text{TE}/T_2}
$$ (eq-se-signal)

The first factor describes the incomplete recovery of $M_z$ during TR; the
second describes the irreversible T2 decay between excitation and echo.

Two important limiting cases:

- **T1-weighted** ($\text{TR} \ll T_1$, short TE): contrast is dominated by
  differences in T1.
- **T2-weighted** ($\text{TR} \gg T_1$, long TE): the exponential pre-factor
  approaches 1 and contrast reflects T2 differences.
- **Proton-density weighted** ($\text{TR} \gg T_1$, short TE): both
  relaxation factors approach 1 and signal ∝ $M_0 \propto \rho_H$.

### Phase Evolution and Refocusing

In the rotating frame, a spin at position $x$ with a resonance offset
$\Delta\omega(x)$ accumulates phase $\phi = \Delta\omega(x)\,t$ after
the 90° pulse. The 180° pulse at time $\tau$ inverts all phases:
$\phi \to -\phi$. Subsequently the phase evolves as
$\phi = -\Delta\omega(x)(\tau - [t - \tau]) = \Delta\omega(x)(2\tau - t)$,
reaching zero at $t = \text{TE} = 2\tau$ — the echo. Because $\Delta\omega$
cancels exactly, the echo amplitude depends only on irreversible T2 dephasing.

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Tissue parameters (approximate values at 3 T)
tissues = {
    'White Matter': {'T1': 840,  'T2': 70,  'color': '#1f77b4'},
    'Gray Matter':  {'T1': 1300, 'T2': 90,  'color': '#ff7f0e'},
    'CSF':          {'T1': 4500, 'T2': 2000,'color': '#2ca02c'},
    'Muscle':       {'T1': 1400, 'T2': 32,  'color': '#d62728'},
}

fig = make_subplots(rows=1, cols=2,
                    subplot_titles=['Signal vs TE (TR = 3000 ms)',
                                    'Signal vs TR (TE = 15 ms)'],
                    shared_yaxes=False)

TE_arr = np.linspace(0, 400, 500)
TR_arr = np.linspace(100, 6000, 500)
TR_fixed = 3000
TE_fixed = 15

for name, p in tissues.items():
    T1, T2, col = p['T1'], p['T2'], p['color']
    sig_TE = (1 - np.exp(-TR_fixed / T1)) * np.exp(-TE_arr / T2)
    sig_TR = (1 - np.exp(-TR_arr  / T1)) * np.exp(-TE_fixed / T2)

    fig.add_trace(go.Scatter(x=TE_arr, y=sig_TE, name=name,
                             line=dict(color=col, width=2),
                             legendgroup=name, showlegend=True), row=1, col=1)
    fig.add_trace(go.Scatter(x=TR_arr, y=sig_TR, name=name,
                             line=dict(color=col, width=2),
                             legendgroup=name, showlegend=False), row=1, col=2)

fig.update_xaxes(title_text='Echo Time TE (ms)', row=1, col=1)
fig.update_xaxes(title_text='Repetition Time TR (ms)', row=1, col=2)
fig.update_yaxes(title_text='Normalised Signal (a.u.)', row=1, col=1)

fig.update_layout(
    height=420,
    title_text='Spin-Echo Signal Behaviour',
    legend=dict(orientation='h', yanchor='bottom', y=-0.25, xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows how increasing TE attenuates signal via T2 decay; tissues
with longer T2 (CSF) retain signal for longer. The right panel shows how
increasing TR allows greater T1 recovery, with slowly recovering tissues (CSF)
benefiting most from long TR.

## Data Fitting

T2 can be estimated from a multi-echo experiment by fitting
[Eq. %s](eq-se-signal) to data collected at several TE values (with TR fixed
and long enough to ensure full T1 recovery).

```{code-cell} python
:tags: [hide-input]

import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objects as go

rng = np.random.default_rng(42)

# Ground-truth parameters
T2_true = 80.0   # ms
M0_true = 1.0
TR = 5000.0       # ms >> T1, so factor ≈ 1

TE_data = np.array([12, 24, 36, 60, 90, 120, 160, 200, 260, 320])
noise_level = 0.03

# Synthetic data
signal_true = M0_true * np.exp(-TE_data / T2_true)
signal_noisy = signal_true + rng.normal(0, noise_level, size=TE_data.shape)
signal_noisy = np.clip(signal_noisy, 1e-6, None)

# Fitting
def se_model(TE, M0, T2):
    return M0 * np.exp(-TE / T2)

popt, pcov = curve_fit(se_model, TE_data, signal_noisy, p0=[1.0, 60.0],
                       bounds=([0, 1], [2, 1000]))
M0_fit, T2_fit = popt
perr = np.sqrt(np.diag(pcov))

TE_dense = np.linspace(0, 350, 500)
signal_fit = se_model(TE_dense, *popt)

# Plot
fig = go.Figure()
fig.add_trace(go.Scatter(
    x=TE_dense, y=signal_fit,
    name=f'Fit: T2 = {T2_fit:.1f} ± {perr[1]:.1f} ms',
    line=dict(color='royalblue', width=2.5)
))
fig.add_trace(go.Scatter(
    x=TE_data, y=signal_noisy,
    mode='markers',
    name=f'Synthetic data (T2_true = {T2_true} ms)',
    marker=dict(color='tomato', size=9, symbol='circle')
))
fig.add_trace(go.Scatter(
    x=TE_dense, y=np.exp(-TE_dense / T2_true),
    name='Ground truth',
    line=dict(color='gray', dash='dot', width=1.5)
))

fig.update_layout(
    title=f'Spin-Echo T2 Fitting (TR = {TR:.0f} ms)',
    xaxis_title='Echo Time TE (ms)',
    yaxis_title='Signal (a.u.)',
    template='plotly_white',
    height=380,
    legend=dict(x=0.98, y=0.98, xanchor='right', yanchor='top'),
)
fig.show()
```

The synthetic dataset (red markers) was generated at ten echo times with
Gaussian noise ($\sigma = 0.03$). A two-parameter non-linear least-squares
fit recovers $T_2$ within the statistical uncertainty determined by the
noise level and the echo time spacing. In practice, absolute signal
magnitudes are used (magnitude images), which introduces a Rician noise
floor that can bias short-T2 estimates at low SNR.

## References

```{bibliography}
:filter: docname in docnames
```
