---
title: "1.2 Gradient Echo"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 1.2 Gradient Echo

The gradient echo (GRE) sequence replaces the 180° refocusing pulse of the
spin echo with a simple gradient reversal. This allows much shorter repetition
times (TR) and is the workhorse of fast 3D MRI. The absence of the 180° pulse
means the echo is sensitive to T2\* rather than T2, making GRE useful for
applications such as susceptibility-weighted imaging and BOLD fMRI.

## Pulse Sequence

A basic GRE sequence applies a flip angle $\alpha$ (typically $< 90°$), waits
for the echo time TE, then reads out the echo and repeats after TR:

1. **$\alpha$° excitation pulse** — partially tips $M_z$ into the transverse
   plane.
2. **Readout pre-phaser** (negative $G_x$ lobe) — dephases spins to set the
   starting point of k-space.
3. **Frequency-encode gradient** (positive $G_x$ lobe) — rephases spins,
   forming the **gradient echo** at the centre of the readout window (time TE).
4. **Spoiling** (gradient or RF) — eliminates residual transverse magnetisation
   before the next TR.

```{code-cell} python
---
tags: [hide-input]
mystnb:
  figure:
    name: fig-gre-sequence
    caption: |
      Spoiled gradient-echo pulse sequence. The small flip angle α excitation is
      followed by a negative Gx pre-phaser; the positive readout lobe creates a
      gradient echo at time TE. A spoiler gradient at the end of TR destroys
      residual transverse magnetisation.
    align: center
---

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(4, 1, figsize=(10, 5.5), sharex=True)
fig.subplots_adjust(hspace=0.05)

t = np.linspace(0, 1, 1000)
TR = 1.0
TE = 0.25
alpha_label = 'α (<90°)'

def rect(t, t0, width, height=1.0):
    r = np.zeros_like(t)
    r[(t >= t0) & (t <= t0 + width)] = height
    return r

def sinc_pulse(t, t0, width, height=1.0):
    x = (t - (t0 + width / 2)) / (width / 6)
    v = height * np.sinc(x)
    v[(t < t0) | (t > t0 + width)] = 0
    return v

# RF
rf = sinc_pulse(t, 0.02, 0.06, 0.55)
axes[0].fill_between(t, 0, rf, color='royalblue', alpha=0.8)
axes[0].axhline(0, color='k', linewidth=0.8)
axes[0].annotate(alpha_label, (0.05, 0.58), ha='center', fontsize=9)
axes[0].set_ylabel('RF', fontsize=9)
axes[0].set_ylim(-0.15, 0.85)
axes[0].yaxis.set_tick_params(labelleft=False)

# Gz
gz = rect(t, 0.01, 0.08, 0.7) + rect(t, 0.09, 0.04, -0.4)
axes[1].fill_between(t, 0, gz, where=gz > 0, color='forestgreen', alpha=0.7)
axes[1].fill_between(t, 0, gz, where=gz < 0, color='tomato', alpha=0.7)
axes[1].axhline(0, color='k', linewidth=0.8)
axes[1].set_ylabel('Gz', fontsize=9)
axes[1].yaxis.set_tick_params(labelleft=False)

# Gy
gy = rect(t, 0.13, 0.06, 0.55)
axes[2].fill_between(t, 0, gy, color='darkorange', alpha=0.7)
axes[2].axhline(0, color='k', linewidth=0.8)
axes[2].set_ylabel('Gy', fontsize=9)
axes[2].yaxis.set_tick_params(labelleft=False)

# Gx: pre-phaser (negative) then readout (positive) — echo at TE
pre_end = TE - 0.08
gx = rect(t, 0.13, 0.06, -0.5) + rect(t, TE - 0.08, 0.16, 0.5)
# Spoiler at end
gx += rect(t, TR - 0.10, 0.08, 0.6)
axes[3].fill_between(t, 0, gx, where=gx > 0, color='purple', alpha=0.7)
axes[3].fill_between(t, 0, gx, where=gx < 0, color='plum', alpha=0.7)
# ADC
axes[3].fill_between(t, -0.15, -0.05,
                     where=((t >= TE - 0.08) & (t <= TE + 0.08)),
                     color='gold', alpha=0.9, label='ADC')
axes[3].axhline(0, color='k', linewidth=0.8)
axes[3].set_ylabel('Gx / ADC', fontsize=9)
axes[3].yaxis.set_tick_params(labelleft=False)
axes[3].annotate('Spoiler', (TR - 0.06, 0.65), ha='center', fontsize=8,
                 color='purple')

for ax in axes:
    ax.axvline(TE, color='gray', linestyle=':', linewidth=0.8, alpha=0.6)

axes[3].annotate('', xy=(TE, -0.25), xytext=(0, -0.25),
                 arrowprops=dict(arrowstyle='<->', color='k'))
axes[3].text(TE / 2, -0.35, 'TE', ha='center', fontsize=9)
axes[3].annotate('', xy=(TR, -0.25), xytext=(0, -0.25),
                 arrowprops=dict(arrowstyle='<->', color='dimgray'))
axes[3].text(TR / 2, -0.45, 'TR', ha='center', fontsize=9, color='dimgray')

axes[3].set_xlim(0, TR)
axes[3].set_xlabel('Time', fontsize=9)
axes[3].set_ylim(-0.55, 0.85)
axes[0].set_title('Gradient-Echo (Spoiled GRE) Pulse Sequence Diagram',
                  fontsize=11, pad=8)
plt.tight_layout()
plt.show()
```

## Mathematical Model

### Ernst Equation

When RF spoiling is applied (residual $M_{xy}$ is randomised each TR), the
longitudinal magnetisation reaches a steady state described by the **Ernst
equation** {cite}`ernst1966`:

$$
S(\alpha, \text{TR}, \text{TE}) =
  M_0 \sin\alpha \;\frac{1 - E_1}{1 - E_1 \cos\alpha}\;
  e^{-\text{TE}/T_2^*}
$$ (eq-ernst)

where $E_1 = e^{-\text{TR}/T_1}$.

### Ernst Angle

The **Ernst angle** $\alpha_E$ maximises the steady-state signal for a given
TR and T1 {cite}`ernst1966`:

$$
\cos\alpha_E = e^{-\text{TR}/T_1} = E_1
$$ (eq-ernst-angle)

For short TR ($\text{TR} \ll T_1$), $E_1 \to 1$ and $\alpha_E \to 0°$
(very small flip angle). For long TR ($\text{TR} \gg T_1$), $E_1 \to 0$ and
$\alpha_E \to 90°$, recovering the fully-relaxed case.

### T2\* Decay

Unlike the spin echo, the GRE echo decays with T2\*:

$$
\frac{1}{T_2^*} = \frac{1}{T_2} + \gamma\,\Delta B_0
$$

where $\Delta B_0$ represents local field inhomogeneity. At air–tissue
interfaces and near iron deposits, T2\* can be much shorter than T2,
causing signal voids in GRE images.

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

fig = make_subplots(rows=1, cols=2,
                    subplot_titles=['Signal vs Flip Angle (Ernst curve)',
                                    'Signal vs TR'])

# Panel 1: Ernst curve for different T1
alpha_arr = np.linspace(0.5, 90, 500)
alpha_rad = np.radians(alpha_arr)

T1_vals = [400, 840, 1300, 2000]
colors = ['#e41a1c', '#1f77b4', '#ff7f0e', '#2ca02c']
TR_fixed = 30  # ms
TE_fixed = 5

for T1, col in zip(T1_vals, colors):
    E1 = np.exp(-TR_fixed / T1)
    S = np.sin(alpha_rad) * (1 - E1) / (1 - E1 * np.cos(alpha_rad))
    ernst_angle = np.degrees(np.arccos(E1))
    fig.add_trace(go.Scatter(
        x=alpha_arr, y=S,
        name=f'T1 = {T1} ms (αE = {ernst_angle:.0f}°)',
        line=dict(color=col, width=2)
    ), row=1, col=1)

# Panel 2: Signal vs TR at Ernst angle
TR_arr = np.linspace(5, 500, 500)
for T1, col in zip(T1_vals, colors):
    E1 = np.exp(-TR_arr / T1)
    alpha_ernst = np.arccos(E1)
    S = np.sin(alpha_ernst) * (1 - E1) / (1 - E1 * np.cos(alpha_ernst))
    S *= np.exp(-TE_fixed / 50)  # T2* = 50 ms
    fig.add_trace(go.Scatter(
        x=TR_arr, y=S,
        name=f'T1 = {T1} ms',
        line=dict(color=col, width=2),
        showlegend=False
    ), row=1, col=2)

fig.update_xaxes(title_text='Flip Angle α (degrees)', row=1, col=1)
fig.update_xaxes(title_text='TR (ms)', row=1, col=2)
fig.update_yaxes(title_text='Normalised Signal (a.u.)', row=1, col=1)

fig.update_layout(
    height=420,
    title_text='Gradient-Echo Steady-State Signal (Spoiled GRE)',
    legend=dict(orientation='h', yanchor='bottom', y=-0.30,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows the Ernst curve for several T1 values (TR = 30 ms). The
peak of each curve identifies the Ernst angle; shorter-T1 tissues have a
larger Ernst angle at this TR. The right panel shows how the optimal
(Ernst-angle) signal grows with TR as more T1 recovery accumulates.

## Data Fitting

The Ernst equation can be rearranged to estimate T1 from a set of GRE
acquisitions at different flip angles (Variable Flip Angle — VFA). This
linearisation is covered in depth in Chapter 3 (T1 Mapping); here we
demonstrate a direct non-linear fit.

```{code-cell} python
:tags: [hide-input]

import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objects as go

rng = np.random.default_rng(7)

T1_true = 1000.0   # ms
M0_true = 1.0
TR = 30.0          # ms
TE = 5.0           # ms
T2s = 40.0         # ms (T2*)

alpha_data = np.array([2, 4, 6, 8, 10, 15, 20, 30, 40, 60, 80])
noise = 0.005

def ernst(alpha_deg, M0, T1):
    E1 = np.exp(-TR / T1)
    alpha = np.radians(alpha_deg)
    return M0 * np.sin(alpha) * (1 - E1) / (1 - E1 * np.cos(alpha)) \
           * np.exp(-TE / T2s)

signal_true  = ernst(alpha_data, M0_true, T1_true)
signal_noisy = signal_true + rng.normal(0, noise, size=len(alpha_data))
signal_noisy = np.clip(signal_noisy, 0, None)

popt, pcov = curve_fit(ernst, alpha_data, signal_noisy,
                       p0=[0.9, 800], bounds=([0, 50], [2, 5000]))
perr = np.sqrt(np.diag(pcov))
M0_fit, T1_fit = popt

alpha_dense = np.linspace(0.5, 85, 400)

fig = go.Figure()
fig.add_trace(go.Scatter(x=alpha_dense, y=ernst(alpha_dense, *popt),
                         name=f'Fit: T1 = {T1_fit:.0f} ± {perr[1]:.0f} ms',
                         line=dict(color='royalblue', width=2.5)))
fig.add_trace(go.Scatter(x=alpha_data, y=signal_noisy,
                         mode='markers',
                         name=f'Synthetic data (T1_true = {T1_true:.0f} ms)',
                         marker=dict(color='tomato', size=9)))
fig.add_trace(go.Scatter(x=alpha_dense, y=ernst(alpha_dense, M0_true, T1_true),
                         name='Ground truth',
                         line=dict(color='gray', dash='dot', width=1.5)))

# Mark Ernst angle
alpha_E = np.degrees(np.arccos(np.exp(-TR / T1_fit)))
S_E = ernst(alpha_E, *popt)
fig.add_trace(go.Scatter(x=[alpha_E], y=[S_E], mode='markers',
                         marker=dict(color='gold', size=12, symbol='star'),
                         name=f'Ernst angle = {alpha_E:.0f}°'))

fig.update_layout(
    title=f'GRE T1 Fitting via Variable Flip Angle (TR = {TR} ms)',
    xaxis_title='Flip Angle α (degrees)',
    yaxis_title='Signal (a.u.)',
    template='plotly_white',
    height=380,
    legend=dict(x=0.98, y=0.98, xanchor='right', yanchor='top'),
)
fig.show()
```

A non-linear least-squares fit to the Ernst equation recovers T1 from the
synthetic VFA dataset. The gold star marks the fitted Ernst angle. In
practice, B1 inhomogeneity (Chapter 2) can bias VFA T1 estimates, which
is why B1 mapping is an important prerequisite for reliable T1 mapping.

## References

```{bibliography}
:filter: docname in docnames
```
