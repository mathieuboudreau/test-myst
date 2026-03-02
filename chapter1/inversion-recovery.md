---
title: "1.3 Inversion Recovery"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 1.3 Inversion Recovery

Inversion recovery (IR) is one of the most powerful contrast mechanisms in
MRI. By prepending a 180° inversion pulse to any readout module, it encodes
the longitudinal relaxation time T1 directly into signal magnitude, allowing
selective nulling of specific tissues and quantitative T1 mapping.

## Pulse Sequence

The IR sequence consists of three events:

1. **180° inversion pulse** — inverts the longitudinal magnetisation:
   $M_z(0^+) = -M_0$.
2. **Inversion time TI** — $M_z$ recovers toward $M_0$ via T1 relaxation.
3. **90° readout pulse** — samples $M_z(\text{TI})$, then TE and TR complete
   the cycle.

```{code-cell} python
---
tags: [hide-input]
mystnb:
  figure:
    name: fig-ir-sequence
    caption: |
      Inversion recovery pulse sequence. A 180° inversion pulse at $t = 0$ is
      followed by the inversion time TI, then a 90° readout with standard
      spin-echo or GRE encoding. The signal depends on how much $M_z$ has
      recovered during TI.
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
TI = 0.3
TE_loc = TI + 0.35   # readout at TI + short TE

def rect(t, t0, width, height=1.0):
    r = np.zeros_like(t)
    r[(t >= t0) & (t <= t0 + width)] = height
    return r

def sinc_pulse(t, t0, width, height=1.0):
    x = (t - (t0 + width / 2)) / (width / 6)
    v = height * np.sinc(x)
    v[(t < t0) | (t > t0 + width)] = 0
    return v

# RF: 180° inversion then 90° readout
rf = sinc_pulse(t, 0.02, 0.07, 1.0) + sinc_pulse(t, TI, 0.06, 0.55)
axes[0].fill_between(t, 0, rf, color='royalblue', alpha=0.8)
axes[0].axhline(0, color='k', linewidth=0.8)
axes[0].annotate('180°', (0.055, 1.05), ha='center', fontsize=9)
axes[0].annotate('90°', (TI + 0.03, 0.60), ha='center', fontsize=9)
axes[0].set_ylabel('RF', fontsize=9)
axes[0].set_ylim(-0.2, 1.3)
axes[0].yaxis.set_tick_params(labelleft=False)

# Gz
gz = rect(t, 0.01, 0.09, 0.7) + rect(t, 0.10, 0.04, -0.4) + \
     rect(t, TI - 0.01, 0.08, 0.7) + rect(t, TI + 0.07, 0.04, -0.4)
axes[1].fill_between(t, 0, gz, where=gz > 0, color='forestgreen', alpha=0.7)
axes[1].fill_between(t, 0, gz, where=gz < 0, color='tomato', alpha=0.7)
axes[1].axhline(0, color='k', linewidth=0.8)
axes[1].set_ylabel('Gz', fontsize=9)
axes[1].yaxis.set_tick_params(labelleft=False)

# Gy
gy = rect(t, TI + 0.12, 0.05, 0.55)
axes[2].fill_between(t, 0, gy, color='darkorange', alpha=0.7)
axes[2].axhline(0, color='k', linewidth=0.8)
axes[2].set_ylabel('Gy', fontsize=9)
axes[2].yaxis.set_tick_params(labelleft=False)

# Gx
gx = rect(t, TI + 0.12, 0.05, -0.5) + rect(t, TE_loc - 0.07, 0.14, 0.5)
axes[3].fill_between(t, 0, gx, where=gx > 0, color='purple', alpha=0.7)
axes[3].fill_between(t, 0, gx, where=gx < 0, color='plum', alpha=0.7)
axes[3].fill_between(t, -0.15, -0.05,
                     where=((t >= TE_loc - 0.07) & (t <= TE_loc + 0.07)),
                     color='gold', alpha=0.9)
axes[3].axhline(0, color='k', linewidth=0.8)
axes[3].set_ylabel('Gx / ADC', fontsize=9)
axes[3].yaxis.set_tick_params(labelleft=False)

for ax in axes:
    ax.axvline(TI, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)

axes[3].annotate('', xy=(TI, -0.25), xytext=(0, -0.25),
                 arrowprops=dict(arrowstyle='<->', color='k'))
axes[3].text(TI / 2, -0.35, 'TI', ha='center', fontsize=9)
axes[3].annotate('', xy=(TR, -0.25), xytext=(0, -0.25),
                 arrowprops=dict(arrowstyle='<->', color='dimgray'))
axes[3].text(TR / 2, -0.45, 'TR', ha='center', fontsize=9, color='dimgray')

axes[3].set_xlim(0, TR)
axes[3].set_xlabel('Time', fontsize=9)
axes[3].set_ylim(-0.55, 0.8)
axes[0].set_title('Inversion Recovery Pulse Sequence Diagram', fontsize=11, pad=8)
plt.tight_layout()
plt.show()
```

## Mathematical Model

### Signal Equation

After the 180° inversion, $M_z$ recovers toward $M_0$:

$$
M_z(\text{TI}) = M_0 \left(1 - 2e^{-\text{TI}/T_1} + e^{-\text{TR}/T_1}\right)
$$ (eq-ir-mz)

For the common approximation $\text{TR} \gg T_1$:

$$
M_z(\text{TI}) \approx M_0 \left(1 - 2e^{-\text{TI}/T_1}\right)
$$

Because magnitude images are typically reconstructed, the detected signal is:

$$
S(\text{TI}) = \left|M_z(\text{TI})\right| \cdot e^{-\text{TE}/T_2}
$$ (eq-ir-signal)

### Tissue Nulling

The zero-crossing occurs at:

$$
\text{TI}_\text{null} = T_1 \ln 2 \approx 0.693\,T_1
$$

This is exploited clinically:

| Application | Suppressed tissue | Typical TI at 1.5 T |
|-------------|-------------------|---------------------|
| **STIR** (Short TI IR) | Fat ($T_1 \approx 270$ ms) | $\approx 187$ ms |
| **FLAIR** (Fluid-Attenuated IR) | CSF ($T_1 \approx 4000$ ms) | $\approx 2800$ ms |

### Three-Parameter Model

When TR is not sufficiently long, the general three-parameter model should
be used for T1 fitting:

$$
S(\text{TI}) = a \left|1 - b\,e^{-\text{TI}/T_1}\right|
$$ (eq-ir-3param)

where $a = M_0\,e^{-\text{TE}/T_2}$ and $b \leq 2$ accounts for incomplete
inversion and T1 saturation ($b = 1 + e^{-\text{TR}/T_1}$ for ideal inversion).

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

tissues = {
    'White Matter': {'T1': 840,  'color': '#1f77b4'},
    'Gray Matter':  {'T1': 1300, 'color': '#ff7f0e'},
    'Fat':          {'T1': 370,  'color': '#d62728'},
    'CSF':          {'T1': 4500, 'color': '#2ca02c'},
}

fig = make_subplots(
    rows=1, cols=2,
    subplot_titles=['Mz Recovery (TR = 6000 ms)',
                    'Magnitude Signal |Mz| with Null Points'],
)

TI_arr = np.linspace(10, 6000, 1000)
TR = 6000

for name, p in tissues.items():
    T1 = p['T1']
    col = p['color']
    TI_null = T1 * np.log(2)

    Mz = 1 - 2 * np.exp(-TI_arr / T1) + np.exp(-TR / T1)
    MagMz = np.abs(Mz)

    fig.add_trace(go.Scatter(x=TI_arr, y=Mz,
                             name=name, line=dict(color=col, width=2),
                             legendgroup=name), row=1, col=1)
    fig.add_trace(go.Scatter(x=TI_arr, y=MagMz,
                             name=name, line=dict(color=col, width=2),
                             legendgroup=name, showlegend=False), row=1, col=2)
    # Null point marker
    Mz_null = 1 - 2 * np.exp(-TI_null / T1) + np.exp(-TR / T1)
    fig.add_trace(go.Scatter(x=[TI_null], y=[0],
                             mode='markers',
                             marker=dict(color=col, size=10, symbol='x'),
                             showlegend=False,
                             name=f'Null {name}'), row=1, col=2)

fig.add_hline(y=0, line_dash='dash', line_color='gray',
              line_width=1, row=1, col=1)

fig.update_xaxes(title_text='Inversion Time TI (ms)', row=1, col=1)
fig.update_xaxes(title_text='Inversion Time TI (ms)', row=1, col=2)
fig.update_yaxes(title_text='M_z / M_0', row=1, col=1)
fig.update_yaxes(title_text='|M_z| / M_0', row=1, col=2)

fig.update_layout(
    height=420,
    title_text='Inversion Recovery Signal (TR = 6000 ms, TE ≈ 0)',
    legend=dict(orientation='h', yanchor='bottom', y=-0.28,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows the signed longitudinal magnetisation recovery; tissues
with short T1 (fat) recover first, crossing zero at short TI. The right panel
shows the magnitude signal: the null point (×) for each tissue corresponds to
$\text{TI}_\text{null} = T_1 \ln 2$.

## Data Fitting

T1 is estimated by fitting [Eq. %s](eq-ir-3param) to multi-TI data. Using
magnitude images requires careful treatment of the sign ambiguity — phase
images or a sign-restored approach produce more accurate fits.

```{code-cell} python
:tags: [hide-input]

import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objects as go

rng = np.random.default_rng(99)

T1_true = 1000.0
TR = 6000.0
noise = 0.025

TI_data = np.array([50, 100, 200, 350, 500, 750, 1000,
                    1500, 2000, 3000, 4000, 5500])

# Three-parameter model (signed Mz — phase-corrected)
def ir_model_signed(TI, a, b, T1):
    return a * (1 - b * np.exp(-TI / T1))

b_true = 1 + np.exp(-TR / T1_true)   # ≈ 2 for long TR
signal_true  = ir_model_signed(TI_data, 1.0, b_true, T1_true)
signal_noisy = signal_true + rng.normal(0, noise, size=len(TI_data))

popt, pcov = curve_fit(ir_model_signed, TI_data, signal_noisy,
                       p0=[1.0, 1.95, 800],
                       bounds=([0, 1, 50], [2, 2, 8000]))
perr = np.sqrt(np.diag(pcov))
a_fit, b_fit, T1_fit = popt

TI_dense = np.linspace(10, 6000, 600)

fig = go.Figure()
fig.add_trace(go.Scatter(x=TI_dense,
                         y=ir_model_signed(TI_dense, *popt),
                         name=f'Fit: T1 = {T1_fit:.0f} ± {perr[2]:.0f} ms',
                         line=dict(color='royalblue', width=2.5)))
fig.add_trace(go.Scatter(x=TI_data, y=signal_noisy,
                         mode='markers',
                         name=f'Synthetic data (T1_true = {T1_true:.0f} ms)',
                         marker=dict(color='tomato', size=9)))
fig.add_trace(go.Scatter(x=TI_dense,
                         y=ir_model_signed(TI_dense, 1.0, b_true, T1_true),
                         name='Ground truth',
                         line=dict(color='gray', dash='dot', width=1.5)))
fig.add_hline(y=0, line_dash='dot', line_color='gray', line_width=1)
# Null point
TI_null_fit = T1_fit * np.log(b_fit)
fig.add_vline(x=TI_null_fit, line_dash='dash', line_color='orange',
              annotation_text=f'TI_null = {TI_null_fit:.0f} ms',
              annotation_position='top right')

fig.update_layout(
    title=f'Inversion Recovery T1 Fitting (TR = {TR:.0f} ms)',
    xaxis_title='Inversion Time TI (ms)',
    yaxis_title='Signal (phase-corrected, a.u.)',
    template='plotly_white',
    height=380,
    legend=dict(x=0.98, y=0.02, xanchor='right', yanchor='bottom'),
)
fig.show()
```

The three-parameter fit recovers both T1 and the inversion efficiency $b$.
The fitted null point (orange dashed line) agrees well with
$T_1 \ln b / 1 \approx T_1 \ln 2$ when the inversion is ideal. In practice,
incomplete RF inversion (e.g., at high field due to B1 inhomogeneity) reduces
$b$ below 2, which is why the three-parameter model is preferred over the
simpler two-parameter ($b = 2$) form.

## References

```{bibliography}
:filter: docname in docnames
```
