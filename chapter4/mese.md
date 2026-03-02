---
title: "4.1 Multi-Echo Spin Echo"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 4.1 Multi-Echo Spin Echo

The multi-echo spin echo (MESE) sequence is the reference standard for
T2 mapping. It applies a single 90° excitation followed by a train of
180° refocusing pulses (the Carr-Purcell-Meiboom-Gill, or CPMG, echo
train {cite}`carr1954,meiboom1958`), acquiring one k-space line per echo.
Because each 180° pulse refocuses static $\Delta B_0$ dephasing, the
signal at each echo time TE reflects true spin-spin relaxation T2.

## Pulse Sequence

The MESE/CPMG sequence consists of:

1. **90° slice-selective excitation** — tips $M_z$ into the transverse plane.
2. **180° refocusing pulses** — applied at intervals TE/2, creating
   spin echoes at TE, 2TE, 3TE, …, $N \cdot \text{TE}$.
3. **Readout gradients** — a frequency-encode gradient surrounds each
   echo; the phase-encode step is constant (all echoes sample the same
   k-space line).

CPMG phase cycling (90°ₓ – 180°ᵧ) ensures that imperfections in the
refocusing pulse accumulate coherently rather than destroying signal
{cite}`meiboom1958`.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def rect(t, t0, w, h=1.0):
    r = np.zeros_like(t); r[(t >= t0) & (t <= t0 + w)] = h; return r

def sinc_pulse(t, t0, w, h=1.0):
    x = (t - (t0 + w/2)) / (w/6); v = h * np.sinc(x)
    v[(t < t0) | (t > t0 + w)] = 0; return v

def hard_pulse(t, t0, w, h=1.0):
    r = np.zeros_like(t); r[(t >= t0) & (t <= t0 + w)] = h; return r

fig, axes = plt.subplots(4, 1, figsize=(12, 5.2), sharex=True)
fig.subplots_adjust(hspace=0.06, left=0.07, right=0.97, top=0.87, bottom=0.14)
t = np.linspace(0, 1, 5000)

N_echoes = 3
exc_t = 0.04; exc_w = 0.05
ref_w = 0.045; TE_spacing = 0.22

ref_times  = [exc_t + exc_w/2 + TE_spacing/2 * (2*k+1) for k in range(N_echoes)]
echo_times = [exc_t + exc_w/2 + TE_spacing * (k+1)     for k in range(N_echoes)]

gz_row = rect(t, exc_t-0.01, exc_w+0.02, 0.65) + rect(t, exc_t+exc_w+0.005, 0.03, -0.32)
gy_row = np.zeros_like(t); gx_row = np.zeros_like(t); adc = np.zeros(len(t), dtype=bool)

for k, (rt, et) in enumerate(zip(ref_times, echo_times)):
    gz_row  += rect(t, rt-ref_w/2-0.015, ref_w+0.03, 0.65)
    gy_row  += rect(t, et-0.04, 0.025, 0.32 + 0.12*k)
    gx_row  += rect(t, et-0.055, 0.03, -0.42) + rect(t, et-0.02, 0.08, 0.42)
    adc     |= (t >= et-0.02) & (t <= et+0.06)

axes[0].fill_between(t, 0, sinc_pulse(t, exc_t, exc_w, 0.45),
                     color='cornflowerblue', alpha=0.9, label='90° excitation')
for k, rt in enumerate(ref_times):
    axes[0].fill_between(t, 0, hard_pulse(t, rt-ref_w/2, ref_w, 0.88),
                         color='tomato', alpha=0.82,
                         label='180° refocus' if k == 0 else None)
axes[0].axhline(0, color='k', lw=0.7)
axes[0].set_ylabel('RF', fontsize=9); axes[0].yaxis.set_tick_params(labelleft=False)
axes[0].set_ylim(-0.1, 1.05); axes[0].legend(loc='upper right', fontsize=7.5)
axes[0].set_title('Multi-Echo Spin Echo (MESE/CPMG) Pulse Sequence', fontsize=11, pad=5)

axes[1].fill_between(t, 0, gz_row, where=gz_row>0, color='forestgreen', alpha=0.7)
axes[1].fill_between(t, 0, gz_row, where=gz_row<0, color='tomato', alpha=0.7)
axes[1].axhline(0, color='k', lw=0.7)
axes[1].set_ylabel('Gz', fontsize=9); axes[1].yaxis.set_tick_params(labelleft=False)

axes[2].fill_between(t, 0, gy_row, color='darkorange', alpha=0.7)
axes[2].axhline(0, color='k', lw=0.7)
axes[2].set_ylabel('Gy', fontsize=9); axes[2].yaxis.set_tick_params(labelleft=False)

axes[3].fill_between(t, 0, gx_row, where=gx_row>0, color='purple', alpha=0.7)
axes[3].fill_between(t, 0, gx_row, where=gx_row<0, color='plum', alpha=0.7)
axes[3].fill_between(t, -0.18, -0.06, where=adc, color='gold', alpha=0.9)
axes[3].axhline(0, color='k', lw=0.7)
axes[3].set_ylabel('Gx / ADC', fontsize=9); axes[3].yaxis.set_tick_params(labelleft=False)
axes[3].set_xlabel('Time', fontsize=9); axes[3].set_xlim(0, 1); axes[3].set_ylim(-0.28, 0.68)

exc_centre = exc_t + exc_w/2
for k, et in enumerate(echo_times):
    axes[3].annotate('', xy=(et, -0.23), xytext=(exc_centre, -0.23),
                     arrowprops=dict(arrowstyle='<->', color=f'C{k}', lw=1.0))
    axes[3].text((exc_centre + et)/2, -0.27, f'TE{"₁₂₃"[k]}',
                 ha='center', fontsize=8.5, color=f'C{k}')
plt.close()
```

:::{figure} mese_pulse_sequence.png
:name: fig-mese-sequence
:align: center
MESE/CPMG pulse sequence showing a 90° excitation followed by three
180° refocusing pulses. Spin echoes form at TE₁, 2TE₁, 3TE₁.
Each echo acquires one k-space line; the full 2D dataset is built by
repeating the sequence with different phase-encode amplitudes.
:::

## Mathematical Model

### Monoexponential Decay

For a single-compartment tissue, the signal magnitude at echo $n$ is:

$$
S(n \cdot \text{TE}) = S_0\,e^{-n \cdot \text{TE}/T_2}
$$ (eq-mese-mono)

where $S_0 = M_0\sin\alpha$ is the initial transverse magnetisation and
TE is the inter-echo spacing. Fitting [Eq. %s](eq-mese-mono) to data
acquired at $N$ echo times yields T2 and S0.

The logarithm linearises the model:

$$
\ln S(n \cdot \text{TE}) = \ln S_0 - \frac{n \cdot \text{TE}}{T_2}
$$

so T2 is the negative reciprocal of the slope of $\ln S$ vs TE.

### Biexponential (Multi-Component) Model

Many tissues contain multiple water pools with distinct T2 values. For
two components:

$$
S(\text{TE}) = S_{\text{fast}}\,e^{-\text{TE}/T_{2,\text{fast}}}
             + S_{\text{slow}}\,e^{-\text{TE}/T_{2,\text{slow}}}
$$ (eq-mese-biexp)

In white matter, the fast component (T2 ≈ 10–20 ms) represents myelin
water; the slow component (T2 ≈ 70–80 ms) represents intra/extra-cellular
water. The ratio $S_\text{fast}/(S_\text{fast}+S_\text{slow})$ is the
**myelin water fraction (MWF)**, a marker of myelin integrity
{cite}`mackay1994`.

### Stimulated Echoes and EPG

Real MESE acquisitions suffer from stimulated echo contamination when
the 180° pulses are imperfect. The **extended phase graph (EPG)**
algorithm {cite}`weigel2015` models the full echo-train signal
including stimulated echoes, enabling accurate T2 estimation even with
non-ideal refocusing.

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

TE_arr = np.linspace(10, 300, 500)   # ms

def mese_mono(TE, T2, S0=1.0):
    return S0 * np.exp(-TE / T2)

def mese_biexp(TE, T2f, T2s, ff, S0=1.0):
    return S0 * (ff * np.exp(-TE / T2f) + (1 - ff) * np.exp(-TE / T2s))

fig = make_subplots(rows=1, cols=2,
                    subplot_titles=['Monoexponential T2 Decay by Tissue',
                                    'Biexponential Decay — Myelin Water vs Bulk Water'])

tissues = [('WM  T2=75 ms',  75,  '#1f77b4'),
           ('GM  T2=100 ms', 100, '#d62728'),
           ('CSF T2=1500 ms',1500,'#2ca02c')]

for name, T2, col in tissues:
    fig.add_trace(go.Scatter(x=TE_arr, y=mese_mono(TE_arr, T2),
                             name=name, line=dict(color=col, width=2.2)), row=1, col=1)

# Right: biexponential WM (myelin water + IE water)
TE_biexp = np.linspace(5, 200, 500)
# WM typical parameters
for label, T2f, T2s, ff, col in [
        ('WM full signal', 15, 75, 0.15, '#1f77b4'),
        ('Myelin water component', 15, 75, 0.15, 'tomato'),
        ('IE water component', 15, 75, 0.0, 'steelblue')]:
    if label == 'WM full signal':
        y = mese_biexp(TE_biexp, T2f, T2s, ff)
    elif label == 'Myelin water component':
        y = 1.0 * 0.15 * np.exp(-TE_biexp / 15)
    else:
        y = 1.0 * 0.85 * np.exp(-TE_biexp / 75)
    dash = 'solid' if 'full' in label else 'dash'
    fig.add_trace(go.Scatter(x=TE_biexp, y=y, name=label,
                             line=dict(color=col, dash=dash, width=2)), row=1, col=2)

for col in [1, 2]:
    fig.update_xaxes(title_text='TE (ms)', row=1, col=col)
fig.update_yaxes(title_text='Signal (a.u.)', row=1, col=1)
fig.update_yaxes(title_text='Signal (a.u.)', row=1, col=2)

fig.update_layout(
    height=430,
    title_text='MESE T2 Decay Simulation',
    legend=dict(orientation='h', yanchor='bottom', y=-0.30,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows monoexponential T2 decay for WM, GM, and CSF.
The large T2 difference between CSF and parenchyma is the basis for
T2-weighted clinical imaging. The right panel decomposes the WM signal
into its myelin water (fast, dashed) and intra/extra-cellular (slow,
dashed) components; the solid blue line is the measured sum.

## Data Fitting

```{code-cell} python
:tags: [hide-input]

import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objects as go

rng = np.random.default_rng(7)

# 8-echo MESE at 15 ms spacing
TE_pts = np.arange(1, 9) * 15.0   # ms: 15, 30, 45, ..., 120
T2_true = 80.0; S0_true = 1.0; noise = 0.025

S_true  = S0_true * np.exp(-TE_pts / T2_true)
S_noisy = S_true + rng.normal(0, noise, len(TE_pts))
S_noisy = np.maximum(S_noisy, 0)

# 1. Log-linear fit
log_S = np.log(np.maximum(S_noisy, 1e-6))
coeffs = np.polyfit(TE_pts, log_S, 1)
T2_linear = -1.0 / coeffs[0]
S0_linear = np.exp(coeffs[1])

# 2. Non-linear least-squares
def mono_model(TE, S0, T2):
    return S0 * np.exp(-TE / T2)

popt, pcov = curve_fit(mono_model, TE_pts, S_noisy, p0=[1.0, 60.0],
                       bounds=([0, 5], [5, 2000]), maxfev=4000)
S0_nl, T2_nl = popt; T2_err = np.sqrt(np.diag(pcov))[1]

TE_dense = np.linspace(5, 130, 300)

fig = go.Figure()
fig.add_trace(go.Scatter(x=TE_dense, y=S0_true * np.exp(-TE_dense / T2_true),
                         name='Ground truth', mode='lines',
                         line=dict(color='gray', dash='dot', width=1.5)))
fig.add_trace(go.Scatter(x=TE_pts, y=S_noisy, name='Noisy data',
                         mode='markers', marker=dict(color='royalblue', size=9)))
fig.add_trace(go.Scatter(x=TE_dense, y=S0_linear * np.exp(-TE_dense / T2_linear),
                         name=f'Log-linear fit  T2={T2_linear:.0f} ms',
                         line=dict(color='darkorange', dash='dash', width=2)))
fig.add_trace(go.Scatter(x=TE_dense, y=mono_model(TE_dense, *popt),
                         name=f'NLS fit  T2={T2_nl:.0f} ± {T2_err:.0f} ms',
                         line=dict(color='tomato', width=2.2)))

fig.update_layout(
    title=f'MESE T2 Fitting  (True T2 = {T2_true:.0f} ms)',
    xaxis_title='TE (ms)', yaxis_title='Signal (a.u.)',
    template='plotly_white', height=380,
    legend=dict(x=0.55, y=0.95),
)
fig.show()
```

Both the log-linear and non-linear least-squares fits recover T2
accurately. The log-linear approach is faster but biased toward early
echoes (high SNR); NLS weights all echoes equally and returns a
confidence interval. For multi-component fitting, NLS is required.

## References

```{bibliography}
:filter: docname in docnames
```
