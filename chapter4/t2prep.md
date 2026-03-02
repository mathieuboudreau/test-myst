---
title: "4.3 T2 Preparation"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 4.3 T2 Preparation

T2 preparation (T2prep) is a magnetisation preparation module that imparts
T2-dependent contrast before any imaging readout {cite}`brittain1995`.
Rather than measuring T2 from a multi-echo spin-echo train, T2prep converts
the T2 decay time into initial longitudinal magnetisation that a subsequent
readout block images. This decoupling of T2 encoding from image acquisition
makes T2prep compatible with fast 3D readouts (bSSFP, EPI, GRE), high-field
imaging, and cardiac applications where scan time is limited.

## Pulse Sequence

The T2prep module consists of three parts:

1. **90°ₓ tip-down** — rotates $M_z$ into the transverse plane.
2. **CPMG refocusing train** — $N$ adiabatic 180°ᵧ pulses with total
   duration $\tau = \text{T2prep}$ allow T2 dephasing to accumulate
   while refocusing B0 and B1 inhomogeneity.
3. **90°₋ₓ tip-up** — returns the T2-decayed transverse magnetisation
   back to the longitudinal axis.

A spoiler gradient then crushes residual transverse magnetisation. The
prepared $M_z$ is:

$$
M_z^{\text{prep}} = M_0\,e^{-\tau/T_2}
$$

A conventional imaging readout (GRE, bSSFP, EPI) then encodes this
prepared magnetisation. By repeating the experiment with different T2prep
durations $\tau$, T2 is mapped.

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

t_90x = 0.02; w_90x = 0.04; t_90mx = 0.34; w_90mx = 0.04
t_180_list = [0.10, 0.18, 0.26]; w_180 = 0.04; t_sp = 0.40
t_rd_start = 0.52; n_rd = 4; rd_gap = 0.10; rd_w = 0.04

rf_t2p = sinc_pulse(t, t_90x, w_90x, 0.44)
gz_t2p = rect(t, t_90x-0.01, w_90x+0.02, 0.65) + rect(t, t_90x+w_90x+0.01, 0.025, -0.32)
gy_t2p = np.zeros_like(t)
gx_t2p = rect(t, t_sp, 0.08, 0.60)
adc_t2p = np.zeros(len(t), dtype=bool)

for t180 in t_180_list:
    rf_t2p  += hard_pulse(t, t180, w_180, 0.88)
    gz_t2p  += rect(t, t180-0.01, w_180+0.02, 0.65)

rf_t2p += sinc_pulse(t, t_90mx, w_90mx, 0.44)
gz_t2p += rect(t, t_90mx-0.01, w_90mx+0.02, 0.65)

for k in range(n_rd):
    tp = t_rd_start + k * rd_gap
    rf_t2p  += sinc_pulse(t, tp, rd_w, 0.35)
    gz_t2p  += rect(t, tp-0.005, rd_w+0.01, 0.65) + rect(t, tp+rd_w+0.005, 0.02, -0.32)
    gy_t2p  += rect(t, tp+rd_w+0.01, 0.02, 0.32 + 0.08*k)
    gx_t2p  += rect(t, tp+rd_w+0.01, 0.025, -0.38) + rect(t, tp+rd_w+0.04, 0.06, 0.38)
    adc_t2p |= (t >= tp+rd_w+0.04) & (t <= tp+rd_w+0.10)

axes[0].fill_between(t, 0, sinc_pulse(t, t_90x, w_90x, 0.44),
                     color='cornflowerblue', alpha=0.90, label='90°ₓ')
for t180 in t_180_list:
    axes[0].fill_between(t, 0, hard_pulse(t, t180, w_180, 0.88),
                         color='tomato', alpha=0.82,
                         label='180°ᵧ' if t180 == t_180_list[0] else None)
axes[0].fill_between(t, 0, sinc_pulse(t, t_90mx, w_90mx, 0.44),
                     color='navy', alpha=0.80, label='90°₋ₓ')
for k in range(n_rd):
    tp = t_rd_start + k * rd_gap
    axes[0].fill_between(t, 0, sinc_pulse(t, tp, rd_w, 0.35),
                         color='mediumseagreen', alpha=0.80,
                         label='α readout' if k == 0 else None)
axes[0].axhline(0, color='k', lw=0.7)
axes[0].set_ylabel('RF', fontsize=9); axes[0].yaxis.set_tick_params(labelleft=False)
axes[0].set_ylim(-0.1, 1.05); axes[0].legend(loc='upper right', fontsize=7, ncol=2)
axes[0].set_title('T2-Prepared Sequence — T2prep Module + GRE Readout', fontsize=11, pad=5)

axes[1].fill_between(t, 0, gz_t2p, where=gz_t2p>0, color='forestgreen', alpha=0.7)
axes[1].fill_between(t, 0, gz_t2p, where=gz_t2p<0, color='tomato', alpha=0.7)
axes[1].axhline(0, color='k', lw=0.7); axes[1].set_ylabel('Gz', fontsize=9)
axes[1].yaxis.set_tick_params(labelleft=False)

axes[2].fill_between(t, 0, gy_t2p, color='darkorange', alpha=0.7)
axes[2].axhline(0, color='k', lw=0.7); axes[2].set_ylabel('Gy', fontsize=9)
axes[2].yaxis.set_tick_params(labelleft=False)

axes[3].fill_between(t, 0, gx_t2p, where=gx_t2p>0, color='purple', alpha=0.7)
axes[3].fill_between(t, 0, gx_t2p, where=gx_t2p<0, color='plum', alpha=0.7)
axes[3].fill_between(t, -0.18, -0.06, where=adc_t2p, color='gold', alpha=0.9)
axes[3].axhline(0, color='k', lw=0.7)
axes[3].set_ylabel('Gx / ADC', fontsize=9); axes[3].yaxis.set_tick_params(labelleft=False)
axes[3].set_xlabel('Time', fontsize=9); axes[3].set_xlim(0, 1); axes[3].set_ylim(-0.28, 0.75)

axes[3].annotate('', xy=(t_sp+0.06, -0.23), xytext=(t_90x, -0.23),
                 arrowprops=dict(arrowstyle='<->', color='k', lw=1.2))
axes[3].text((t_90x + t_sp+0.06)/2, -0.27, 'T2prep (τ)',
             ha='center', fontsize=9)
axes[3].axvline(t_rd_start - 0.01, color='gray', lw=0.8, ls='--')
axes[3].text(t_rd_start + 0.15, 0.62, 'Readout block', ha='center', fontsize=8, color='gray')
axes[3].text(0.22, 0.62, 'T2prep module', ha='center', fontsize=8)
plt.close()
```

:::{figure} t2prep_pulse_sequence.png
:name: fig-t2prep-sequence
:align: center
T2prep pulse sequence. The T2prep module (left of dashed line) consists
of a 90°ₓ tip-down, three adiabatic 180°ᵧ refocusing pulses, a 90°₋ₓ
tip-up, and a spoiler. The T2-weighted longitudinal magnetisation is
then imaged by a GRE readout block. Repeating with different durations
τ maps T2.
:::

## Mathematical Model

### Prepared Magnetisation

Immediately after the T2prep module, the stored longitudinal magnetisation is:

$$
M_z^{\text{prep}}(\tau) = M_0\,e^{-\tau/T_2}
$$ (eq-t2prep-mz)

This ignores T1 recovery during the module (valid when τ $\ll T_1$).
Including the finite T1 recovery during the preparation period gives:

$$
M_z^{\text{prep}}(\tau) \approx M_0\left[e^{-\tau/T_2}
  + (1 - e^{-\tau/T_1})\right] \approx M_0\,e^{-\tau/T_2}
$$

for τ $\ll T_1$, confirming that [Eq. %s](eq-t2prep-mz) is accurate
for T2prep durations below ~200 ms.

### Signal Model

With a GRE readout of flip angle $\alpha$ immediately after preparation:

$$
S(\tau) = S_0\,e^{-\tau/T_2}, \qquad S_0 = M_0\sin\alpha
$$ (eq-t2prep-signal)

Acquiring $K$ images at different T2prep durations
$\tau_1 < \tau_2 < \cdots < \tau_K$ and fitting [Eq. %s](eq-t2prep-signal)
yields T2. Typical T2prep durations for brain tissue at 3 T are
0, 24, 48, 72, 96 ms.

### Contrast

The T2prep contrast mechanism is directly analogous to inversion recovery
for T1 (Chapter 3). The null condition — T2prep duration where a tissue
with T2 = T2₀ is nulled — is:

$$
\tau_\text{null} = T_2\,\ln\!\left(\frac{1}{f}\right)
$$

where $f$ is the fraction of equilibrium signal retained by the readout
(typically $f \approx 1$ for a small flip angle). This selective nulling
can suppress fat or specific tissue types.

### B1 Robustness

The adiabatic 180° pulses in the T2prep CPMG train are highly insensitive
to B1 inhomogeneity, provided the adiabaticity condition is met
{cite}`brittain1995`. This makes T2prep more robust than simple SE or MESE
at 3 T and 7 T where B1 varies substantially.

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

tau_arr = np.linspace(0, 200, 500)   # ms

def t2prep_signal(tau, T2, S0=1.0):
    return S0 * np.exp(-tau / T2)

fig = make_subplots(rows=1, cols=2,
                    subplot_titles=['T2prep Signal vs Duration by Tissue',
                                    'T2-Dependent Contrast at Fixed T2prep'])

tissues = [('WM  T2=75 ms',  75,  '#1f77b4'),
           ('GM  T2=100 ms', 100, '#d62728'),
           ('CSF T2=1500 ms',1500,'#2ca02c')]

for name, T2, col in tissues:
    fig.add_trace(go.Scatter(x=tau_arr, y=t2prep_signal(tau_arr, T2),
                             name=name, line=dict(color=col, width=2.2)), row=1, col=1)

# Right: contrast between WM and GM as function of T2prep duration
S_wm = t2prep_signal(tau_arr, 75)
S_gm = t2prep_signal(tau_arr, 100)
CNR  = np.abs(S_gm - S_wm) / 0.05   # normalised to noise σ=0.05

# Optimal T2prep
idx_opt = np.argmax(CNR)
tau_opt = tau_arr[idx_opt]

fig.add_trace(go.Scatter(x=tau_arr, y=S_wm, name='WM signal',
                         line=dict(color='#1f77b4', width=2)), row=1, col=2)
fig.add_trace(go.Scatter(x=tau_arr, y=S_gm, name='GM signal',
                         line=dict(color='#d62728', width=2)), row=1, col=2)
fig.add_trace(go.Scatter(x=tau_arr, y=(S_gm - S_wm),
                         name='GM − WM difference',
                         line=dict(color='purple', dash='dash', width=2)), row=1, col=2)
fig.add_vline(x=tau_opt, line_dash='dot', line_color='purple',
              annotation_text=f'τ_opt={tau_opt:.0f} ms',
              annotation_position='top right', row=1, col=2)

for col_idx in [1, 2]:
    fig.update_xaxes(title_text='T2prep duration τ (ms)', row=1, col=col_idx)
fig.update_yaxes(title_text='Signal (a.u.)', row=1, col=1)
fig.update_yaxes(title_text='Signal (a.u.)', row=1, col=2)
fig.add_hline(y=0, line_dash='dot', line_color='lightgray', row=1, col=2)

fig.update_layout(
    height=430,
    title_text='T2prep Simulation',
    legend=dict(orientation='h', yanchor='bottom', y=-0.28,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows T2prep signal decay for three tissues. CSF decays
much more slowly than parenchyma, providing high T2-weighted contrast at
long τ. The right panel shows the WM–GM signal difference (dashed purple)
which peaks near τ ≈ 87 ms; this optimal T2prep duration maximises
WM/GM contrast.

## Data Fitting

```{code-cell} python
:tags: [hide-input]

import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objects as go

rng = np.random.default_rng(17)

# Typical cardiac/brain T2prep protocol
tau_pts = np.array([0, 24, 48, 72, 96])   # ms
T2_true = 90.0; S0_true = 1.0; noise = 0.030

S_true  = S0_true * np.exp(-tau_pts / T2_true)
S_noisy = S_true + rng.normal(0, noise, len(tau_pts))
S_noisy = np.maximum(S_noisy, 0)

def t2prep_model(tau, S0, T2):
    return S0 * np.exp(-tau / T2)

popt, pcov = curve_fit(t2prep_model, tau_pts, S_noisy, p0=[1.0, 70.0],
                       bounds=([0, 10], [5, 3000]), maxfev=4000)
S0_fit, T2_fit = popt; T2_err = np.sqrt(np.diag(pcov))[1]

tau_dense = np.linspace(0, 120, 300)

fig = go.Figure()
fig.add_trace(go.Scatter(x=tau_dense, y=S0_true * np.exp(-tau_dense/T2_true),
                         name='Ground truth', mode='lines',
                         line=dict(color='gray', dash='dot', width=1.5)))
fig.add_trace(go.Scatter(x=tau_pts, y=S_noisy, name='Noisy data',
                         mode='markers', marker=dict(color='royalblue', size=10)))
fig.add_trace(go.Scatter(x=tau_dense, y=t2prep_model(tau_dense, *popt),
                         name=f'NLS fit  T2={T2_fit:.0f} ± {T2_err:.0f} ms',
                         line=dict(color='tomato', width=2.2)))

fig.update_layout(
    title=f'T2prep T2 Fitting  (True T2 = {T2_true:.0f} ms, 5 acquisitions)',
    xaxis_title='T2prep duration τ (ms)',
    yaxis_title='Signal (a.u.)',
    template='plotly_white', height=380,
    legend=dict(x=0.55, y=0.95),
)
fig.show()
```

T2prep achieves T2 mapping with only 5 acquisitions, much fewer than
the 8–16 echoes typically used in MESE. The sparsely sampled T2prep
protocol is well suited to cardiac MRI where breath-holds limit scan
time. The constraint that B1-insensitive adiabatic pulses are used in
the module means T2prep results are more robust to transmit field
inhomogeneity than MESE at high field strength.

## References

```{bibliography}
:filter: docname in docnames
```
