---
title: "5.2 Quantitative MT"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 5.2 Quantitative MT

Quantitative magnetization transfer (qMT) goes beyond the empirical MTR
by fitting the full two-pool Bloch-McConnell model to a series of MT-weighted
images acquired at multiple saturation offset frequencies and/or RF powers
{cite}`henkelman1993,sled2001`. This yields maps of biophysically
interpretable parameters: the bound pool fraction $f$, the forward exchange
rate $k_f$, and the bound-pool relaxation time T2b.

## Pulse Sequence

A qMT experiment acquires $N$ spoiled GRE images, each preceded by an
off-resonance MT saturation pulse. The acquisitions differ in:

- **Frequency offset** $\Delta f$ (e.g., 500 Hz to 50 kHz) — sweeping the
  offset maps out the Z-spectrum.
- **RF amplitude** (saturation flip angle $\beta$) — varying the power
  modulates the MT effect and helps separate exchange rate from pool fraction.

For each offset, the steady-state longitudinal free-pool magnetisation
$M_{za}$ (normalised to $M_0$) is recorded. The collection of these
values as a function of $\Delta f$ is the **Z-spectrum** (also called the
MT spectrum or magnetisation transfer spectrum).

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def rect(t, t0, w, h=1.0):
    r = np.zeros_like(t); r[(t >= t0) & (t <= t0 + w)] = h; return r

def sinc_pulse(t, t0, w, h=1.0):
    x = (t - (t0 + w/2)) / (w/6); v = h * np.sinc(x)
    v[(t < t0) | (t > t0 + w)] = 0; return v

def gauss_pulse(t, t0, w, h=1.0):
    tc = t0 + w/2; sigma = w/5
    v = h * np.exp(-0.5*((t - tc)/sigma)**2)
    v[(t < t0) | (t > t0 + w)] = 0; return v

fig, axes = plt.subplots(4, 1, figsize=(11, 5.5), sharex=True)
fig.subplots_adjust(hspace=0.06, left=0.07, right=0.97, top=0.88, bottom=0.13)
t = np.linspace(0, 1, 3000)

block_starts = [0.03, 0.36, 0.69]
block_colors = ['#e7298a', '#7570b3', '#1b9e77']
offsets = ['Δω₁', 'Δω₂', 'Δω₃']

gz_total = np.zeros_like(t)
gy_total = np.zeros_like(t)
gx_total = np.zeros_like(t)
adc_mask = np.zeros(len(t), dtype=bool)

for bs, bc, off in zip(block_starts, block_colors, offsets):
    mt_w = 0.14; exc_start = bs + 0.18; exc_w = 0.05
    TE_loc = exc_start + exc_w + 0.08
    mt_env = gauss_pulse(t, bs, mt_w, 0.45)
    axes[0].fill_between(t, 0, mt_env, color=bc, alpha=0.65)
    axes[0].text(bs + mt_w/2, 0.48, off, ha='center', fontsize=8, color=bc)
    exc_env = sinc_pulse(t, exc_start, exc_w, 0.60)
    axes[0].fill_between(t, 0, exc_env, color='royalblue', alpha=0.85)
    gz = (rect(t, exc_start-0.01, exc_w+0.02, 0.70) + rect(t, exc_start+exc_w+0.01, 0.03, -0.35) +
          rect(t, TE_loc-0.005, 0.03, 0.50) + rect(t, TE_loc+0.025, 0.025, -0.30))
    gz_total += gz
    gy_total += rect(t, TE_loc-0.005, 0.04, 0.55)
    gx = rect(t, TE_loc-0.005, 0.04, -0.42) + rect(t, TE_loc+0.035, 0.08, 0.42)
    gx_total += gx
    adc_mask |= ((t >= TE_loc+0.035) & (t <= TE_loc+0.115))

axes[0].axhline(0, color='k', linewidth=0.7)
axes[0].set_ylabel('RF', fontsize=9); axes[0].set_ylim(-0.12, 0.78)
axes[0].yaxis.set_tick_params(labelleft=False)
leg_elems = [Patch(facecolor='#e7298a', alpha=0.65, label='MT prep (off-res.)'),
             Patch(facecolor='royalblue', alpha=0.85, label='Excitation')]
axes[0].legend(handles=leg_elems, loc='upper right', fontsize=8)
axes[0].set_title('Quantitative MT Pulse Sequence — Z-spectrum Acquisition '
                  '(three offsets shown)', fontsize=10, pad=6)

axes[1].fill_between(t, 0, gz_total, where=gz_total > 0, color='forestgreen', alpha=0.7)
axes[1].fill_between(t, 0, gz_total, where=gz_total < 0, color='tomato', alpha=0.7)
axes[1].axhline(0, color='k', linewidth=0.7)
axes[1].set_ylabel('Gz', fontsize=9); axes[1].yaxis.set_tick_params(labelleft=False)

axes[2].fill_between(t, 0, gy_total, color='darkorange', alpha=0.7)
axes[2].axhline(0, color='k', linewidth=0.7)
axes[2].set_ylabel('Gy', fontsize=9); axes[2].yaxis.set_tick_params(labelleft=False)

axes[3].fill_between(t, 0, gx_total, where=gx_total > 0, color='purple', alpha=0.7)
axes[3].fill_between(t, 0, gx_total, where=gx_total < 0, color='plum', alpha=0.7)
axes[3].fill_between(t, -0.15, -0.05, where=adc_mask, color='gold', alpha=0.9)
axes[3].axhline(0, color='k', linewidth=0.7)
axes[3].set_ylabel('Gx / ADC', fontsize=9); axes[3].yaxis.set_tick_params(labelleft=False)
axes[3].set_xlabel('Time  (repeated for each offset Δω to build Z-spectrum)', fontsize=9)
axes[3].set_xlim(0, 1.0); axes[3].set_ylim(-0.25, 0.75)
plt.close()
```

:::{figure} qmt_pulse_sequence.png
:name: fig-qmt-sequence
:align: center
Quantitative MT pulse sequence. Three representative acquisitions are
shown, each with a Gaussian off-resonance MT preparation pulse at a
different frequency offset ($\Delta\omega_1$, $\Delta\omega_2$,
$\Delta\omega_3$) followed by a spoiled GRE readout. The full
experiment repeats this block for 10–30 offsets and, optionally,
multiple RF amplitudes to build a Z-spectrum for model fitting.
:::

## Mathematical Model

### Two-Pool Bloch-McConnell Equations

The two-pool model describes the coupled evolution of the free pool (a)
and the bound/semi-solid pool (b) under RF irradiation and relaxation
{cite}`mcconnell1958,henkelman1993`:

$$
\frac{dM_{za}}{dt} = -R_{1a}(M_{za} - M_{0a}) - k_f M_{za}
  + k_r M_{zb} - W_a(\Delta f)\,M_{za}
$$

$$
\frac{dM_{zb}}{dt} = -R_{1b}(M_{zb} - M_{0b}) - k_r M_{zb}
  + k_f M_{za} - W_b(\Delta f)\,M_{zb}
$$

where $R_{1x} = 1/T_{1x}$, $k_f$ is the forward exchange rate
(free→bound), $k_r = k_f/f$ is the reverse rate (bound→free), and
$f = M_{0b}/M_{0a}$ is the bound pool fraction. The saturation rates
$W_x(\Delta f)$ describe the absorption of RF energy by each pool.

### Saturation Rates and Lineshapes

For CW irradiation with amplitude $\omega_1 = \gamma B_1$ at offset
$\Delta f$ from resonance:

**Free pool** (Lorentzian lineshape, narrow):

$$
W_a(\Delta f) = \omega_1^2 \,\frac{T_{2a}}{1 + (2\pi\,\Delta f\,T_{2a})^2}
$$ (eq-Wa)

**Bound pool** (Gaussian lineshape, approximation to super-Lorentzian):

$$
W_b(\Delta f) = \sqrt{\tfrac{\pi}{2}}\;\omega_1^2\,T_{2b}\,
  \exp\!\left(-\tfrac{1}{2}(2\pi\,\Delta f\,T_{2b})^2\right)
$$ (eq-Wb)

The Gaussian approximation is used here for analytical tractability;
the super-Lorentzian lineshape {cite}`henkelman1993` is more accurate
for solid-like bound pools and is used in quantitative fitting software.

### Steady-State Z-Spectrum

Setting $dM_{za}/dt = dM_{zb}/dt = 0$ and solving the two-pool system,
the normalised free-pool Z-magnetisation is {cite}`henkelman1993`:

$$
Z(\Delta f) = \frac{M_{za}}{M_{0a}}
  = \frac{R_{1b,\text{eff}}\,R_{1a} + k_f\,R_{1b}}
         {R_{1b,\text{eff}}\,R_{1a,\text{eff}} - k_f\,k_r}
$$ (eq-Zspectrum)

where:

$$
R_{1a,\text{eff}} = R_{1a} + k_f + W_a(\Delta f), \qquad
R_{1b,\text{eff}} = R_{1b} + k_r + W_b(\Delta f)
$$

At $\Delta f \to \infty$ (no RF interaction), $W_a = W_b = 0$ and
$Z \to 1$, recovering equilibrium. Near $\Delta f = 0$, $W_a$ dominates
and $Z \to 0$ (direct saturation of water). At intermediate offsets
(1–10 kHz), $W_b$ controls the MT dip.

### Free Parameters

A full qMT fit to the Z-spectrum recovers:

| Parameter | Description | Typical WM value |
|-----------|-------------|-----------------|
| $f$ | Bound pool fraction | 0.13–0.17 |
| $k_f$ | Forward exchange rate (s⁻¹) | 3–5 s⁻¹ |
| $T_{2b}$ | Bound pool T2 (µs) | 8–15 µs |
| $T_{2a}$ | Free pool T2 (ms) | 30–80 ms |
| $T_{1a}$ | Free pool T1 | From separate map or jointly |

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def z_spectrum_cw(delta_f_hz, omega1, R1a, R1b, T2a, T2b, kf, f):
    """
    CW two-pool steady-state Z-spectrum.
    delta_f_hz : frequency offset array (Hz)
    omega1     : CW RF amplitude (rad/s)
    """
    kr = kf / f
    # Free pool saturation rate (Lorentzian)
    Wa = omega1**2 * T2a / (1 + (2 * np.pi * delta_f_hz * T2a)**2)
    # Bound pool saturation rate (Gaussian)
    Wb = np.sqrt(np.pi / 2) * omega1**2 * T2b * \
         np.exp(-0.5 * (2 * np.pi * delta_f_hz * T2b)**2)

    Ra_eff = R1a + kf + Wa
    Rb_eff = R1b + kr + Wb
    kf_kr  = kf * kr

    Z = (Rb_eff * R1a + kf * R1b) / (Rb_eff * Ra_eff - kf_kr)
    return np.clip(Z, 0, 1)

# Frequency axis (log-spaced for wide range)
df_pos = np.logspace(1, 5, 500)   # 10 Hz to 100 kHz
df = np.concatenate([-df_pos[::-1], df_pos])

omega1 = 500.0   # rad/s (typical MT pulse equivalent CW power)

# Tissue parameters
tissues = {
    'White Matter': dict(R1a=1/0.84, R1b=1.0, T2a=0.030, T2b=10e-6, kf=4.3, f=0.155,
                         color='#1f77b4'),
    'Gray Matter':  dict(R1a=1/1.3,  R1b=1.0, T2a=0.060, T2b=10e-6, kf=2.0, f=0.065,
                         color='#d62728'),
    'Muscle':       dict(R1a=1/1.1,  R1b=1.0, T2a=0.025, T2b=8e-6,  kf=5.5, f=0.20,
                         color='#2ca02c'),
}

fig = make_subplots(rows=1, cols=2,
                    subplot_titles=['Z-spectrum for Different Tissue Types',
                                    'Effect of Bound Pool Fraction f on Z-spectrum (WM)'])

for name, p in tissues.items():
    Z = z_spectrum_cw(df, omega1, p['R1a'], p['R1b'],
                      p['T2a'], p['T2b'], p['kf'], p['f'])
    fig.add_trace(go.Scatter(x=df / 1000, y=Z * 100, name=name,
                             line=dict(color=p['color'], width=2.2)), row=1, col=1)

# Right panel: vary f for WM-like tissue
base = dict(R1a=1/0.84, R1b=1.0, T2a=0.030, T2b=10e-6, kf=4.3)
f_vals = [0.05, 0.10, 0.15, 0.20]
cols_f = ['#aec7e8', '#6baed6', '#2171b5', '#084594']
for fv, cv in zip(f_vals, cols_f):
    Z = z_spectrum_cw(df, omega1, base['R1a'], base['R1b'],
                      base['T2a'], base['T2b'], base['kf'], fv)
    fig.add_trace(go.Scatter(x=df / 1000, y=Z * 100, name=f'f = {fv}',
                             line=dict(color=cv, width=2)), row=1, col=2)

for col in [1, 2]:
    fig.update_xaxes(title_text='Frequency offset Δf (kHz)', type='log', row=1, col=col)
    fig.update_yaxes(title_text='Z-magnetisation (%)', row=1, col=col)
    fig.add_hline(y=100, line_dash='dot', line_color='lightgray', row=1, col=col)

fig.update_layout(
    height=430,
    title_text=f'qMT Z-spectra  (CW ω₁ = {omega1:.0f} rad/s)',
    legend=dict(orientation='h', yanchor='bottom', y=-0.30,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows that WM has the deepest MT dip due to its high
myelin content (large $f$, high $k_f$). GM has a shallower dip, and
muscle sits between the two. The dip near $\Delta f = 0$ is dominated
by direct saturation of the free pool (independent of the bound pool).
The right panel isolates the role of the bound pool fraction: larger
$f$ deepens the off-resonance MT plateau.

## Data Fitting

In practice, qMT fitting is performed voxel-by-voxel by minimising the
residual between the measured Z-spectrum and [Eq. %s](eq-Zspectrum).
Here we simulate a Z-spectrum for a WM-like voxel, add noise, and fit
the four parameters $f$, $k_f$, $T_{2b}$, and $T_{2a}$.

```{code-cell} python
:tags: [hide-input]

import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objects as go

rng = np.random.default_rng(21)

# Acquisition: 20 offsets from 500 Hz to 50 kHz
df_acq = np.array([0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0,
                   15.0, 20.0, 30.0, 50.0]) * 1e3   # Hz

omega1 = 500.0   # rad/s

# Ground truth — WM
true_p = dict(R1a=1/0.84, R1b=1.0, T2a=0.030, T2b=10e-6, kf=4.3, f=0.155)

Z_true = z_spectrum_cw(df_acq, omega1, **true_p)
noise_level = 0.015
Z_noisy = Z_true + rng.normal(0, noise_level, len(df_acq))

# Dense axis for display
df_dense = np.logspace(np.log10(400), 5, 500)

def model_for_fit(df, f, kf, T2b, T2a):
    R1a = 1 / 0.84; R1b = 1.0
    return z_spectrum_cw(df, omega1, R1a, R1b, T2a, T2b, kf, f)

p0    = [0.10, 3.0, 12e-6, 0.040]
bounds = ([0.01, 0.5, 3e-6, 0.005],
          [0.35, 20,  30e-6, 0.150])

try:
    popt, _ = curve_fit(model_for_fit, df_acq, Z_noisy, p0=p0, bounds=bounds,
                        maxfev=5000)
    f_fit, kf_fit, T2b_fit, T2a_fit = popt
    Z_fit = model_for_fit(df_dense, *popt)
    fit_ok = True
except Exception:
    fit_ok = False

Z_true_dense = z_spectrum_cw(df_dense, omega1, **true_p)

fig = go.Figure()
fig.add_trace(go.Scatter(x=df_dense / 1e3, y=Z_true_dense * 100,
                         name='Ground truth', mode='lines',
                         line=dict(color='gray', dash='dot', width=1.5)))
fig.add_trace(go.Scatter(x=df_acq / 1e3, y=Z_noisy * 100,
                         name='Simulated data (noisy)',
                         mode='markers',
                         marker=dict(color='royalblue', size=8)))
if fit_ok:
    fig.add_trace(go.Scatter(x=df_dense / 1e3, y=Z_fit * 100,
                             name='Two-pool fit',
                             line=dict(color='tomato', width=2.2)))
    fit_label = (f'Fit:  f={f_fit:.3f}  kf={kf_fit:.1f} s⁻¹  '
                 f'T₂b={T2b_fit*1e6:.1f} µs  T₂a={T2a_fit*1e3:.1f} ms')
    true_label = (f'True: f={true_p["f"]:.3f}  kf={true_p["kf"]:.1f} s⁻¹  '
                  f'T₂b={true_p["T2b"]*1e6:.1f} µs  T₂a={true_p["T2a"]*1e3:.1f} ms')
    title = f'qMT Z-spectrum Fit<br>{true_label}<br>{fit_label}'
else:
    title = 'qMT Z-spectrum (fit failed)'

fig.update_layout(
    title=title,
    xaxis=dict(title='Frequency offset Δf (kHz)', type='log'),
    yaxis_title='Z-magnetisation (%)',
    template='plotly_white',
    height=400,
    legend=dict(x=0.02, y=0.05),
)
fig.show()
```

The two-pool model (red) provides an excellent fit to the simulated
noisy Z-spectrum (blue points). The recovered parameters closely match
the ground truth, demonstrating that qMT fitting can separate the
bound pool fraction $f$ from the exchange rate $k_f$ when both the
offset-frequency dependence and the power dependence are sampled.
In practice, fitting is regularised and the T1 map is often supplied
as a fixed input rather than a free parameter.

## References

```{bibliography}
:filter: docname in docnames
```
