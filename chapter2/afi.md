---
title: "2.2 Actual Flip-angle Imaging"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 2.2 Actual Flip-angle Imaging

Actual flip-angle imaging (AFI) is a rapid 3D B1 mapping technique introduced
by Yarnykh in 2007 {cite}`yarnykh2007`. By interleaving two steady-state
signals with different repetition times (TR1 and TR2) within a single
sequence, AFI achieves B1 maps with reduced sensitivity to T1 and shorter
scan times compared to the double angle method.

## Pulse Sequence

AFI applies a single RF pulse of flip angle $\alpha$ and then collects two
signals after two different delays — TR1 and TR2 — before the next RF pulse:

1. **RF pulse** (flip angle $\alpha$) — excites the spin system.
2. **TR1 acquisition** — signal $S_1$ is acquired after a short delay TR1.
3. **TR2 delay** — the magnetisation continues to relax for time TR2.
4. **Next RF pulse** — the cycle repeats.

RF spoiling is applied between excitations to avoid coherent transverse
contributions. The ratio $n = \text{TR2}/\text{TR1}$ is typically set to
$n = 5$–10.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(4, 1, figsize=(11, 5.5), sharex=True)
fig.subplots_adjust(hspace=0.05)

t = np.linspace(0, 1, 2000)
TR1 = 0.2
TR2 = 0.8   # TR2/TR1 = 4 (n=4 for illustration)
TR_total = TR1 + TR2

def rect(t, t0, w, h=1.0):
    r = np.zeros_like(t)
    r[(t >= t0) & (t <= t0 + w)] = h
    return r

def sinc_pulse(t, t0, w, h=1.0):
    x = (t - (t0 + w/2)) / (w/6)
    v = h * np.sinc(x)
    v[(t < t0) | (t > t0 + w)] = 0
    return v

# Two RF pulses per cycle (representative)
pulse_times = [0.01, TR_total + 0.01]
rf = sinc_pulse(t, 0.01, 0.06, 0.65) + sinc_pulse(t, TR_total + 0.01, 0.06, 0.65)
axes[0].fill_between(t, 0, rf, color='royalblue', alpha=0.8)
axes[0].axhline(0, color='k', linewidth=0.8)
for pt in pulse_times:
    axes[0].annotate('α', (pt + 0.03, 0.70), ha='center', fontsize=9)
axes[0].set_ylabel('RF', fontsize=9)
axes[0].set_ylim(-0.15, 0.95)
axes[0].yaxis.set_tick_params(labelleft=False)

# Gz
gz = rect(t, 0.005, 0.08, 0.7) + rect(t, 0.085, 0.04, -0.35) + \
     rect(t, TR_total + 0.005, 0.08, 0.7) + rect(t, TR_total + 0.085, 0.04, -0.35)
axes[1].fill_between(t, 0, gz, where=gz > 0, color='forestgreen', alpha=0.7)
axes[1].fill_between(t, 0, gz, where=gz < 0, color='tomato', alpha=0.7)
axes[1].axhline(0, color='k', linewidth=0.8)
axes[1].set_ylabel('Gz', fontsize=9)
axes[1].yaxis.set_tick_params(labelleft=False)

# Gy
gy = rect(t, 0.13, 0.05, 0.55) + rect(t, TR_total + 0.13, 0.05, 0.55)
axes[2].fill_between(t, 0, gy, color='darkorange', alpha=0.7)
axes[2].axhline(0, color='k', linewidth=0.8)
axes[2].set_ylabel('Gy', fontsize=9)
axes[2].yaxis.set_tick_params(labelleft=False)

# Gx + ADC for S1 and S2 in first cycle
TE1 = 0.17; TE2_loc = TR1 + 0.07

gx = (rect(t, 0.13, 0.04, -0.4) + rect(t, TE1 - 0.04, 0.08, 0.4) +
      rect(t, TR1 + 0.02, 0.04, -0.4) + rect(t, TE2_loc - 0.04, 0.08, 0.4) +
      rect(t, TR_total + 0.13, 0.04, -0.4) +
      rect(t, TR_total + TE1 - 0.04, 0.08, 0.4))
axes[3].fill_between(t, 0, gx, where=gx > 0, color='purple', alpha=0.7)
axes[3].fill_between(t, 0, gx, where=gx < 0, color='plum', alpha=0.7)
# ADC windows
for (ta, tb, label) in [(TE1 - 0.04, TE1 + 0.04, 'S₁'),
                         (TE2_loc - 0.04, TE2_loc + 0.04, 'S₂')]:
    axes[3].fill_between(t, -0.15, -0.05,
                         where=((t >= ta) & (t <= tb)), color='gold', alpha=0.9)
    axes[3].text((ta + tb) / 2, -0.23, label, ha='center', fontsize=8)
axes[3].axhline(0, color='k', linewidth=0.8)
axes[3].set_ylabel('Gx / ADC', fontsize=9)
axes[3].yaxis.set_tick_params(labelleft=False)

# TR1 and TR2 brackets
axes[3].annotate('', xy=(TR1, -0.35), xytext=(0, -0.35),
                 arrowprops=dict(arrowstyle='<->', color='royalblue', lw=1.2))
axes[3].text(TR1 / 2, -0.43, 'TR₁', ha='center', fontsize=9, color='royalblue')
axes[3].annotate('', xy=(TR_total, -0.35), xytext=(TR1, -0.35),
                 arrowprops=dict(arrowstyle='<->', color='tomato', lw=1.2))
axes[3].text((TR1 + TR_total) / 2, -0.43, 'TR₂', ha='center',
             fontsize=9, color='tomato')

axes[3].set_xlim(0, 2 * TR_total)
axes[3].set_xlabel('Time', fontsize=9)
axes[3].set_ylim(-0.52, 0.75)
axes[0].set_title(
    'AFI Pulse Sequence (one TR cycle shown; n = TR₂/TR₁ = 4)', fontsize=11, pad=8)
plt.tight_layout()
plt.close()
```

:::{figure} afi_pulse_sequence.png
:name: fig-afi-sequence
:align: center
AFI pulse sequence. Each cycle consists of one RF pulse followed by two
acquisition windows: $S_1$ acquired after TR1 and $S_2$ acquired after TR2.
The ratio $n = \text{TR2}/\text{TR1}$ controls the sensitivity and T1
independence of the method.
:::

## Mathematical Model

### Steady-State Signals

Under RF spoiling, the two steady-state signals after the $\alpha$-pulse are
{cite}`yarnykh2007`:

$$
S_1 = M_0 \sin\alpha \;\frac{1 - E_2 + E_2(1-E_1)\cos\alpha}
                             {1 - E_1 E_2 \cos^2\alpha}
$$ (eq-afi-s1)

$$
S_2 = M_0 \sin\alpha \;\frac{1 - E_1 + E_1(1-E_2)\cos\alpha}
                             {1 - E_1 E_2 \cos^2\alpha}
$$ (eq-afi-s2)

where $E_1 = e^{-\text{TR1}/T_1}$ and $E_2 = e^{-\text{TR2}/T_1}$.

### B1 Estimation

The ratio $r = S_1/S_2$ simplifies considerably when $\text{TR1} \ll T_1$
and $n = \text{TR2}/\text{TR1} \gg 1$. In the general case, the actual flip
angle is recovered by solving {cite}`yarnykh2007`:

$$
r = \frac{S_1}{S_2} \approx \frac{1 + (n-1)e^{-\text{TR1}/T_1}}
                                   {n + (1-n)e^{-\text{TR1}/T_1}}
$$ (eq-afi-r)

Note that in this approximation the ratio $r$ depends on T1 and TR1 but
**not** on $M_0$ or T2\*. Inverting:

$$
\alpha_\text{actual} = \arccos\!\left(\frac{r\,n - 1}{n - r}\right)
$$ (eq-afi-alpha)

The B1 scaling factor is $F = \alpha_\text{actual} / \alpha_\text{nominal}$.

:::{note}
The T1-independence of the ratio $r$ (Eq. [%s](eq-afi-r)) holds only
approximately. For finite $n$ and TR1 comparable to T1, a residual T1
dependence remains. Choosing $n \geq 5$ and TR1 $\leq 0.5\,T_1$ gives
B1 errors below 2% for typical brain tissues {cite}`yarnykh2007`.
:::

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

alpha_nom = np.radians(60)   # nominal flip angle
F_arr = np.linspace(0.4, 1.6, 300)
alpha_arr = F_arr * alpha_nom

def afi_ratio(alpha, TR1, TR2, T1):
    E1 = np.exp(-TR1 / T1)
    E2 = np.exp(-TR2 / T1)
    denom = 1 - E1 * E2 * np.cos(alpha)**2
    S1 = np.sin(alpha) * (1 - E2 + E2*(1 - E1)*np.cos(alpha)) / denom
    S2 = np.sin(alpha) * (1 - E1 + E1*(1 - E2)*np.cos(alpha)) / denom
    return S1 / S2

def afi_recover_alpha(r, n):
    arg = np.clip((r * n - 1) / (n - r), -1, 1)
    return np.arccos(arg)

fig = make_subplots(rows=1, cols=2,
                    subplot_titles=['Signal Ratio vs True B1 Factor',
                                    'B1 Map Error vs T1 (TR1=20 ms, n=5)'])

TR1 = 20     # ms
n_vals = [3, 5, 8]
T1_ref = 1000

colors_n = ['#d62728', '#ff7f0e', '#2ca02c']
for n, col in zip(n_vals, colors_n):
    TR2 = n * TR1
    r = afi_ratio(alpha_arr, TR1, TR2, T1_ref)
    alpha_rec = afi_recover_alpha(r, n)
    F_rec = alpha_rec / alpha_nom
    fig.add_trace(go.Scatter(x=F_arr, y=r,
                             name=f'n = {n}',
                             line=dict(color=col, width=2)), row=1, col=1)

# T1 sensitivity at fixed n=5
T1_vals = np.linspace(200, 3000, 300)
n_fixed = 5
TR2_fixed = n_fixed * TR1
F_test = 1.0   # test at true B1 = 1
alpha_test = F_test * alpha_nom

pct_errors = []
for T1 in T1_vals:
    r_T1 = afi_ratio(alpha_test, TR1, TR2_fixed, T1)
    F_est = afi_recover_alpha(r_T1, n_fixed) / alpha_nom
    pct_errors.append((F_est - F_test) / F_test * 100)

fig.add_trace(go.Scatter(x=T1_vals, y=pct_errors,
                         name='n = 5, F_true = 1.0',
                         line=dict(color='royalblue', width=2.5),
                         showlegend=True), row=1, col=2)

# DAM comparison (same TR, long TR approximation)
TR_dam = TR1 + TR2_fixed  # same total time as one AFI cycle
E1_dam_arr = np.exp(-TR_dam / T1_vals)
S_alpha_dam  = np.sin(alpha_test) * (1 - E1_dam_arr) / (1 - E1_dam_arr * np.cos(alpha_test))
S_2alpha_dam = np.sin(2*alpha_test) * (1 - E1_dam_arr) / (1 - E1_dam_arr * np.cos(2*alpha_test))
ratio_dam = S_2alpha_dam / S_alpha_dam
F_dam = np.arccos(ratio_dam / 2) / alpha_nom
pct_err_dam = (F_dam - F_test) / F_test * 100

fig.add_trace(go.Scatter(x=T1_vals, y=pct_err_dam,
                         name=f'DAM (TR = {TR_dam} ms)',
                         line=dict(color='tomato', width=2, dash='dash')),
             row=1, col=2)

fig.add_hline(y=0, line_dash='dot', line_color='gray', row=1, col=2)

fig.update_xaxes(title_text='True B1 Factor F', row=1, col=1)
fig.update_xaxes(title_text='T1 (ms)', row=1, col=2)
fig.update_yaxes(title_text='Signal Ratio S₁/S₂', row=1, col=1)
fig.update_yaxes(title_text='B1 Map Error (%)', row=1, col=2)

fig.update_layout(
    height=420,
    title_text=f'AFI Simulation (α_nom = 60°, TR1 = {TR1} ms)',
    legend=dict(orientation='h', yanchor='bottom', y=-0.3,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows how the AFI signal ratio varies with the true B1 factor
for different values of $n$. Larger $n$ increases the dynamic range of the
ratio and improves sensitivity. The right panel compares the T1-induced bias
of AFI ($n = 5$, TR1 = 20 ms) against a DAM acquisition with the same total
TR: AFI shows dramatically reduced T1 sensitivity across all tissue T1 values.

## Data Fitting

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go

rng = np.random.default_rng(42)

N = 200
alpha_nom = np.radians(60)
TR1 = 20.0    # ms
n = 5
TR2 = n * TR1
T1_tissue = 1000.0   # ms (constant for simplicity)
noise = 0.012

x = np.linspace(0, np.pi, N)
F_true = 0.65 + 0.45 * np.cos(x - np.pi/3)**2

def afi_signals(alpha_act, TR1, TR2, T1):
    E1 = np.exp(-TR1 / T1)
    E2 = np.exp(-TR2 / T1)
    denom = 1 - E1 * E2 * np.cos(alpha_act)**2
    S1 = np.sin(alpha_act) * (1 - E2 + E2*(1-E1)*np.cos(alpha_act)) / denom
    S2 = np.sin(alpha_act) * (1 - E1 + E1*(1-E2)*np.cos(alpha_act)) / denom
    return S1, S2

alpha_true = F_true * alpha_nom
S1_true, S2_true = afi_signals(alpha_true, TR1, TR2, T1_tissue)
S1_n = S1_true + rng.normal(0, noise, N)
S2_n = S2_true + rng.normal(0, noise, N)

r = S1_n / S2_n
arg = np.clip((r * n - 1) / (n - r), -1, 1)
F_afi = np.arccos(arg) / alpha_nom

pixel = np.arange(N)
rmse = np.sqrt(np.mean((F_afi - F_true)**2))

fig = go.Figure()
fig.add_trace(go.Scatter(x=pixel, y=F_true, name='True B1 factor',
                         line=dict(color='gray', dash='dot', width=1.5)))
fig.add_trace(go.Scatter(x=pixel, y=F_afi,
                         name='AFI estimate (noisy)',
                         line=dict(color='royalblue', width=2)))

fig.update_layout(
    title=f'AFI B1 Map Reconstruction (α_nom=60°, n={n}, TR1={TR1:.0f} ms, RMSE={rmse:.3f})',
    xaxis_title='Pixel Index',
    yaxis_title='B1 Scaling Factor F',
    template='plotly_white',
    height=360,
    legend=dict(x=0.98, y=0.98, xanchor='right', yanchor='top'),
)
fig.show()
```

The AFI B1 map is reconstructed pixel-by-pixel using [Eq. %s](eq-afi-alpha).
Noise in the source images propagates through the ratio and the $\arccos$
operation. Gaussian smoothing of the B1 map is often applied in practice,
since B1 varies smoothly over the scale of the RF wavelength and is not
expected to show sharp spatial features.

## Example Brain Maps

A simulated AFI acquisition on a 3 T brain phantom demonstrates the
two interleaved SPGR images and the resulting B1 map.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

rng = np.random.default_rng(1)

def brain_phantom(n=180):
    y_i, x_i = np.mgrid[:n, :n]
    cx, cy = n*0.5, n*0.5
    dx = (x_i - cx) / (n*0.46); dy = (y_i - cy) / (n*0.42)
    def ell(ax, ay, dx0=0., dy0=0.):
        return ((dx+dx0)/ax)**2 + ((dy+dy0)/ay)**2 < 1.0
    mh=ell(1.00,1.00); mb=ell(0.93,0.91); mw=ell(0.76,0.74)
    mvl=ell(0.17,0.11,-0.30,0); mvr=ell(0.17,0.11,+0.30,0)
    mbl=ell(0.12,0.10,-0.38,-0.20); mbr=ell(0.12,0.10,+0.38,-0.20)
    mtl=ell(0.11,0.09,-0.13,-0.09); mtr_=ell(0.11,0.09,+0.13,-0.09)
    t=np.zeros((n,n),dtype=int)
    t[mh]=1; t[mb]=3; t[mw]=4; t[mbl|mbr]=5; t[mtl|mtr_]=6; t[mvl|mvr]=7
    B1=np.where(mb, 1.0+0.22*np.exp(-4.0*(dx**2+dy**2)), 0.0)
    return t, mb, B1

T1v={0:0,1:0,3:1300,4:840, 5:1200,6:1100,7:4500}
PDv={0:0,1:0,3:0.85,4:0.75,5:0.80,6:0.78,7:1.00}
tissue, bmask, B1m = brain_phantom()
T1m = np.vectorize(T1v.get)(tissue).astype(float)
PDm = np.vectorize(PDv.get)(tissue).astype(float)

TR1, TR2, fa_nom = 20.0, 100.0, 60.0   # ms, ms, degrees
n_afi = TR1 / TR2
fa_true = np.radians(fa_nom * B1m)

E1 = np.exp(-TR1 / np.where(T1m > 0, T1m, 1))
E2 = np.exp(-TR2 / np.where(T1m > 0, T1m, 1))

S1 = PDm * np.sin(fa_true) * (1 - E2) / (1 - E1 * E2 * np.cos(fa_true) + 1e-10)
S2 = PDm * np.sin(fa_true) * (1 - E1) / (1 - E1 * E2 * np.cos(fa_true) + 1e-10)

def rician(img, s):
    re = img + rng.normal(0, s, img.shape)
    im = rng.normal(0, s, img.shape)
    return np.sqrt(re**2 + im**2)

img_s1 = rician(S1 * bmask, 0.02)
img_s2 = rician(S2 * bmask, 0.02)

r = img_s1 / np.where(bmask & (img_s2 > 0.005), img_s2, np.nan)
arg = np.clip((r * n_afi - 1) / (n_afi - r + 1e-10), -1, 1)
B1_afi = np.where(bmask, np.degrees(np.arccos(arg)) / fa_nom, np.nan)

titles = [f'S₁  (TR₁={int(TR1)} ms)', f'S₂  (TR₂={int(TR2)} ms)', 'AFI B1 map']
fig = make_subplots(rows=1, cols=3, subplot_titles=titles, horizontal_spacing=0.06)

for i, img in enumerate([img_s1, img_s2], 1):
    fig.add_trace(go.Heatmap(z=np.flipud(img), colorscale='gray',
                             zmin=0, zmax=0.5, showscale=False), row=1, col=i)

fig.add_trace(go.Heatmap(z=np.flipud(B1_afi), colorscale='RdBu',
                         zmin=0.7, zmax=1.3, showscale=True,
                         colorbar=dict(title='B1 factor', len=0.75, thickness=15)),
              row=1, col=3)

fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False)
fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False)
fig.update_layout(
    title=dict(text=f'AFI Brain Maps  (α_nom = {int(fa_nom)}°, TR₁ = {int(TR1)} ms, TR₂ = {int(TR2)} ms)', x=0.5),
    height=320, template='plotly_white',
)
fig.show()
```

The S₁ image (short TR) has strong T1 weighting, while S₂ (long TR)
is closer to PD-weighted. Both images have the same flip angle and
hence the same B1 spatial pattern embedded in their signal level.
Forming the ratio $r = S_1/S_2$ and inverting Eq. [%s](eq-afi-ratio)
isolates the B1 field, which matches the known centre-bright phantom
input.

## References

```{bibliography}
:filter: docname in docnames
```
