---
title: "3.2 Variable Flip Angle"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 3.2 Variable Flip Angle

The variable flip angle (VFA) method — also called DESPOT1 (driven equilibrium
single pulse observation of T1) — maps T1 from a series of spoiled gradient
echo (SPGR) images acquired at different flip angles with a fixed TR
{cite}`fram1987,deoni2003`. Because each acquisition is a standard 3D SPGR
sequence, VFA is inherently fast: a T1 map of the whole brain can be obtained
in a few minutes compared to 20–40 minutes for multi-TI inversion recovery.

## Pulse Sequence

VFA acquires $N$ spoiled GRE images ($N \geq 2$, typically 2–6) with:

- **Identical TR** (typically 5–30 ms) and TE.
- **Different flip angles** $\alpha_1, \alpha_2, \ldots, \alpha_N$,
  spanning a range that brackets the Ernst angle for the expected T1
  (e.g., 2°–20° for T1 ≈ 1 s with TR = 20 ms).

Between acquisitions, spoiler gradients destroy transverse coherence,
ensuring each acquisition reaches a new steady state independently.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def rect(t, t0, w, h=1.0):
    r = np.zeros_like(t); r[(t >= t0) & (t <= t0 + w)] = h; return r

def sinc_pulse(t, t0, w, h=1.0):
    x = (t - (t0 + w/2)) / (w/6); v = h * np.sinc(x)
    v[(t < t0) | (t > t0 + w)] = 0; return v

fig = plt.figure(figsize=(11, 5))
gs = gridspec.GridSpec(4, 3, hspace=0.06, wspace=0.32,
                       left=0.06, right=0.97, top=0.88, bottom=0.13)
t = np.linspace(0, 1, 1000)
TR = 1.0; TE = 0.28

fa_configs = [(0.12, '#1f77b4', 'α₁  (small)'),
              (0.42, '#ff7f0e', 'α₂  (medium)'),
              (0.72, '#d62728', 'α₃  (large)')]

for col_idx, (fa, color, label) in enumerate(fa_configs):
    axes = [fig.add_subplot(gs[i, col_idx]) for i in range(4)]
    plt.setp([a.get_xticklabels() for a in axes[:-1]], visible=False)

    rf = sinc_pulse(t, 0.04, 0.07, fa)
    axes[0].fill_between(t, 0, rf, color=color, alpha=0.85)
    axes[0].axhline(0, color='k', linewidth=0.7)
    axes[0].annotate(label.split('(')[0].strip(), (0.075, fa + 0.04),
                     ha='center', fontsize=9)
    axes[0].set_ylabel('RF', fontsize=8); axes[0].set_ylim(-0.1, 0.90)
    axes[0].yaxis.set_tick_params(labelleft=False)
    axes[0].set_title(label, fontsize=9, pad=4)

    gz = rect(t, 0.03, 0.09, 0.70) + rect(t, 0.12, 0.035, -0.35)
    axes[1].fill_between(t, 0, gz, where=gz > 0, color='forestgreen', alpha=0.7)
    axes[1].fill_between(t, 0, gz, where=gz < 0, color='tomato', alpha=0.7)
    axes[1].axhline(0, color='k', linewidth=0.7)
    axes[1].set_ylabel('Gz', fontsize=8); axes[1].yaxis.set_tick_params(labelleft=False)

    gy = rect(t, TE - 0.03, 0.05, 0.55)
    axes[2].fill_between(t, 0, gy, color='darkorange', alpha=0.7)
    axes[2].axhline(0, color='k', linewidth=0.7)
    axes[2].set_ylabel('Gy', fontsize=8); axes[2].yaxis.set_tick_params(labelleft=False)

    gx = rect(t, TE - 0.03, 0.05, -0.42) + rect(t, TE + 0.03, 0.12, 0.42) + \
         rect(t, TR - 0.09, 0.07, 0.55)
    axes[3].fill_between(t, 0, gx, where=gx > 0, color='purple', alpha=0.7)
    axes[3].fill_between(t, 0, gx, where=gx < 0, color='plum', alpha=0.7)
    axes[3].fill_between(t, -0.15, -0.05,
                         where=((t >= TE+0.03) & (t <= TE+0.15)), color='gold', alpha=0.9)
    axes[3].axhline(0, color='k', linewidth=0.7)
    axes[3].set_ylabel('Gx / ADC', fontsize=8); axes[3].yaxis.set_tick_params(labelleft=False)
    axes[3].set_xlabel('Time (one TR)', fontsize=8)
    axes[3].set_xlim(0, TR); axes[3].set_ylim(-0.25, 0.75)

fig.suptitle('Variable Flip Angle T1 Mapping — Three Spoiled-GRE Acquisitions '
             '(same TR, different α)', fontsize=11, y=0.97)
plt.close()
```

:::{figure} vfa_pulse_sequence.png
:name: fig-vfa-sequence
:align: center
VFA pulse sequence. Three spoiled GRE acquisitions are shown with
the same TR but increasing flip angle ($\alpha_1 < \alpha_2 < \alpha_3$).
The steady-state signal depends on both T1 and $\alpha$ through the
Ernst equation; sampling multiple flip angles provides enough
information to separate T1 from $M_0$.
:::

## Mathematical Model

### Ernst Equation

The steady-state spoiled GRE signal is {cite}`ernst1966`:

$$
S(\alpha) = M_0\sin\alpha\,\frac{1 - E_1}{1 - E_1\cos\alpha},
\qquad E_1 = e^{-\text{TR}/T_1}
$$ (eq-vfa-ernst)

At the **Ernst angle** $\alpha_E = \arccos(E_1)$ the signal is
maximised for a given T1 and TR.

### Linearisation

Dividing both sides of [Eq. %s](eq-vfa-ernst) by $\sin\alpha$ and
recognising $S/\tan\alpha = S\cos\alpha/\sin\alpha$, one obtains the
Deoni linearisation {cite}`deoni2003`:

$$
\frac{S}{\sin\alpha} = E_1\,\frac{S}{\tan\alpha}
  + M_0(1 - E_1)
$$ (eq-vfa-linear)

Plotting $y = S/\sin\alpha$ versus $x = S/\tan\alpha$ yields a straight
line with slope $E_1 = e^{-\text{TR}/T_1}$ and intercept $M_0(1-E_1)$.
T1 is recovered from the slope:

$$
T_1 = -\frac{\text{TR}}{\ln(E_1)}
$$ (eq-vfa-T1)

This linearisation is computationally trivial (a least-squares fit to a
line) but requires adequate spread of flip angles to sample the Ernst
curve on both sides of the Ernst angle.

### B1 Sensitivity

The actual flip angle is $\alpha_\text{actual} = F \cdot \alpha_\text{nominal}$,
where $F(\mathbf{r})$ is the B1 scaling factor (Chapter 2). Substituting
into [Eq. %s](eq-vfa-ernst):

$$
S(\alpha) = M_0\sin(F\alpha)\,\frac{1 - E_1}{1 - E_1\cos(F\alpha)}
$$

Ignoring B1 inhomogeneity (assuming $F = 1$ everywhere) introduces a
systematic T1 bias. For a 20% B1 deviation ($F = 0.8$), T1 errors can
exceed 50% {cite}`stollberger1996`. **B1 correction is therefore
mandatory for quantitative VFA T1 mapping.**

B1 correction is applied by replacing $\alpha_\text{nominal}$ with
$F \cdot \alpha_\text{nominal}$ at every voxel before performing the
linear fit.

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def ernst(alpha_deg, T1, TR, M0=1.0):
    alpha = np.radians(alpha_deg)
    E1 = np.exp(-TR / T1)
    return M0 * np.sin(alpha) * (1 - E1) / (1 - E1 * np.cos(alpha))

alpha_arr = np.linspace(1, 90, 500)
TR = 20   # ms

fig = make_subplots(rows=1, cols=2,
                    subplot_titles=['Ernst Curves for Different T1 Values',
                                    'T1 Bias Due to B1 Inhomogeneity (TR=20 ms)'])

tissues = [('WM  T1=840 ms',  840,  '#1f77b4'),
           ('GM  T1=1300 ms', 1300, '#d62728'),
           ('CSF T1=4500 ms', 4500, '#2ca02c')]

for name, T1, col in tissues:
    S = ernst(alpha_arr, T1, TR)
    fig.add_trace(go.Scatter(x=alpha_arr, y=S, name=name,
                             line=dict(color=col, width=2.2)), row=1, col=1)
    # Ernst angle
    E1 = np.exp(-TR / T1)
    ae = np.degrees(np.arccos(E1))
    Se = ernst(ae, T1, TR)
    fig.add_trace(go.Scatter(x=[ae], y=[Se], mode='markers',
                             marker=dict(symbol='star', color=col, size=10),
                             showlegend=False), row=1, col=1)

# Right: T1 bias vs true B1 for typical FA choices
FA1, FA2 = 5, 15   # degrees, typical VFA pair for T1~1s with TR=20ms
T1_test = np.array([840, 1300])
T1_names = ['WM (T1=840 ms)', 'GM (T1=1300 ms)']
T1_colors = ['#1f77b4', '#d62728']
B1_arr = np.linspace(0.5, 1.5, 300)

for T1, name, col in zip(T1_test, T1_names, T1_colors):
    E1_true = np.exp(-TR / T1)
    T1_estimates = []
    for b1 in B1_arr:
        # Signal acquired with nominal FA but actual FA = b1 * FA
        S1 = ernst(b1 * FA1, T1, TR)
        S2 = ernst(b1 * FA2, T1, TR)
        # Fit assuming B1=1 (no correction)
        y1 = S1 / np.sin(np.radians(FA1)); x1 = S1 / np.tan(np.radians(FA1))
        y2 = S2 / np.sin(np.radians(FA2)); x2 = S2 / np.tan(np.radians(FA2))
        if abs(x2 - x1) > 1e-10:
            E1_est = (y2 - y1) / (x2 - x1)
            E1_est = np.clip(E1_est, 1e-6, 1 - 1e-6)
            T1_est = -TR / np.log(E1_est)
        else:
            T1_est = np.nan
        T1_estimates.append((T1_est - T1) / T1 * 100)
    fig.add_trace(go.Scatter(x=B1_arr, y=T1_estimates, name=name,
                             line=dict(color=col, width=2.2)), row=1, col=2)

fig.add_hline(y=0, line_dash='dot', line_color='gray', row=1, col=2)
fig.add_vline(x=1.0, line_dash='dot', line_color='lightgray', row=1, col=2)

fig.update_xaxes(title_text='Flip Angle α (°)', row=1, col=1)
fig.update_xaxes(title_text='True B1 Factor F', row=1, col=2)
fig.update_yaxes(title_text='Signal (a.u.)', row=1, col=1)
fig.update_yaxes(title_text='T1 Bias (%)', row=1, col=2)
fig.update_layout(
    height=430,
    title_text=f'VFA Simulation  (TR = {TR} ms;  ★ = Ernst angle;  α₁={FA1}°, α₂={FA2}°)',
    legend=dict(orientation='h', yanchor='bottom', y=-0.28,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows the Ernst curves for WM, GM, and CSF at TR = 20 ms.
Stars mark the Ernst angle where signal is maximised for each tissue.
Optimal VFA sampling brackets the Ernst angle. The right panel
quantifies the T1 bias when B1 correction is omitted: a B1 factor of
0.8 (20% under-estimation of flip angle) causes ~80% T1 over-estimation
in WM and ~60% in GM, demonstrating that VFA without B1 correction is
unreliable in inhomogeneous B1 environments.

## Data Fitting

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go

rng = np.random.default_rng(3)

N = 200; TR = 20.0
FA_nom_deg = np.array([2, 5, 10, 15, 20])
FA_nom     = np.radians(FA_nom_deg)

x = np.linspace(0, np.pi, N)
T1_true = 700 + 700 * np.sin(x / 2)**2      # smooth profile 700–1400 ms
B1_true = 0.75 + 0.50 * np.cos(x - np.pi/3)**2   # B1 varies 0.75–1.25
M0_true = 0.9 + 0.2 * np.sin(x)

noise = 0.010

def vfa_fit_linear(S_arr, FA_arr, TR):
    y = S_arr / np.sin(FA_arr)
    x = S_arr / np.tan(FA_arr)
    coeffs = np.polyfit(x, y, 1)   # slope = E1, intercept = M0(1-E1)
    E1 = np.clip(coeffs[0], 1e-6, 1 - 1e-6)
    T1 = -TR / np.log(E1)
    M0 = coeffs[1] / (1 - E1)
    return T1, M0

T1_no_b1  = np.zeros(N)
T1_with_b1 = np.zeros(N)

for i in range(N):
    FA_actual = B1_true[i] * FA_nom
    E1 = np.exp(-TR / T1_true[i])
    S_clean = M0_true[i] * np.sin(FA_actual) * (1 - E1) / (1 - E1 * np.cos(FA_actual))
    S_noisy = S_clean + rng.normal(0, noise, len(FA_nom))
    # Without B1 correction (assume FA_actual = FA_nominal)
    T1_no_b1[i],  _ = vfa_fit_linear(S_noisy, FA_nom, TR)
    # With B1 correction (use known B1 to correct flip angles)
    T1_with_b1[i], _ = vfa_fit_linear(S_noisy, FA_actual, TR)

pixel = np.arange(N)
rmse_no   = np.sqrt(np.mean((T1_no_b1  - T1_true)**2))
rmse_corr = np.sqrt(np.mean((T1_with_b1 - T1_true)**2))

fig = go.Figure()
fig.add_trace(go.Scatter(x=pixel, y=T1_true,
                         name='True T1 (ms)',
                         line=dict(color='gray', dash='dot', width=1.5)))
fig.add_trace(go.Scatter(x=pixel, y=T1_no_b1,
                         name=f'VFA no B1 corr. (RMSE={rmse_no:.0f} ms)',
                         line=dict(color='tomato', width=2)))
fig.add_trace(go.Scatter(x=pixel, y=T1_with_b1,
                         name=f'VFA with B1 corr. (RMSE={rmse_corr:.0f} ms)',
                         line=dict(color='royalblue', width=2)))

fig.update_layout(
    title=f'VFA T1 Reconstruction  ({len(FA_nom_deg)} FAs: {list(FA_nom_deg)}°,  TR={TR:.0f} ms)',
    xaxis_title='Pixel index',
    yaxis_title='T1 (ms)',
    template='plotly_white',
    height=380,
    legend=dict(x=0.01, y=0.99, xanchor='left', yanchor='top'),
)
fig.show()
```

The red curve (no B1 correction) shows substantial systematic error
wherever B1 deviates from 1, distorting the spatial T1 profile.
After B1 correction (blue), the reconstruction closely follows the
true T1 with RMSE driven primarily by noise. This comparison
underscores that B1 maps are not optional in VFA — they are a
fundamental prerequisite.

## Example Brain Maps

A three-flip-angle VFA acquisition is simulated with the centre-bright
3 T B1 field. Two T1 maps are reconstructed — one ignoring B1, one
correcting for it — to show the spatial bias introduced by B1
inhomogeneity.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

rng = np.random.default_rng(4)

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
n   = tissue.shape[0]

FA_deg = np.array([3., 10., 20.])
FA_rad = np.radians(FA_deg)
TR = 20.0

E1m = np.exp(-TR / np.where(T1m > 0, T1m, 1))

def spgr(fa_rad, E1, PD):
    return PD * np.sin(fa_rad) * (1 - E1) / (1 - E1 * np.cos(fa_rad) + 1e-10)

def rician(img, s):
    re = img + rng.normal(0, s, img.shape)
    im = rng.normal(0, s, img.shape)
    return np.sqrt(re**2 + im**2)

# Acquire: actual FA = B1 * nominal FA
imgs = [rician(spgr(f * B1m, E1m, PDm) * bmask, 0.020) for f in FA_rad]

def vfa_t1(imgs, FA_rad_used, TR):
    """Vectorised Deoni linearisation."""
    y = np.stack([imgs[i] / np.sin(FA_rad_used[i]) for i in range(len(imgs))])
    x = np.stack([imgs[i] / np.tan(FA_rad_used[i]) for i in range(len(imgs))])
    xm = x.mean(0); ym = y.mean(0)
    xc = x - xm;    yc = y - ym
    E1 = np.clip((xc * yc).sum(0) / np.maximum((xc**2).sum(0), 1e-9), 1e-6, 1 - 1e-6)
    return np.where(bmask, -TR / np.log(E1), np.nan)

# Without B1 correction: assume FA = nominal
T1_no_b1   = np.clip(vfa_t1(imgs, FA_rad, TR), 100, 6000)
# With B1 correction: use measured B1m to correct flip angles
T1_with_b1 = np.clip(vfa_t1(imgs, FA_rad * B1m[np.newaxis, ...], TR), 100, 6000)

titles = [f'α = {int(fa)}°' for fa in FA_deg] + ['VFA T1 (no B1 corr.)', 'VFA T1 (B1 corrected)']
fig = make_subplots(rows=1, cols=5, subplot_titles=titles, horizontal_spacing=0.04)
for i, img in enumerate(imgs, 1):
    fig.add_trace(go.Heatmap(z=np.flipud(img), colorscale='gray',
                             zmin=0, zmax=0.50, showscale=False), row=1, col=i)
for i, data in enumerate([T1_no_b1, T1_with_b1], 4):
    fig.add_trace(go.Heatmap(z=np.flipud(data), colorscale='Magma',
                             zmin=400, zmax=4500, showscale=(i == 5),
                             colorbar=dict(title='T1 (ms)', len=0.75, thickness=15)),
                 row=1, col=i)
fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False)
fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False)
fig.update_layout(
    title=dict(text=f'VFA Brain Maps  (TR = {int(TR)} ms, α = {list(FA_deg.astype(int))}°)', x=0.5),
    height=280, template='plotly_white',
)
fig.show()
```

The uncorrected T1 map (fourth panel) shows a concentric ring artefact:
T1 is under-estimated at the brain centre (where B1 > 1 inflates the
effective flip angle) and over-estimated at the periphery (where B1 < 1).
After B1 correction (fifth panel) the T1 values are uniform within each
tissue, consistent with the known phantom ground truth.

## References

```{bibliography}
:filter: docname in docnames
```
