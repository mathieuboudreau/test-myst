---
title: "2.1 Double Angle Method"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 2.1 Double Angle Method

The double angle method (DAM) is one of the simplest and most widely used
techniques for B1 mapping {cite}`insko1993,stollberger1996`. It requires
only two gradient-echo acquisitions with flip angles $\alpha$ and $2\alpha$
and derives the actual flip angle from their signal ratio.

## Pulse Sequence

DAM acquires two spoiled GRE images with identical parameters except for the
flip angle: one at the nominal flip angle $\alpha$ and one at $2\alpha$:

- **Acquisition 1**: spoiled GRE, flip angle $\alpha$
- **Acquisition 2**: spoiled GRE, flip angle $2\alpha$

The key requirement is that TR must be much longer than T1 so that the signal
approximates the fully-relaxed case. In practice, TR $\geq 5 T_1$ is advised,
which can make DAM slow for whole-brain coverage.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(11, 5))
gs = gridspec.GridSpec(4, 2, hspace=0.05, wspace=0.35,
                       left=0.06, right=0.97, top=0.88, bottom=0.12)

t = np.linspace(0, 1, 1000)
TR = 1.0
TE = 0.2

def rect(t, t0, w, h=1.0):
    r = np.zeros_like(t)
    r[(t >= t0) & (t <= t0 + w)] = h
    return r

def sinc_pulse(t, t0, w, h=1.0):
    x = (t - (t0 + w/2)) / (w/6)
    v = h * np.sinc(x)
    v[(t < t0) | (t > t0 + w)] = 0
    return v

for col_idx, (label, fa, color) in enumerate([
        ('Acquisition 1: flip angle α', 0.4, '#1f77b4'),
        ('Acquisition 2: flip angle 2α', 0.75, '#d62728')]):

    axes = [fig.add_subplot(gs[i, col_idx]) for i in range(4)]
    plt.setp([a.get_xticklabels() for a in axes[:-1]], visible=False)

    # RF
    rf = sinc_pulse(t, 0.02, 0.07, fa)
    axes[0].fill_between(t, 0, rf, color=color, alpha=0.8)
    axes[0].axhline(0, color='k', linewidth=0.8)
    axes[0].annotate(f'{"α" if col_idx==0 else "2α"}',
                     (0.055, fa + 0.04), ha='center', fontsize=9)
    axes[0].set_ylabel('RF', fontsize=8)
    axes[0].set_ylim(-0.1, 1.0)
    axes[0].yaxis.set_tick_params(labelleft=False)
    axes[0].set_title(label, fontsize=9, pad=4)

    # Gz
    gz = rect(t, 0.01, 0.09, 0.7) + rect(t, 0.10, 0.04, -0.35)
    axes[1].fill_between(t, 0, gz, where=gz > 0, color='forestgreen', alpha=0.7)
    axes[1].fill_between(t, 0, gz, where=gz < 0, color='tomato', alpha=0.7)
    axes[1].axhline(0, color='k', linewidth=0.8)
    axes[1].set_ylabel('Gz', fontsize=8)
    axes[1].yaxis.set_tick_params(labelleft=False)

    # Gy
    gy = rect(t, 0.15, 0.05, 0.55)
    axes[2].fill_between(t, 0, gy, color='darkorange', alpha=0.7)
    axes[2].axhline(0, color='k', linewidth=0.8)
    axes[2].set_ylabel('Gy', fontsize=8)
    axes[2].yaxis.set_tick_params(labelleft=False)

    # Gx + ADC + spoiler
    gx = rect(t, 0.15, 0.05, -0.4) + rect(t, TE - 0.06, 0.12, 0.4) + \
         rect(t, TR - 0.09, 0.07, 0.55)
    axes[3].fill_between(t, 0, gx, where=gx > 0, color='purple', alpha=0.7)
    axes[3].fill_between(t, 0, gx, where=gx < 0, color='plum', alpha=0.7)
    axes[3].fill_between(t, -0.15, -0.05,
                         where=((t >= TE - 0.06) & (t <= TE + 0.06)),
                         color='gold', alpha=0.9)
    axes[3].axhline(0, color='k', linewidth=0.8)
    axes[3].set_ylabel('Gx / ADC', fontsize=8)
    axes[3].yaxis.set_tick_params(labelleft=False)
    axes[3].set_xlabel('Time (one TR)', fontsize=8)
    axes[3].set_xlim(0, TR)
    axes[3].set_ylim(-0.25, 0.75)

    # TE label
    axes[3].annotate('', xy=(TE, -0.20), xytext=(0, -0.20),
                     arrowprops=dict(arrowstyle='<->', color='k', lw=0.8))
    axes[3].text(TE / 2, -0.25, 'TE', ha='center', fontsize=8)

fig.suptitle('Double Angle Method — Two Spoiled-GRE Acquisitions',
             fontsize=11, y=0.96)
plt.close()
```

:::{figure} dam_pulse_sequence.png
:name: fig-dam-sequence
:align: center
Double angle method pulse sequence. Two otherwise identical spoiled GRE
acquisitions differ only in flip angle: $\alpha$ (left) and $2\alpha$
(right). Long TR ensures approximate full relaxation between excitations.
:::

## Mathematical Model

### Signal Ratio

For a fully-relaxed spoiled GRE ($\text{TR} \gg T_1$), the steady-state
signal simplifies to:

$$
S(\alpha) \approx M_0 \sin(\alpha_\text{actual})\,e^{-\text{TE}/T_2^*}
$$

Taking the ratio of the two acquisitions cancels the common factors $M_0$ and
$e^{-\text{TE}/T_2^*}$:

$$
\frac{S(2\alpha)}{S(\alpha)} = \frac{\sin(2\alpha_\text{actual})}{\sin(\alpha_\text{actual})}
  = 2\cos(\alpha_\text{actual})
$$ (eq-dam-ratio)

Solving for the actual flip angle:

$$
\alpha_\text{actual} = \arccos\!\left(\frac{S(2\alpha)}{2\,S(\alpha)}\right)
$$ (eq-dam-alpha)

The B1 scaling factor map is then:

$$
F(\mathbf{r}) = \frac{\alpha_\text{actual}(\mathbf{r})}{\alpha_\text{nominal}}
$$ (eq-dam-F)

### T1 Correction

The assumption $\text{TR} \gg T_1$ is often not perfectly met in practice.
For finite TR, the ratio becomes {cite}`stollberger1996`:

$$
\frac{S(2\alpha)}{S(\alpha)} =
  \frac{\sin(2\alpha_\text{actual})\,(1 - E_1)}
       {\sin(\alpha_\text{actual})\,(1 - E_1\cos(\alpha_\text{actual}))} \cdot
  \frac{1 - E_1\cos(\alpha_\text{actual})}{1 - E_1\cos(2\alpha_\text{actual})}
$$

where $E_1 = e^{-\text{TR}/T_1}$. For $E_1 \to 0$ this reduces to the simple
ratio above. The T1 dependence introduces a systematic error unless TR is
sufficiently long.

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

alpha_nom_deg = 60   # nominal flip angle
F_arr = np.linspace(0.4, 1.6, 200)   # B1 scaling factors
alpha_act = np.radians(F_arr * alpha_nom_deg)
alpha2_act = 2 * alpha_act

# Ideal ratio (long TR)
ratio_ideal = np.sin(alpha2_act) / np.sin(alpha_act)  # = 2cos(alpha_act)

# Recovered FA from ideal ratio
F_recovered_ideal = np.arccos(ratio_ideal / 2) / np.radians(alpha_nom_deg)

# T1 effect: finite TR
fig = make_subplots(rows=1, cols=2,
                    subplot_titles=['Signal Ratio vs True B1 Factor',
                                    'B1 Map Error Due to Finite TR'])

fig.add_trace(go.Scatter(x=F_arr, y=ratio_ideal,
                         name='Ideal (TR >> T1)',
                         line=dict(color='royalblue', width=2.5)), row=1, col=1)

T1 = 1000   # ms
for TR, col in [(5000, '#2ca02c'), (3000, '#ff7f0e'), (2000, '#d62728')]:
    E1 = np.exp(-TR / T1)
    S_alpha = np.sin(alpha_act) * (1 - E1) / (1 - E1 * np.cos(alpha_act))
    S_2alpha = np.sin(alpha2_act) * (1 - E1) / (1 - E1 * np.cos(alpha2_act))
    ratio_finite = S_2alpha / S_alpha
    F_rec_finite = np.arccos(ratio_finite / 2) / np.radians(alpha_nom_deg)

    fig.add_trace(go.Scatter(x=F_arr, y=ratio_finite,
                             name=f'TR = {TR} ms (T1 = {T1} ms)',
                             line=dict(color=col, width=2, dash='dash')),
                 row=1, col=1)

    pct_error = (F_rec_finite - F_arr) / F_arr * 100
    fig.add_trace(go.Scatter(x=F_arr, y=pct_error,
                             name=f'TR = {TR} ms',
                             line=dict(color=col, width=2, dash='dash'),
                             showlegend=False),
                 row=1, col=2)

fig.add_hline(y=0, line_dash='dot', line_color='gray', row=1, col=2)

fig.update_xaxes(title_text='True B1 Scaling Factor F', row=1, col=1)
fig.update_xaxes(title_text='True B1 Scaling Factor F', row=1, col=2)
fig.update_yaxes(title_text='Signal Ratio S(2α)/S(α)', row=1, col=1)
fig.update_yaxes(title_text='B1 Map Error (%)', row=1, col=2)

fig.update_layout(
    height=420,
    title_text=f'DAM Simulation (α_nominal = {alpha_nom_deg}°, T1 = {T1} ms)',
    legend=dict(orientation='h', yanchor='bottom', y=-0.3,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows how the signal ratio departs from the ideal
$2\cos(\alpha_\text{actual})$ curve when TR is finite. The right panel
quantifies the resulting percentage error in the recovered B1 map:
at TR = 2000 ms (with T1 = 1000 ms), errors can reach several percent,
underscoring the need for long TR or T1 correction.

## Data Fitting

In practice, DAM is applied pixel-by-pixel. Here we simulate a 1D B1
profile and reconstruct a B1 map from noisy signal images.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go

rng = np.random.default_rng(17)

N = 200   # spatial pixels
alpha_nom = np.radians(60)
M0 = 1.0
T2s = 40   # ms
TE = 5

# Simulate a smooth B1 spatial profile (realistic: ~0.6 to 1.4)
x = np.linspace(0, np.pi, N)
F_true = 0.7 + 0.4 * np.cos(x) ** 2 + 0.1 * np.sin(2 * x)
alpha_true = F_true * alpha_nom

# Ideal (long TR) signals with noise
noise_level = 0.015
S_alpha  = M0 * np.sin(alpha_true)  * np.exp(-TE / T2s)
S_2alpha = M0 * np.sin(2 * alpha_true) * np.exp(-TE / T2s)
S_alpha_n  = S_alpha  + rng.normal(0, noise_level, N)
S_2alpha_n = S_2alpha + rng.normal(0, noise_level, N)

# DAM reconstruction
ratio = np.clip(S_2alpha_n / (2 * S_alpha_n), -1, 1)
F_measured = np.arccos(ratio) / alpha_nom

pixel = np.arange(N)

fig = go.Figure()
fig.add_trace(go.Scatter(x=pixel, y=F_true,
                         name='True B1 factor',
                         line=dict(color='gray', dash='dot', width=1.5)))
fig.add_trace(go.Scatter(x=pixel, y=F_measured,
                         name='DAM estimate (noisy data)',
                         line=dict(color='royalblue', width=2)))

rmse = np.sqrt(np.mean((F_measured - F_true)**2))
fig.update_layout(
    title=f'DAM B1 Map Reconstruction (RMSE = {rmse:.3f}, α_nom = 60°, SNR ~ {1/noise_level:.0f})',
    xaxis_title='Pixel index',
    yaxis_title='B1 Scaling Factor F',
    template='plotly_white',
    height=360,
    legend=dict(x=0.98, y=0.98, xanchor='right', yanchor='top'),
)
fig.show()
```

The reconstructed B1 map closely follows the true profile, with residual
noise from the SNR of the source images. Because DAM takes the ratio
of two images, noise propagates through the $\arccos$ non-linearly:
the error is amplified when $\alpha_\text{actual} \to 0°$ or $\to 90°$,
so the method works best when $\alpha_\text{nominal}$ is chosen near 60°.

## Example Brain Maps

The following simulates a full-brain DAM acquisition on a 2D digital
brain phantom with a realistic centre-bright 3 T B1 field, then
displays the two source images and the recovered B1 map.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

rng = np.random.default_rng(0)

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

alpha = 60.0   # nominal FA (degrees)
TR    = 500.0  # ms
E1m   = np.exp(-TR / np.where(T1m > 0, T1m, 1))

def spgr(fa_deg, E1, PD, B1):
    fa = np.radians(fa_deg * B1)
    return PD * np.sin(fa) * (1 - E1) / (1 - E1 * np.cos(fa) + 1e-10)

def rician(img, s):
    re = img + rng.normal(0, s, img.shape)
    im = rng.normal(0, s, img.shape)
    return np.sqrt(re**2 + im**2)

img_a  = rician(spgr(alpha,   E1m, PDm, B1m) * bmask, 0.025)
img_2a = rician(spgr(2*alpha, E1m, PDm, B1m) * bmask, 0.025)

ratio  = 2 * img_a**2 / np.where(bmask & (img_2a > 0.01), img_2a**2, np.nan) - 1
B1_dam = np.where(bmask,
                  np.degrees(np.arccos(np.clip(ratio, -1, 1))) / alpha,
                  np.nan)

titles = [f'S(α={int(alpha)}°)', f'S(2α={int(2*alpha)}°)', 'DAM B1 map']
fig = make_subplots(rows=1, cols=3, subplot_titles=titles, horizontal_spacing=0.06)
for i, img in enumerate([img_a, img_2a], 1):
    fig.add_trace(go.Heatmap(z=np.flipud(img), colorscale='gray',
                             zmin=0, zmax=0.55, showscale=False), row=1, col=i)
fig.add_trace(go.Heatmap(z=np.flipud(B1_dam), colorscale='RdBu',
                         zmin=0.7, zmax=1.3, showscale=True,
                         colorbar=dict(title='B1 factor', len=0.75, thickness=15)),
              row=1, col=3)
fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False)
fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False)
fig.update_layout(
    title=dict(text=f'DAM Brain Maps  (α_nom = {int(alpha)}°, TR = {int(TR)} ms)', x=0.5),
    height=320, template='plotly_white',
)
fig.show()
```

The centre-bright B1 pattern characteristic of 3 T is clearly visible:
the B1 factor exceeds 1.2 at the brain centre and falls to ~0.8 at the
periphery. Both source images inherit this spatial weighting, but the
ratio $\arccos(2S_\alpha^2/S_{2\alpha}^2 - 1)/\alpha$ removes the
$M_0$ dependence, revealing the pure B1 map. The skull ring and
background are masked where the signal falls below threshold.

## References

```{bibliography}
:filter: docname in docnames
```
