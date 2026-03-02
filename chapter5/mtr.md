---
title: "5.1 MT Ratio"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 5.1 MT Ratio

The magnetization transfer ratio (MTR) is the simplest and most widely
used MT metric in clinical MRI {cite}`wolff1989`. It requires only two
acquisitions — one with and one without an off-resonance MT saturation
pulse — and derives a semi-quantitative index of MT effect from their
signal ratio.

## Pulse Sequence

MTR acquires two spoiled GRE images with identical parameters except
for the presence of the MT preparation pulse:

- **Reference scan** ($M_0$): standard spoiled GRE with no MT pulse.
- **MT-saturated scan** ($M_\text{sat}$): spoiled GRE preceded by an
  off-resonance Gaussian (or sinc-Gaussian) RF pulse applied at frequency
  offset $\pm\Delta\omega$ from the water resonance.

The MT pulse is applied far enough off-resonance (~1–5 kHz) to avoid
direct saturation of free-pool water while effectively saturating the
broad macromolecular resonance.

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

def gauss_pulse(t, t0, w, h=1.0):
    tc = t0 + w/2; sigma = w/5
    v = h * np.exp(-0.5*((t - tc)/sigma)**2)
    v[(t < t0) | (t > t0 + w)] = 0; return v

fig = plt.figure(figsize=(11, 5))
gs = gridspec.GridSpec(4, 2, hspace=0.06, wspace=0.35,
                       left=0.07, right=0.97, top=0.88, bottom=0.13)
t = np.linspace(0, 1, 1000)
TR = 1.0; TE = 0.35

for col_idx, (label, has_mt) in enumerate([
        ('Reference scan  (no MT pulse)', False),
        ('MT-saturated scan  (with MT pulse)', True)]):

    axes = [fig.add_subplot(gs[i, col_idx]) for i in range(4)]
    plt.setp([a.get_xticklabels() for a in axes[:-1]], visible=False)

    if has_mt:
        mt = gauss_pulse(t, 0.02, 0.18, 0.55)
        axes[0].fill_between(t, 0, mt, color='coral', alpha=0.75, label='MT pulse (off-res.)')
        axes[0].annotate('+Δω', (0.11, 0.58), ha='center', fontsize=8, color='coral')
    rf = sinc_pulse(t, 0.25, 0.07, 0.60)
    axes[0].fill_between(t, 0, rf, color='royalblue', alpha=0.85, label='Excitation')
    axes[0].axhline(0, color='k', linewidth=0.7)
    axes[0].annotate('α', (0.285, 0.63), ha='center', fontsize=9)
    axes[0].set_ylabel('RF', fontsize=8); axes[0].set_ylim(-0.1, 0.85)
    axes[0].yaxis.set_tick_params(labelleft=False)
    axes[0].set_title(label, fontsize=9, pad=4)
    if has_mt:
        axes[0].legend(loc='upper right', fontsize=7, framealpha=0.8)

    gz = rect(t, 0.23, 0.10, 0.70) + rect(t, 0.33, 0.04, -0.35)
    axes[1].fill_between(t, 0, gz, where=gz > 0, color='forestgreen', alpha=0.7)
    axes[1].fill_between(t, 0, gz, where=gz < 0, color='tomato', alpha=0.7)
    axes[1].axhline(0, color='k', linewidth=0.7)
    axes[1].set_ylabel('Gz', fontsize=8); axes[1].yaxis.set_tick_params(labelleft=False)

    gy = rect(t, 0.38, 0.05, 0.55)
    axes[2].fill_between(t, 0, gy, color='darkorange', alpha=0.7)
    axes[2].axhline(0, color='k', linewidth=0.7)
    axes[2].set_ylabel('Gy', fontsize=8); axes[2].yaxis.set_tick_params(labelleft=False)

    gx = rect(t, 0.38, 0.05, -0.42) + rect(t, TE - 0.07, 0.14, 0.42) + \
         rect(t, TR - 0.09, 0.07, 0.55)
    axes[3].fill_between(t, 0, gx, where=gx > 0, color='purple', alpha=0.7)
    axes[3].fill_between(t, 0, gx, where=gx < 0, color='plum', alpha=0.7)
    axes[3].fill_between(t, -0.15, -0.05,
                         where=((t >= TE-0.07) & (t <= TE+0.07)), color='gold', alpha=0.9)
    axes[3].axhline(0, color='k', linewidth=0.7)
    axes[3].set_ylabel('Gx / ADC', fontsize=8)
    axes[3].yaxis.set_tick_params(labelleft=False)
    axes[3].set_xlabel('Time (one TR)', fontsize=8)
    axes[3].set_xlim(0, TR); axes[3].set_ylim(-0.25, 0.75)
    sig_lbl = 'M₀' if not has_mt else 'M_sat'
    axes[3].text(TE, -0.22, sig_lbl, ha='center', fontsize=9, color='goldenrod')

fig.suptitle('MT Ratio Pulse Sequence — Two Spoiled-GRE Acquisitions', fontsize=11, y=0.97)
plt.close()
```

:::{figure} mtr_pulse_sequence.png
:name: fig-mtr-sequence
:align: center
MTR pulse sequence. The reference scan (left) acquires a standard
spoiled GRE image. The MT-saturated scan (right) adds a Gaussian
off-resonance preparation pulse at offset $+\Delta\omega$ before the
excitation. The ratio of the two signals yields the MTR map.
:::

## Mathematical Model

### Definition

The MTR is defined as the fractional signal reduction due to the
MT pulse {cite}`wolff1989`:

$$
\text{MTR} = \frac{M_0 - M_\text{sat}}{M_0} \times 100\ \%
$$ (eq-mtr-def)

where $M_0$ is the signal without MT saturation and $M_\text{sat}$ is
the signal with the MT pulse. MTR ranges from 0 % (no MT effect) to
theoretically 100 % (complete saturation transferred to free pool).

### Signal Model

For a spoiled GRE with TR and flip angle $\alpha$, the steady-state
signal without MT is the Ernst equation:

$$
S_0 = M_0 \sin\alpha\,\frac{1 - E_1}{1 - E_1\cos\alpha}
$$

where $E_1 = e^{-\text{TR}/T_1}$. With an MT preparation pulse that
delivers fractional saturation $\delta$ to the longitudinal
magnetisation at the start of each TR, the signal becomes:

$$
S_\text{sat} = M_0 \sin\alpha\,
  \frac{(1-\delta)(1 - E_1)}{1 - (1-\delta)E_1\cos\alpha}
$$ (eq-mtr-sat)

Substituting into [Eq. %s](eq-mtr-def):

$$
\text{MTR} = 1 - \frac{S_\text{sat}}{S_0}
  = 1 - \frac{(1-\delta)(1 - E_1\cos\alpha)}
             {1 - (1-\delta)E_1\cos\alpha}
$$ (eq-mtr-model)

Several key properties follow from this expression:

- MTR increases with increasing $\delta$ (more MT saturation).
- MTR depends on T1 through $E_1 = e^{-\text{TR}/T_1}$: tissues with
  longer T1 show higher MTR at the same $\delta$, because recovered
  magnetisation cannot "refill" the saturated pool as quickly.
- MTR depends on the acquisition parameters TR and $\alpha$, making
  it **not a pure tissue property** — values are not directly comparable
  across different protocols or scanners.

### Limitations of MTR

Because MTR conflates MT physics (δ) with T1 relaxation and protocol
parameters, it is described as *semi-quantitative*. In particular:

- TR and flip angle affect the observed MTR independently of tissue MT.
- B1 inhomogeneity alters the MT pulse efficiency and the excitation
  flip angle simultaneously.
- B0 shifts can move the MT pulse off its intended offset $\Delta\omega$.

Quantitative methods (qMT, MTsat) address these limitations by
separating the underlying tissue parameters from protocol effects.

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def mtr_model(delta, TR, T1, alpha_deg):
    """MTR from the Ernst-based two-parameter model."""
    alpha = np.radians(alpha_deg)
    E1 = np.exp(-TR / T1)
    S0  = np.sin(alpha) * (1 - E1) / (1 - E1 * np.cos(alpha))
    Smt = np.sin(alpha) * (1 - delta) * (1 - E1) / (1 - (1 - delta) * E1 * np.cos(alpha))
    return (1 - Smt / S0) * 100

delta_arr = np.linspace(0, 0.5, 300)

fig = make_subplots(rows=1, cols=2,
                    subplot_titles=['MTR vs MT Saturation (δ) for Different Tissues',
                                    'MTR vs Flip Angle (δ = 0.20)'])

# Left: MTR vs delta for three tissue types (different T1)
tissues = [('White Matter (T1=840 ms)',   840,  '#1f77b4'),
           ('Gray Matter (T1=1300 ms)',   1300, '#d62728'),
           ('CSF (T1=4500 ms)',           4500, '#2ca02c')]
TR = 30; alpha = 10

for name, T1, col in tissues:
    mtr = mtr_model(delta_arr, TR, T1, alpha)
    fig.add_trace(go.Scatter(x=delta_arr, y=mtr, name=name,
                             line=dict(color=col, width=2.2)), row=1, col=1)

# Right: MTR vs flip angle for fixed delta
alpha_arr = np.linspace(1, 30, 200)
delta_fixed = 0.20
for name, T1, col in tissues:
    mtr_fa = [mtr_model(delta_fixed, TR, T1, a) for a in alpha_arr]
    fig.add_trace(go.Scatter(x=alpha_arr, y=mtr_fa, name=name,
                             line=dict(color=col, width=2.2),
                             showlegend=False), row=1, col=2)

fig.update_xaxes(title_text='MT Saturation δ', row=1, col=1)
fig.update_xaxes(title_text='Flip Angle α (°)', row=1, col=2)
fig.update_yaxes(title_text='MTR (%)', row=1, col=1)
fig.update_yaxes(title_text='MTR (%)', row=1, col=2)
fig.update_layout(
    height=430,
    title_text=f'MTR Simulation  (TR = {TR} ms, α = {alpha}° left; δ = {delta_fixed} right)',
    legend=dict(orientation='h', yanchor='bottom', y=-0.28,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows how MTR increases with the MT saturation
parameter $\delta$. At a given $\delta$, WM shows higher MTR than GM
due to its shorter T1 — this demonstrates that MTR reflects both MT
physics and T1 relaxation. The right panel shows the dependence of
MTR on flip angle: MTR is not constant with $\alpha$, illustrating the
protocol-dependence that limits inter-site comparability.

## Data Fitting

MTR is computed pixel-by-pixel as a simple arithmetic operation. Here
we simulate a 1D profile through tissue with varying bound-pool fraction
(proportional to $\delta$) and reconstruct the MTR map from noisy data.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go

rng = np.random.default_rng(7)
N = 250; TR = 30; alpha_deg = 10; alpha = np.radians(alpha_deg); T1 = 1000.0

x = np.linspace(0, np.pi, N)
# True delta profile: high in WM-like region, low in GM-like region
delta_true = 0.06 + 0.22 * np.exp(-((x - 1.0) / 0.5)**2) + \
             0.10 * np.exp(-((x - 2.5) / 0.4)**2)

E1 = np.exp(-TR / T1)
noise = 0.008

S0_true  = np.sin(alpha) * (1 - E1) / (1 - E1 * np.cos(alpha))
Smt_true = np.sin(alpha) * (1 - delta_true) * (1 - E1) / \
           (1 - (1 - delta_true) * E1 * np.cos(alpha))

S0_n  = S0_true  + rng.normal(0, noise, N)
Smt_n = Smt_true + rng.normal(0, noise, N)

MTR_true     = (1 - Smt_true / S0_true) * 100
MTR_measured = np.clip((1 - Smt_n / S0_n) * 100, 0, 100)

pixel = np.arange(N)
rmse  = np.sqrt(np.mean((MTR_measured - MTR_true)**2))

fig = go.Figure()
fig.add_trace(go.Scatter(x=pixel, y=MTR_true,
                         name='True MTR',
                         line=dict(color='gray', dash='dot', width=1.5)))
fig.add_trace(go.Scatter(x=pixel, y=MTR_measured,
                         name='Measured MTR (noisy)',
                         line=dict(color='royalblue', width=2)))
fig.update_layout(
    title=f'MTR Map Reconstruction  (TR={TR} ms, α={alpha_deg}°, RMSE={rmse:.2f}%)',
    xaxis_title='Pixel index',
    yaxis_title='MTR (%)',
    template='plotly_white',
    height=360,
    legend=dict(x=0.98, y=0.98, xanchor='right', yanchor='top'),
)
fig.show()
```

The MTR map faithfully tracks the underlying $\delta$ profile. Because
MTR is a ratio, common-mode noise largely cancels, giving robust
estimates even at modest SNR. However, the absolute MTR values depend
on the chosen TR and flip angle, underscoring the need for standardised
protocols when comparing MTR across sites or time points.

## Example Brain Maps

Unsaturated and MT-saturated SPGR images are simulated and the
pixelwise MTR map is computed.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

rng = np.random.default_rng(9)

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
    return t, mb

T1v  ={0:0,1:0,3:1300,4:840, 5:1200,6:1100,7:4500}
PDv  ={0:0,1:0,3:0.85,4:0.75,5:0.80, 6:0.78,7:1.00}
MTRv ={0:0,1:0,3:35,  4:44,  5:40,   6:38,  7:2}    # percent
tissue, bmask = brain_phantom()
T1m   = np.vectorize(T1v.get)(tissue).astype(float)
PDm   = np.vectorize(PDv.get)(tissue).astype(float)
MTRtm = np.vectorize(MTRv.get)(tissue).astype(float) / 100.0

TR, fa_deg = 28.0, 15.0
fa  = np.radians(fa_deg)
E1m = np.exp(-TR / np.where(T1m > 0, T1m, 1))
S0_ref = PDm * np.sin(fa) * (1 - E1m) / (1 - E1m * np.cos(fa) + 1e-10)

def rician(img, s):
    re = img + rng.normal(0, s, img.shape)
    im = rng.normal(0, s, img.shape)
    return np.sqrt(re**2 + im**2)

img_s0  = rician(S0_ref * bmask, 0.020)
img_smt = rician(S0_ref * (1 - MTRtm) * bmask, 0.020)
MTR_map = np.where(bmask & (img_s0 > 0.01),
                   (1 - img_smt / img_s0) * 100, np.nan)

fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
axes[0].imshow(img_s0,  cmap='gray', vmin=0, vmax=0.5, interpolation='bilinear')
axes[0].set_title('S₀ (no saturation)', fontsize=10); axes[0].axis('off')
axes[1].imshow(img_smt, cmap='gray', vmin=0, vmax=0.5, interpolation='bilinear')
axes[1].set_title('S_MT (with saturation)', fontsize=10); axes[1].axis('off')
im = axes[2].imshow(MTR_map, cmap='viridis', vmin=0, vmax=60,
                    interpolation='bilinear')
axes[2].set_title('MTR map', fontsize=10); axes[2].axis('off')
div = make_axes_locatable(axes[2])
cax = div.append_axes('right', size='5%', pad=0.04)
plt.colorbar(im, cax=cax, label='MTR (%)')

fig.suptitle(f'MTR Brain Maps  (TR = {int(TR)} ms, α = {int(fa_deg)}°)',
             fontsize=11, y=1.02)
plt.tight_layout()
plt.show()
```

The MT-saturated image (S_MT) has lower overall signal than S₀, with
WM showing the greatest attenuation (highest MTR ≈ 44%) due to its
dense myelin-associated macromolecular pool. CSF has virtually no MT
effect (MTR ≈ 2%). The MTR map clearly delineates WM from GM and CSF,
making it sensitive to myelin-related pathology.

## References

```{bibliography}
:filter: docname in docnames
```
