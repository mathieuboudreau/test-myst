---
title: "3.3 MP2RAGE"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 3.3 MP2RAGE

MP2RAGE (magnetization-prepared 2 rapid acquisition gradient echoes)
is a 3D whole-brain T1 mapping technique introduced by Marques et al.
in 2010 {cite}`marques2010`. It combines a magnetization-prepared
inversion recovery with two GRE readout trains at different inversion
times. By forming a specific ratio of the two images, MP2RAGE eliminates
the receive-coil sensitivity profile and is substantially less sensitive
to B1 transmit inhomogeneity than VFA — making it the method of choice
for T1 mapping at 7 T.

## Pulse Sequence

Each MP2RAGE repetition consists of four stages:

1. **180° inversion pulse** — inverts all longitudinal magnetisation.
2. **Wait TI1** — partial inversion recovery.
3. **GRE1 readout block** — $N_1$ excitations at flip angle $\alpha_1$,
   acquiring image 1 centred at inversion time TI1.
4. **Wait TI2 − TI1** — further recovery.
5. **GRE2 readout block** — $N_2$ excitations at flip angle $\alpha_2$,
   acquiring image 2 centred at inversion time TI2.
6. **Wait TD** — relaxation back toward equilibrium before the next
   inversion.

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

def rect_pulse(t, t0, w, h=1.0):
    r = np.zeros_like(t); r[(t >= t0) & (t <= t0 + w)] = h; return r

fig, axes = plt.subplots(4, 1, figsize=(12, 5.5), sharex=True)
fig.subplots_adjust(hspace=0.06, left=0.07, right=0.97, top=0.88, bottom=0.14)
t = np.linspace(0, 1, 4000)

t_inv = 0.02
t_gre1_start = 0.14
t_gre2_start = 0.55
n_gre = 5; gre_width = 0.03; gre_gap = 0.038

gz_all = np.zeros_like(t); gy_all = np.zeros_like(t)
gx_all = np.zeros_like(t); adc_mask = np.zeros(len(t), dtype=bool)

inv = rect_pulse(t, t_inv, 0.04, 0.90)
axes[0].fill_between(t, 0, inv, color='tomato', alpha=0.82, label='180° inv.')
axes[0].annotate('180°', (t_inv + 0.02, 0.93), ha='center', fontsize=8, color='tomato')
gz_all += rect(t, t_inv - 0.005, 0.05, 0.72)

for block_idx, (t_start, blk_color, blk_label, alpha_lbl) in enumerate([
        (t_gre1_start, '#1f77b4', 'GRE₁ block (TI₁)', 'α₁'),
        (t_gre2_start, '#2ca02c', 'GRE₂ block (TI₂)', 'α₂')]):
    for k in range(n_gre):
        tp = t_start + k * gre_gap
        rf_k = sinc_pulse(t, tp, gre_width, 0.50)
        axes[0].fill_between(t, 0, rf_k, color=blk_color, alpha=0.80)
        gz_all += rect(t, tp - 0.005, gre_width + 0.01, 0.65)
        gz_all += rect(t, tp + gre_width + 0.005, 0.015, -0.32)
        if k % 2 == 0:
            gy_all += rect(t, tp + gre_width + 0.005, 0.02,
                           0.50 * (0.3 + 0.7 * k / (n_gre - 1)))
        te = tp + gre_width + 0.025
        gx_all += rect(t, tp + gre_width + 0.005, 0.02, -0.40)
        gx_all += rect(t, te, 0.04, 0.40)
        adc_mask |= ((t >= te) & (t <= te + 0.04))
    axes[0].text(t_start + 2.5 * gre_gap, 0.55, blk_label,
                 ha='center', fontsize=8, color=blk_color,
                 bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7))

gz_all += rect(t, 0.93, 0.05, 0.62)

axes[0].axhline(0, color='k', linewidth=0.7)
axes[0].set_ylabel('RF', fontsize=9); axes[0].set_ylim(-0.1, 1.12)
axes[0].yaxis.set_tick_params(labelleft=False)
axes[0].legend(loc='upper right', fontsize=8)
axes[0].set_title('MP2RAGE Pulse Sequence — One Full Inversion Recovery', fontsize=11, pad=6)

axes[1].fill_between(t, 0, gz_all, where=gz_all > 0, color='forestgreen', alpha=0.7)
axes[1].fill_between(t, 0, gz_all, where=gz_all < 0, color='tomato', alpha=0.7)
axes[1].axhline(0, color='k', linewidth=0.7)
axes[1].set_ylabel('Gz', fontsize=9); axes[1].yaxis.set_tick_params(labelleft=False)

axes[2].fill_between(t, 0, gy_all, color='darkorange', alpha=0.7)
axes[2].axhline(0, color='k', linewidth=0.7)
axes[2].set_ylabel('Gy', fontsize=9); axes[2].yaxis.set_tick_params(labelleft=False)

axes[3].fill_between(t, 0, gx_all, where=gx_all > 0, color='purple', alpha=0.7)
axes[3].fill_between(t, 0, gx_all, where=gx_all < 0, color='plum', alpha=0.7)
axes[3].fill_between(t, -0.15, -0.05, where=adc_mask, color='gold', alpha=0.9)
axes[3].axhline(0, color='k', linewidth=0.7)
axes[3].set_ylabel('Gx / ADC', fontsize=9); axes[3].yaxis.set_tick_params(labelleft=False)
axes[3].set_xlabel('Time (one TR_MP2RAGE)', fontsize=9)
axes[3].set_xlim(0, 1.0); axes[3].set_ylim(-0.25, 0.70)

ti1_c = t_gre1_start + 2 * gre_gap; ti2_c = t_gre2_start + 2 * gre_gap
axes[3].annotate('', xy=(ti1_c, -0.20), xytext=(t_inv, -0.20),
                 arrowprops=dict(arrowstyle='<->', color='#1f77b4', lw=1.2))
axes[3].text((t_inv + ti1_c)/2, -0.25, 'TI₁', ha='center', fontsize=9, color='#1f77b4')
axes[3].annotate('', xy=(ti2_c, -0.20), xytext=(t_inv, -0.20),
                 arrowprops=dict(arrowstyle='<->', color='#2ca02c', lw=1.2))
axes[3].text((t_inv + ti2_c)/2 + 0.12, -0.25, 'TI₂', ha='center', fontsize=9, color='#2ca02c')
axes[3].annotate('', xy=(0.98, -0.20), xytext=(ti2_c + 0.15, -0.20),
                 arrowprops=dict(arrowstyle='<->', color='gray', lw=1.0))
axes[3].text((ti2_c + 0.15 + 0.98)/2, -0.25, 'TD', ha='center', fontsize=9, color='gray')

plt.close()
```

:::{figure} mp2rage_pulse_sequence.png
:name: fig-mp2rage-sequence
:align: center
MP2RAGE pulse sequence. After a 180° inversion, two GRE readout blocks
(blue: TI1, green: TI2) are interspersed within the inversion recovery.
TI1 and TI2 are chosen so that the two images straddle the zero-crossing
of the recovery curve. The recovery delay TD allows the magnetisation to
return toward equilibrium before the next inversion.
:::

## Mathematical Model

### Magnetisation Evolution

In the simplified model where each readout block is approximated as a
single instantaneous flip, the free-pool longitudinal magnetisation at
the centre of each readout is {cite}`marques2010`:

$$
\text{GRE}_1 \propto M_0\!\left(1 - 2\,e^{-\text{TI}_1/T_1}\right)\sin\alpha_1
$$ (eq-mp2rage-gre1)

$$
\text{GRE}_2 \propto M_0\!\left(1 - 2\,e^{-\text{TI}_2/T_1}\right)\sin\alpha_2
$$ (eq-mp2rage-gre2)

Here TI1 < TI2 are chosen such that $1 - 2e^{-\text{TI}_1/T_1} < 0$ and
$1 - 2e^{-\text{TI}_2/T_1} > 0$ for the tissue of interest — i.e., TI1
is before the null crossing and TI2 is after it.

### The MP2RAGE Contrast

The key insight is to form the **complex combination** of the two images
{cite}`marques2010`:

$$
\text{MP2RAGE} =
\frac{\operatorname{Re}\!\left(\text{GRE}_1^* \cdot \text{GRE}_2\right)}
     {|\text{GRE}_1|^2 + |\text{GRE}_2|^2}
$$ (eq-mp2rage-contrast)

For the simplified two-image model with real-valued signals this reduces to:

$$
\text{MP2RAGE}(T_1) =
\frac{\sin\alpha_1\sin\alpha_2\,(1 - 2E_1)(1 - 2E_2)}
     {\sin^2\!\alpha_1\,(1-2E_1)^2 + \sin^2\!\alpha_2\,(1-2E_2)^2}
$$ (eq-mp2rage-ratio)

where $E_i = e^{-\text{TI}_i/T_1}$.

Several important properties follow:

- **Receive-field cancellation**: $M_0$ and the coil sensitivity factor
  appear in both numerator and denominator and cancel exactly.
- **B1 robustness**: the ratio depends on $\alpha_1$ and $\alpha_2$ only
  weakly through the $\sin\alpha$ factors; for uniform $\alpha$ the
  dependence on $B_1$ nearly cancels.
- **Monotonic T1 mapping**: Eq. [%s](eq-mp2rage-ratio) is a monotonic
  function of T1, enabling a look-up table (LUT) inversion to obtain T1.

The MP2RAGE ratio ranges from $-0.5$ (short T1, both images same sign)
to $+0.5$ (long T1, images opposite sign), with 0 at the T1 value where
$1 - 2e^{-\text{TI}_1/T_1} = 0$ (the null tissue).

:::{note}
The full Marques model accounts for the perturbation of the inversion
recovery curve by the $N_1$ readout pulses in the GRE1 block, which
modifies the effective magnetisation presented to GRE2. Here the simplified
single-excitation model is used for conceptual clarity; the qualitative
T1-dependence and M0-cancellation properties are preserved.
:::

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def mp2rage_contrast(T1, TI1, TI2, alpha1_deg, alpha2_deg):
    """Simplified MP2RAGE ratio (single-excitation model)."""
    E1 = np.exp(-TI1 / T1); E2 = np.exp(-TI2 / T1)
    a1 = np.sin(np.radians(alpha1_deg)); a2 = np.sin(np.radians(alpha2_deg))
    g1 = a1 * (1 - 2 * E1); g2 = a2 * (1 - 2 * E2)
    num = g1 * g2; den = g1**2 + g2**2
    return np.where(den > 1e-12, num / den, 0.0)

# Standard 3T protocol (Marques 2010)
TI1 = 800; TI2 = 2700; alpha1 = 4; alpha2 = 5   # ms, degrees

T1_arr = np.linspace(200, 4500, 800)
mp2r = mp2rage_contrast(T1_arr, TI1, TI2, alpha1, alpha2)

fig = make_subplots(rows=1, cols=2,
                    subplot_titles=['MP2RAGE Contrast vs T1 (LUT)',
                                    'Effect of TI Choice on T1 Sensitivity'])

fig.add_trace(go.Scatter(x=T1_arr, y=mp2r,
                         name=f'TI1={TI1} ms, TI2={TI2} ms',
                         line=dict(color='royalblue', width=2.5)), row=1, col=1)
fig.add_hline(y=0, line_dash='dot', line_color='gray', row=1, col=1)

# Mark brain tissue T1 values
for name, T1_v, col in [('WM', 840, '#1f77b4'),
                         ('GM', 1300, '#d62728'),
                         ('CSF', 4500, '#2ca02c')]:
    mp2_v = mp2rage_contrast(T1_v, TI1, TI2, alpha1, alpha2)
    fig.add_trace(go.Scatter(x=[T1_v], y=[mp2_v], mode='markers',
                             name=name,
                             marker=dict(color=col, size=10, symbol='diamond')),
                 row=1, col=1)

# Right: different TI pairs show trade-off
ti_pairs = [(700, 2500, 'TI1=700, TI2=2500'),
            (800, 2700, 'TI1=800, TI2=2700 (default)'),
            (900, 3000, 'TI1=900, TI2=3000')]
cols_ti = ['#aec7e8', '#1f77b4', '#08519c']

for (t1, t2, lbl), cv in zip(ti_pairs, cols_ti):
    mp2 = mp2rage_contrast(T1_arr, t1, t2, alpha1, alpha2)
    fig.add_trace(go.Scatter(x=T1_arr, y=mp2, name=lbl,
                             line=dict(color=cv, width=2)), row=1, col=2)

fig.add_hline(y=0, line_dash='dot', line_color='gray', row=1, col=2)

for col in [1, 2]:
    fig.update_xaxes(title_text='T1 (ms)', row=1, col=col)
fig.update_yaxes(title_text='MP2RAGE contrast', row=1, col=1)
fig.update_yaxes(title_text='MP2RAGE contrast', row=1, col=2)

fig.update_layout(
    height=430,
    title_text=f'MP2RAGE Simulation  (α₁={alpha1}°, α₂={alpha2}°)',
    legend=dict(orientation='h', yanchor='bottom', y=-0.32,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows the MP2RAGE look-up table relating the ratio to T1.
The contrast is monotonic and spans $[-0.5, +0.5]$, with the zero
crossing near the null T1 determined by TI1. Brain tissue T1 values
(WM, GM, CSF) map to distinct, well-separated contrast levels. The right
panel shows how varying the TI pair shifts and steepens the LUT, allowing
the sensitivity range to be tuned to the expected tissue T1 distribution.

## Data Fitting

T1 is extracted by inverting the look-up table: for each voxel's
measured MP2RAGE ratio, find the T1 that produces that ratio via the
forward model.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go

rng = np.random.default_rng(42)

TI1 = 800; TI2 = 2700; alpha1 = 4; alpha2 = 5

# Build LUT
T1_lut = np.linspace(100, 5000, 5000)
mp2_lut = mp2rage_contrast(T1_lut, TI1, TI2, alpha1, alpha2)
# Make LUT monotone for inversion (keep only increasing portion ~[-0.5, 0.5])
valid = (mp2_lut > -0.499) & (mp2_lut < 0.499)
T1_lut_v  = T1_lut[valid]
mp2_lut_v = mp2_lut[valid]

def lut_t1(mp2_vals):
    return np.interp(mp2_vals, mp2_lut_v, T1_lut_v)

# Simulate 1D profile
N = 250
x = np.linspace(0, np.pi, N)
T1_true = 700 + 700 * np.sin(x / 2)**2   # 700–1400 ms
M0_true = 1.0 + 0.3 * np.sin(x)

noise_frac = 0.04   # 4 % noise in each image

GRE1_true = M0_true * np.sin(np.radians(alpha1)) * (1 - 2 * np.exp(-TI1 / T1_true))
GRE2_true = M0_true * np.sin(np.radians(alpha2)) * (1 - 2 * np.exp(-TI2 / T1_true))

ref = np.abs(GRE1_true).mean()
GRE1_n = GRE1_true + rng.normal(0, noise_frac * ref, N)
GRE2_n = GRE2_true + rng.normal(0, noise_frac * ref, N)

num_n = GRE1_n * GRE2_n; den_n = GRE1_n**2 + GRE2_n**2
mp2_meas = np.where(den_n > 1e-12, num_n / den_n, 0.0)
mp2_meas = np.clip(mp2_meas, mp2_lut_v.min() + 1e-4, mp2_lut_v.max() - 1e-4)

T1_est = lut_t1(mp2_meas)
rmse = np.sqrt(np.mean((T1_est - T1_true)**2))

pixel = np.arange(N)
fig = go.Figure()
fig.add_trace(go.Scatter(x=pixel, y=T1_true,
                         name='True T1 (ms)',
                         line=dict(color='gray', dash='dot', width=1.5)))
fig.add_trace(go.Scatter(x=pixel, y=T1_est,
                         name=f'MP2RAGE T1 estimate (RMSE={rmse:.0f} ms)',
                         line=dict(color='royalblue', width=2)))

fig.update_layout(
    title=(f'MP2RAGE T1 Map Reconstruction  '
           f'(TI1={TI1} ms, TI2={TI2} ms, {noise_frac*100:.0f}% noise)'),
    xaxis_title='Pixel index',
    yaxis_title='T1 (ms)',
    template='plotly_white',
    height=380,
    legend=dict(x=0.01, y=0.99, xanchor='left', yanchor='top'),
)
fig.show()
```

The MP2RAGE LUT inversion recovers T1 accurately across the full range
of tissue values. Noise in the source images propagates through the
non-linear ratio, but because the ratio suppresses M0, the T1 estimate
is insensitive to receive-coil intensity variations and the method is
robust to B1 inhomogeneity — the primary advantages over IR and VFA.

## Example Brain Maps

An MP2RAGE simulation shows the two GRE images, the M0-insensitive
contrast ratio, and the final T1 map obtained by LUT inversion.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

rng = np.random.default_rng(5)

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

T1v={0:0,1:0,3:1300,4:840, 5:1200,6:1100,7:4500}
PDv={0:0,1:0,3:0.85,4:0.75,5:0.80,6:0.78,7:1.00}
tissue, bmask = brain_phantom()
T1m = np.vectorize(T1v.get)(tissue).astype(float)
PDm = np.vectorize(PDv.get)(tissue).astype(float)
n   = tissue.shape[0]

TI1, TI2, a1_deg, a2_deg = 800, 2700, 4, 5
a1, a2 = np.radians(a1_deg), np.radians(a2_deg)

T1s = np.where(T1m > 0, T1m, 1)
g1 = PDm * np.sin(a1) * (1 - 2*np.exp(-TI1 / T1s))
g2 = PDm * np.sin(a2) * (1 - 2*np.exp(-TI2 / T1s))

def rician(img, s):
    re = img + rng.normal(0, s, img.shape)
    im = rng.normal(0, s, img.shape)
    return np.sqrt(re**2 + im**2)

gre1 = rician(np.abs(g1) * bmask, 0.02)
gre2 = rician(np.abs(g2) * bmask, 0.02)

# Phase-restored combination
g1_signed = gre1 * np.sign(g1)
num = g1_signed * gre2
den = gre1**2 + gre2**2
mp2_map = np.where(bmask & (den > 1e-9), num / den, 0.0)

# LUT inversion
T1_lut  = np.linspace(200, 5000, 3000)
mp2_lut = (np.sin(a1) * np.sin(a2)
           * (1 - 2*np.exp(-TI1/T1_lut)) * (1 - 2*np.exp(-TI2/T1_lut)))
den_lut = (np.sin(a1)**2 * (1 - 2*np.exp(-TI1/T1_lut))**2
         + np.sin(a2)**2 * (1 - 2*np.exp(-TI2/T1_lut))**2)
mp2_lut = np.where(den_lut > 1e-9, mp2_lut / den_lut, 0.0)
T1_mp2  = np.where(bmask, np.interp(mp2_map, mp2_lut, T1_lut), np.nan)

fig, axes = plt.subplots(1, 4, figsize=(14, 3.8))
axes[0].imshow(gre1, cmap='gray', vmin=0, vmax=0.5, interpolation='bilinear')
axes[0].set_title(f'GRE₁  (TI = {TI1} ms)', fontsize=9); axes[0].axis('off')
axes[1].imshow(gre2, cmap='gray', vmin=0, vmax=0.5, interpolation='bilinear')
axes[1].set_title(f'GRE₂  (TI = {TI2} ms)', fontsize=9); axes[1].axis('off')

im_mp2 = axes[2].imshow(np.where(bmask, mp2_map, np.nan),
                         cmap='RdBu_r', vmin=-0.5, vmax=0.5,
                         interpolation='bilinear')
axes[2].set_title('MP2RAGE contrast', fontsize=9); axes[2].axis('off')
div2 = make_axes_locatable(axes[2])
cax2 = div2.append_axes('right', size='5%', pad=0.04)
plt.colorbar(im_mp2, cax=cax2)

im_t1 = axes[3].imshow(T1_mp2, cmap='magma', vmin=400, vmax=4500,
                        interpolation='bilinear')
axes[3].set_title('MP2RAGE T1 map', fontsize=9); axes[3].axis('off')
div3 = make_axes_locatable(axes[3])
cax3 = div3.append_axes('right', size='5%', pad=0.04)
plt.colorbar(im_t1, cax=cax3, label='T1 (ms)')

fig.suptitle(f'MP2RAGE Brain Maps  (TI₁={TI1} ms, TI₂={TI2} ms,'
             f' α₁={a1_deg}°, α₂={a2_deg}°)', fontsize=11, y=1.02)
plt.tight_layout()
plt.show()
```

The GRE₁ image (TI = 800 ms) shows tissues that have passed their
null point as bright (GM, CSF) while WM near its null appears dark.
GRE₂ (TI = 2700 ms) is closer to proton-density-weighted. The MP2RAGE
contrast ratio (third panel) is negative for WM (still inverted at
TI₁) and positive for GM and CSF. After LUT inversion, the T1 map
recovers the correct tissue values without any receive-coil sensitivity
bias.

## References

```{bibliography}
:filter: docname in docnames
```
