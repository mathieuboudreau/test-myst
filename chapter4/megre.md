---
title: "4.2 Multi-Echo Gradient Echo"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 4.2 Multi-Echo Gradient Echo

The multi-echo gradient echo (ME-GRE) sequence measures the apparent
transverse relaxation rate R2* = 1/T2*, which includes both irreversible
T2 dephasing and reversible dephasing from static B0 inhomogeneity.
A single RF excitation is followed by a train of gradient echoes at
increasing echo times. ME-GRE data are the foundation for R2* mapping,
BOLD fMRI, and quantitative susceptibility mapping (QSM).

## Pulse Sequence

After a single slice-selective excitation at flip angle $\alpha$:

1. **Pre-phasing lobe** — dephases the FID in the readout direction.
2. **Alternating readout gradients** — each reversal creates a gradient
   echo at TE₁, TE₂, TE₃, …
3. **Phase encode** — one step applied before the first echo; all
   subsequent echoes share the same k-space line.

Unlike spin echoes, the gradient echoes do **not** refocus $\Delta B_0$
dephasing accumulated between echoes, so the signal decays with T2*.

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

fig, axes = plt.subplots(4, 1, figsize=(12, 5.2), sharex=True)
fig.subplots_adjust(hspace=0.06, left=0.07, right=0.97, top=0.87, bottom=0.14)
t = np.linspace(0, 1, 5000)

N_gre = 4; exc_t = 0.03; exc_w = 0.06; exc_ht = 0.40
TE1 = 0.16; dTE = 0.17

gz_gre = rect(t, exc_t-0.01, exc_w+0.02, 0.65) + rect(t, exc_t+exc_w+0.01, 0.03, -0.32)
gy_gre = rect(t, exc_t+exc_w+0.04, 0.03, 0.52)
gx_gre = rect(t, exc_t+exc_w+0.01, 0.07, -0.50)
adc_gre = np.zeros(len(t), dtype=bool)

echo_tes = [TE1 + k * dTE for k in range(N_gre)]
for k, te in enumerate(echo_tes):
    sign = (-1)**k
    gx_gre  += rect(t, te-0.07, 0.14, sign * 0.50)
    adc_gre |= (t >= te-0.03) & (t <= te+0.03)

axes[0].fill_between(t, 0, sinc_pulse(t, exc_t, exc_w, exc_ht),
                     color='cornflowerblue', alpha=0.9, label='α excitation')
axes[0].axhline(0, color='k', lw=0.7)
axes[0].set_ylabel('RF', fontsize=9); axes[0].yaxis.set_tick_params(labelleft=False)
axes[0].set_ylim(-0.1, 0.55); axes[0].legend(loc='upper right', fontsize=7.5)
axes[0].set_title('Multi-Echo Gradient Echo (ME-GRE) Pulse Sequence', fontsize=11, pad=5)

axes[1].fill_between(t, 0, gz_gre, where=gz_gre>0, color='forestgreen', alpha=0.7)
axes[1].fill_between(t, 0, gz_gre, where=gz_gre<0, color='tomato', alpha=0.7)
axes[1].axhline(0, color='k', lw=0.7); axes[1].set_ylabel('Gz', fontsize=9)
axes[1].yaxis.set_tick_params(labelleft=False)

axes[2].fill_between(t, 0, gy_gre, color='darkorange', alpha=0.7)
axes[2].axhline(0, color='k', lw=0.7); axes[2].set_ylabel('Gy', fontsize=9)
axes[2].yaxis.set_tick_params(labelleft=False)

axes[3].fill_between(t, 0, gx_gre, where=gx_gre>0, color='purple', alpha=0.7)
axes[3].fill_between(t, 0, gx_gre, where=gx_gre<0, color='plum', alpha=0.7)
axes[3].fill_between(t, -0.18, -0.06, where=adc_gre, color='gold', alpha=0.9)
axes[3].axhline(0, color='k', lw=0.7)
axes[3].set_ylabel('Gx / ADC', fontsize=9); axes[3].yaxis.set_tick_params(labelleft=False)
axes[3].set_xlabel('Time (one TR)', fontsize=9)
axes[3].set_xlim(0, 1); axes[3].set_ylim(-0.28, 0.68)

exc_centre = exc_t + exc_w/2
for k, te in enumerate(echo_tes):
    axes[3].annotate('', xy=(te, -0.23), xytext=(exc_centre, -0.23),
                     arrowprops=dict(arrowstyle='<->', color=f'C{k}', lw=1.0))
    axes[3].text((exc_centre + te)/2, -0.27, f'TE{"₁₂₃₄"[k]}',
                 ha='center', fontsize=8.5, color=f'C{k}')
plt.close()
```

:::{figure} megre_pulse_sequence.png
:name: fig-megre-sequence
:align: center
ME-GRE pulse sequence. A single RF excitation is followed by four
gradient echoes formed by reversing the frequency-encode gradient.
Because no 180° refocusing is applied, signal decays with T2*
rather than T2.
:::

## Mathematical Model

### T2* Decay

The complex ME-GRE signal at echo time TE is:

$$
S(\text{TE}) = S_0\,e^{-\text{TE}/T_2^*}\,e^{i\,2\pi\,\Delta f\,\text{TE}}
$$ (eq-megre-complex)

where $\Delta f = \gamma\,\Delta B_0/(2\pi)$ is the local off-resonance
frequency (Hz) and $S_0$ depends on T1, flip angle, and TR through the
Ernst equation. The magnitude:

$$
|S(\text{TE})| = S_0\,e^{-\text{TE}/T_2^*}
$$ (eq-megre-mag)

### Relationship Between T2, T2', and T2*

The effective relaxation rate R2* decomposes as:

$$
R_2^* = R_2 + R_2', \qquad R_2' = \gamma\,|\Delta B_0|
$$ (eq-megre-r2star)

$R_2'$ (pronounced "R2 prime") is the reversible dephasing from $\Delta B_0$.
A spin echo refocuses this contribution; a gradient echo does not. Measuring
both T2 (spin echo) and T2* (gradient echo) allows R2' to be computed and,
through $\Delta B_0$, provides a local susceptibility estimate.

### R2* Mapping

R2* maps are obtained by fitting [Eq. %s](eq-megre-mag) at every voxel.
The linear log-domain fit is:

$$
\ln|S(\text{TE})| = \ln S_0 - R_2^*\,\text{TE}
$$

giving R2* as the slope. A voxel-by-voxel linear regression is fast and
sufficient for most R2* mapping applications. Alternatively, the full
complex signal [Eq. %s](eq-megre-complex) can be fitted simultaneously
for R2* and $\Delta f$, which is required for QSM.

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

TE_arr = np.linspace(2, 80, 400)

fig = make_subplots(rows=1, cols=2,
                    subplot_titles=['T2* Decay Curves by Tissue',
                                    'T2 vs T2* — Effect of B0 Inhomogeneity'])

tissues = [('WM  T2*=35 ms',  35,  '#1f77b4'),
           ('GM  T2*=45 ms',  45,  '#d62728'),
           ('CSF T2*=150 ms', 150, '#2ca02c')]

for name, T2s, col in tissues:
    fig.add_trace(go.Scatter(x=TE_arr, y=np.exp(-TE_arr/T2s),
                             name=name, line=dict(color=col, width=2.2)), row=1, col=1)

# Right: show R2* = R2 + R2' for varying B0 deviation
T2_wm = 75.0   # ms, true T2 of WM
dB0_arr = np.linspace(0, 5, 200)   # uT – gives dF in Hz at 3T (γ=267.5e6 rad/s/T)
gamma_over_2pi = 42.577e6  # Hz/T
dF_arr = gamma_over_2pi * dB0_arr * 1e-6  # Hz
R2_wm = 1000.0 / T2_wm   # Hz (1/s using ms→ s conversion via 1/T2*1000)
# Switch to consistent units: ms
R2_wm_inv_ms = 1.0 / T2_wm
R2p_inv_ms   = (2 * np.pi * dF_arr) / 1000.0   # 1/ms, from rad/s to 1/ms
R2star_inv_ms = R2_wm_inv_ms + R2p_inv_ms
T2star_arr   = 1.0 / R2star_inv_ms   # ms

fig.add_trace(go.Scatter(x=dB0_arr, y=T2star_arr,
                         name='T2* (WM)',
                         line=dict(color='royalblue', width=2.5)), row=1, col=2)
fig.add_hline(y=T2_wm, line_dash='dot', line_color='gray',
              annotation_text=f'T2={T2_wm:.0f} ms',
              annotation_position='right', row=1, col=2)

fig.update_xaxes(title_text='TE (ms)', row=1, col=1)
fig.update_xaxes(title_text='ΔB0 (μT)', row=1, col=2)
fig.update_yaxes(title_text='Signal (a.u.)', row=1, col=1)
fig.update_yaxes(title_text='T2* (ms)', row=1, col=2)
fig.update_layout(
    height=430,
    title_text='ME-GRE / T2* Simulation',
    legend=dict(orientation='h', yanchor='bottom', y=-0.28,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows T2* decay curves. WM has the shortest T2* due to
susceptibility effects from myelin and iron. The right panel shows how
T2* of WM decreases monotonically as B0 inhomogeneity ($\Delta B_0$)
increases; the dotted line is the true T2 floor that T2* approaches
from above.

## Data Fitting

```{code-cell} python
:tags: [hide-input]

import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objects as go

rng = np.random.default_rng(11)

TE_pts = np.array([3, 6, 10, 15, 20, 30, 45, 60])   # ms
T2s_true = 38.0; S0_true = 1.0; dF_true = 8.0   # Hz
noise = 0.025

# Complex signal (magnitude data here)
S_complex = S0_true * np.exp(-TE_pts/T2s_true) * np.exp(1j * 2*np.pi * dF_true * TE_pts * 1e-3)
S_mag = np.abs(S_complex) + rng.normal(0, noise, len(TE_pts))
S_mag = np.maximum(S_mag, 0)

# 1. Log-linear on magnitude
log_S = np.log(np.maximum(S_mag, 1e-6))
coeffs = np.polyfit(TE_pts, log_S, 1)
T2s_lin = -1.0 / coeffs[0]; S0_lin = np.exp(coeffs[1])

# 2. NLS magnitude
def mag_model(TE, S0, T2s):
    return S0 * np.exp(-TE / T2s)
popt, pcov = curve_fit(mag_model, TE_pts, S_mag, p0=[1.0, 30.0],
                       bounds=([0, 1], [5, 500]), maxfev=5000)
S0_nl, T2s_nl = popt; T2s_err = np.sqrt(np.diag(pcov))[1]

TE_dense = np.linspace(1, 70, 400)

fig = go.Figure()
fig.add_trace(go.Scatter(x=TE_dense, y=S0_true * np.exp(-TE_dense/T2s_true),
                         name='Ground truth', mode='lines',
                         line=dict(color='gray', dash='dot', width=1.5)))
fig.add_trace(go.Scatter(x=TE_pts, y=S_mag, name='Noisy magnitude data',
                         mode='markers', marker=dict(color='royalblue', size=9)))
fig.add_trace(go.Scatter(x=TE_dense, y=S0_lin * np.exp(-TE_dense/T2s_lin),
                         name=f'Log-linear  T2*={T2s_lin:.0f} ms',
                         line=dict(color='darkorange', dash='dash', width=2)))
fig.add_trace(go.Scatter(x=TE_dense, y=mag_model(TE_dense, *popt),
                         name=f'NLS fit  T2*={T2s_nl:.0f} ± {T2s_err:.1f} ms',
                         line=dict(color='tomato', width=2.2)))

fig.update_layout(
    title=f'ME-GRE T2* Fitting  (True T2* = {T2s_true:.0f} ms)',
    xaxis_title='TE (ms)', yaxis_title='Signal (a.u.)',
    template='plotly_white', height=380,
    legend=dict(x=0.50, y=0.95),
)
fig.show()
```

Both fitting approaches recover T2* accurately. The key distinction
from MESE fitting is that early echoes here are critical: with T2* ≈
35–45 ms for brain parenchyma at 3 T, a significant fraction of the
signal is lost already by TE = 20 ms, so the earliest echo times
carry the most T2* information.

## Example Brain Maps

A six-echo ME-GRE acquisition is simulated and a vectorised log-linear
fit recovers the T2* map (displayed alongside the derived R2* map).

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

rng = np.random.default_rng(7)

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
T2sv={0:0,1:0,3:45,4:35,5:40,6:38,7:150}
PDv={0:0,1:0,3:0.85,4:0.75,5:0.80,6:0.78,7:1.00}
tissue, bmask = brain_phantom()
T1m  = np.vectorize(T1v.get)(tissue).astype(float)
T2sm = np.vectorize(T2sv.get)(tissue).astype(float)
PDm  = np.vectorize(PDv.get)(tissue).astype(float)
n    = tissue.shape[0]

TE_pts = np.array([3., 8., 15., 25., 40., 60.])
TR_gre, fa_deg = 20.0, 15.0
E1m = np.exp(-TR_gre / np.where(T1m > 0, T1m, 1))
fa  = np.radians(fa_deg)
S0  = PDm * np.sin(fa) * (1 - E1m) / (1 - E1m * np.cos(fa) + 1e-10)

def rician(img, s):
    re = img + rng.normal(0, s, img.shape)
    im = rng.normal(0, s, img.shape)
    return np.sqrt(re**2 + im**2)

imgs = [rician(S0 * np.exp(-te / np.where(T2sm > 0, T2sm, 1)) * bmask, 0.020)
        for te in TE_pts]

# Vectorised log-linear R2* fit
S_log  = np.log(np.maximum(np.stack(imgs).reshape(len(TE_pts), -1), 1e-6))
X      = np.column_stack([np.ones(len(TE_pts)), TE_pts])
Xpinv  = np.linalg.pinv(X)
coeffs = Xpinv @ S_log
R2s_fit = np.where(bmask, np.clip(-coeffs[1].reshape(n, n), 0.001, 0.5), np.nan)
T2s_fit = np.where(bmask, np.clip(1.0 / R2s_fit, 5, 200), np.nan)

show_idx = [0, 1, 3, 5]
fig, axes = plt.subplots(1, 5, figsize=(17, 3.8))
for ax, k in zip(axes[:4], show_idx):
    ax.imshow(imgs[k], cmap='gray', vmin=0, vmax=0.40, interpolation='bilinear')
    ax.set_title(f'TE = {int(TE_pts[k])} ms', fontsize=9); ax.axis('off')

im = axes[4].imshow(T2s_fit, cmap='plasma', vmin=5, vmax=80,
                    interpolation='bilinear')
axes[4].set_title('ME-GRE T2* map', fontsize=9); axes[4].axis('off')
div = make_axes_locatable(axes[4])
cax = div.append_axes('right', size='5%', pad=0.04)
plt.colorbar(im, cax=cax, label='T2* (ms)')

fig.suptitle(f'ME-GRE Brain Maps  (TR = {int(TR_gre)} ms, α = {int(fa_deg)}°)',
             fontsize=11, y=1.02)
plt.tight_layout()
plt.show()
```

The ME-GRE images lose signal much faster than MESE: by TE = 60 ms
all parenchyma is near-zero, while CSF (long T2*) retains signal.
The T2* map shows the expected tissue-contrast hierarchy: WM has
the shortest T2* (~35 ms) owing to susceptibility effects from
myelin, while CSF has the longest (~150 ms). Basal ganglia appear
slightly darker than surrounding WM in the T2* map due to higher
iron content.

## References

```{bibliography}
:filter: docname in docnames
```
