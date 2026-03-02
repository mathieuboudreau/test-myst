---
title: "3.1 Inversion Recovery"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 3.1 Inversion Recovery

Inversion recovery (IR) is the gold-standard T1 mapping technique {cite}`haacke1999`.
A 180° inversion pulse flips the longitudinal magnetisation to $-M_0$; the
free-pool magnetisation then recovers exponentially toward $+M_0$. By acquiring
the signal at multiple inversion times (TI) and fitting an exponential, T1 is
determined directly without relying on steady-state assumptions.

## Pulse Sequence

The IR sequence consists of:

1. **180° inversion pulse** — inverts $M_z$ from $+M_0$ to $-M_0$.
2. **Inversion time TI** — the magnetisation recovers during this delay.
3. **Readout** — a spoiled GRE (or spin-echo) acquisition that samples
   the recovered $M_z$ via a flip angle $\alpha$.
4. **Recovery delay** — TR $-$ TI allows partial re-equilibration before
   the next inversion.

The experiment is repeated $N$ times with different TI values, typically
logarithmically spaced from ~50 ms to 3–4 × T1 (estimated).

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

def rect_pulse(t, t0, w, h=1.0):
    r = np.zeros_like(t); r[(t >= t0) & (t <= t0 + w)] = h; return r

fig = plt.figure(figsize=(11, 5))
gs = gridspec.GridSpec(4, 2, hspace=0.06, wspace=0.35,
                       left=0.07, right=0.97, top=0.88, bottom=0.13)
t = np.linspace(0, 1, 2000)

for col_idx, (label, ti_frac, color) in enumerate([
        ('Short TI  (near null-crossing)', 0.22, '#1f77b4'),
        ('Long TI  (near full recovery)',  0.62, '#d62728')]):

    axes = [fig.add_subplot(gs[i, col_idx]) for i in range(4)]
    plt.setp([a.get_xticklabels() for a in axes[:-1]], visible=False)

    rf_inv  = rect_pulse(t, 0.03, 0.06, 0.85)
    rf_read = sinc_pulse(t, ti_frac, 0.07, 0.55)
    axes[0].fill_between(t, 0, rf_inv,  color='tomato', alpha=0.80, label='180° inv.')
    axes[0].fill_between(t, 0, rf_read, color=color,    alpha=0.85, label='α readout')
    axes[0].axhline(0, color='k', linewidth=0.7)
    axes[0].annotate('180°', (0.06, 0.88), ha='center', fontsize=8, color='tomato')
    axes[0].annotate('α',   (ti_frac + 0.035, 0.58), ha='center', fontsize=9)
    axes[0].set_ylabel('RF', fontsize=8); axes[0].set_ylim(-0.1, 1.05)
    axes[0].yaxis.set_tick_params(labelleft=False)
    axes[0].set_title(label, fontsize=9, pad=4)
    axes[0].legend(loc='upper right', fontsize=7, framealpha=0.8)

    gz = rect(t, 0.02, 0.08, 0.70) + rect(t, ti_frac - 0.01, 0.09, 0.70) + \
         rect(t, ti_frac + 0.08, 0.035, -0.35)
    axes[1].fill_between(t, 0, gz, where=gz > 0, color='forestgreen', alpha=0.7)
    axes[1].fill_between(t, 0, gz, where=gz < 0, color='tomato', alpha=0.7)
    axes[1].axhline(0, color='k', linewidth=0.7)
    axes[1].set_ylabel('Gz', fontsize=8); axes[1].yaxis.set_tick_params(labelleft=False)

    te_loc = ti_frac + 0.14
    gy = rect(t, te_loc - 0.03, 0.05, 0.55)
    axes[2].fill_between(t, 0, gy, color='darkorange', alpha=0.7)
    axes[2].axhline(0, color='k', linewidth=0.7)
    axes[2].set_ylabel('Gy', fontsize=8); axes[2].yaxis.set_tick_params(labelleft=False)

    gx = rect(t, te_loc - 0.03, 0.05, -0.42) + rect(t, te_loc + 0.03, 0.12, 0.42) + \
         rect(t, 0.88, 0.08, 0.55)
    axes[3].fill_between(t, 0, gx, where=gx > 0, color='purple', alpha=0.7)
    axes[3].fill_between(t, 0, gx, where=gx < 0, color='plum', alpha=0.7)
    axes[3].fill_between(t, -0.15, -0.05,
                         where=((t >= te_loc+0.03) & (t <= te_loc+0.15)),
                         color='gold', alpha=0.9)
    axes[3].axhline(0, color='k', linewidth=0.7)
    axes[3].set_ylabel('Gx / ADC', fontsize=8); axes[3].yaxis.set_tick_params(labelleft=False)
    axes[3].set_xlabel('Time (one TR)', fontsize=8)
    axes[3].set_xlim(0, 1.0); axes[3].set_ylim(-0.25, 0.75)
    axes[3].annotate('', xy=(ti_frac, -0.20), xytext=(0.03, -0.20),
                     arrowprops=dict(arrowstyle='<->', color='k', lw=0.8))
    axes[3].text((0.03 + ti_frac) / 2, -0.25, 'TI', ha='center', fontsize=9)

fig.suptitle('Inversion Recovery T1 Mapping — Two Representative TI Values',
             fontsize=11, y=0.97)
plt.close()
```

:::{figure} ir_pulse_sequence.png
:name: fig-ir-sequence
:align: center
IR pulse sequence. Left: short TI — the magnetisation has barely begun
to recover from the inversion and the readout samples $M_z \approx 0$.
Right: long TI — $M_z$ has nearly recovered to $+M_0$. By sampling
multiple TI values the full exponential recovery curve is mapped.
:::

## Mathematical Model

### Ideal Inversion Recovery

After a perfect 180° inversion pulse, the longitudinal magnetisation at
time TI is:

$$
M_z(\text{TI}) = M_0\!\left(1 - 2\,e^{-\text{TI}/T_1}\right)
$$ (eq-ir-ideal)

The factor 2 reflects the full inversion from $+M_0$ to $-M_0$.

For a spoiled GRE readout with flip angle $\alpha$, the acquired signal
amplitude is:

$$
S(\text{TI}) = M_0 \sin\alpha\,\bigl|1 - 2\,e^{-\text{TI}/T_1}\bigr|
$$ (eq-ir-magnitude)

where the absolute value appears for magnitude images. Phase-sensitive
reconstruction retains the sign of $M_z$, enabling the full dynamic
range from $-M_0$ to $+M_0$ and improving fitting precision.

### Three-Parameter Model

In practice TR is finite and the inversion may be imperfect. The
general three-parameter model is {cite}`haacke1999`:

$$
S(\text{TI}) = A\,\bigl|1 - B\,e^{-\text{TI}/T_1}\bigr|
$$ (eq-ir-3p)

where:

- $A = M_0\sin\alpha\,(1 - e^{-\text{TR}/T_1})$ absorbs proton density
  and coil factors.
- $B = 1 + e^{-\text{TR}/T_1}$ accounts for incomplete recovery between
  inversions (for TR $\gg T_1$, $B \to 2$; for shorter TR, $B < 2$).

Fitting [Eq. %s](eq-ir-3p) for the three unknowns $A$, $B$, $T_1$
yields a T1 estimate that is robust to incomplete TR recovery.

### Null Point

The magnetisation crosses zero at the **null time**:

$$
\text{TI}_\text{null} = T_1 \ln(B)
$$ (eq-ir-null)

At this TI, tissues with the corresponding T1 appear dark. IR is widely
used for tissue nulling (e.g., FLAIR nulls CSF with TI $\approx 2200$ ms
at 3 T; STIR nulls fat with TI $\approx 170$ ms at 1.5 T).

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

TI_arr = np.linspace(50, 4000, 600)   # ms

def ir_signal(TI, T1, TR=10000, B1_eff=1.0, A=1.0):
    """Three-parameter IR magnitude signal."""
    B = 1 + np.exp(-TR / T1)
    B_eff = B1_eff * 2   # inversion efficiency scales with B1
    # clamp B_eff to [0, 2]
    B_eff = min(B_eff, 2.0)
    return A * np.abs(1 - B_eff * np.exp(-TI / T1))

fig = make_subplots(rows=1, cols=2,
                    subplot_titles=['IR Curves for Different Tissue Types',
                                    'Effect of Imperfect Inversion (B1 Error)'])

tissues = [('WM  T1=840 ms',  840,  '#1f77b4'),
           ('GM  T1=1300 ms', 1300, '#d62728'),
           ('CSF T1=4500 ms', 4500, '#2ca02c')]

for name, T1, col in tissues:
    S = [ir_signal(ti, T1) for ti in TI_arr]
    fig.add_trace(go.Scatter(x=TI_arr, y=S, name=name,
                             line=dict(color=col, width=2.2)), row=1, col=1)
    # Null point
    TI_null = T1 * np.log(2)
    fig.add_trace(go.Scatter(x=[TI_null], y=[0],
                             mode='markers',
                             marker=dict(symbol='x', color=col, size=9),
                             showlegend=False), row=1, col=1)

# Right panel: B1 effect on IR of GM
T1_gm = 1300
B1_vals = [0.7, 0.85, 1.0, 1.15]
cols_b = ['#d6604d', '#f4a582', '#4393c3', '#2166ac']
for b1, cv in zip(B1_vals, cols_b):
    S_b1 = [ir_signal(ti, T1_gm, B1_eff=b1) for ti in TI_arr]
    fig.add_trace(go.Scatter(x=TI_arr, y=S_b1,
                             name=f'B1 = {b1:.2f}',
                             line=dict(color=cv, width=2)), row=1, col=2)

fig.update_xaxes(title_text='TI (ms)', row=1, col=1)
fig.update_xaxes(title_text='TI (ms)', row=1, col=2)
fig.update_yaxes(title_text='Signal |Mz| (a.u.)', row=1, col=1)
fig.update_yaxes(title_text='Signal |Mz| (a.u.)', row=1, col=2)
fig.update_layout(
    height=430,
    title_text='IR Simulation  (TR = 10 000 ms; ×  marks null-crossing)',
    legend=dict(orientation='h', yanchor='bottom', y=-0.28,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows the magnitude IR curves for three brain tissues.
Crosses mark the tissue-specific null times (Eq. [%s](eq-ir-null)).
WM has the shortest T1 and reaches full recovery soonest. The right
panel illustrates the effect of B1 inhomogeneity: an imperfect
inversion pulse (B1 ≠ 1) shifts the null time and distorts the
recovery curve, introducing a systematic T1 bias if not corrected.

## Data Fitting

T1 is extracted by fitting the three-parameter model
(Eq. [%s](eq-ir-3p)) to the magnitude IR data using a non-linear
least-squares algorithm.

```{code-cell} python
:tags: [hide-input]

import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objects as go

rng = np.random.default_rng(5)

# Simulate 8-point IR acquisition for a WM voxel
TI_pts = np.array([50, 100, 200, 400, 700, 1200, 2000, 3500])  # ms
T1_true = 840.0; A_true = 1.0; B_true = 2.0; TR = 10000.0
noise = 0.025

S_true  = A_true * np.abs(1 - B_true * np.exp(-TI_pts / T1_true))
S_noisy = S_true + rng.normal(0, noise, len(TI_pts))
S_noisy = np.maximum(S_noisy, 0)   # magnitude: non-negative

def ir_model(TI, A, B, T1):
    return A * np.abs(1 - B * np.exp(-TI / T1))

popt, pcov = curve_fit(ir_model, TI_pts, S_noisy, p0=[1.0, 2.0, 800.0],
                       bounds=([0, 0.5, 100], [5, 2.5, 8000]),
                       maxfev=5000)
A_fit, B_fit, T1_fit = popt
T1_err = np.sqrt(np.diag(pcov))[2]

TI_dense = np.linspace(10, 4000, 500)
S_fit_dense = ir_model(TI_dense, *popt)
S_true_dense = A_true * np.abs(1 - B_true * np.exp(-TI_dense / T1_true))

fig = go.Figure()
fig.add_trace(go.Scatter(x=TI_dense, y=S_true_dense,
                         name='Ground truth', mode='lines',
                         line=dict(color='gray', dash='dot', width=1.5)))
fig.add_trace(go.Scatter(x=TI_pts, y=S_noisy,
                         name='Noisy data', mode='markers',
                         marker=dict(color='royalblue', size=9)))
fig.add_trace(go.Scatter(x=TI_dense, y=S_fit_dense,
                         name=f'Fit  T1 = {T1_fit:.0f} ± {T1_err:.0f} ms',
                         line=dict(color='tomato', width=2.2)))

fig.update_layout(
    title=(f'IR T1 Fitting  (True T1 = {T1_true:.0f} ms, '
           f'Fitted T1 = {T1_fit:.0f} ms)'),
    xaxis_title='TI (ms)',
    yaxis_title='Signal (a.u.)',
    template='plotly_white',
    height=380,
    legend=dict(x=0.65, y=0.05),
)
fig.show()
```

The non-linear fit recovers T1 accurately from sparse, noisy data.
Fitting all three parameters $(A, B, T_1)$ is more robust than
assuming $B = 2$ (perfect inversion), particularly when the TR is
short or the B1 field is inhomogeneous. The confidence interval on
T1 (from the covariance matrix) quantifies how SNR and TI spacing
affect precision.

## Example Brain Maps

The following simulates a four-TI inversion recovery acquisition on a
brain phantom and recovers a T1 map using a template-matching look-up table.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

rng = np.random.default_rng(3)

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

TI_pts = np.array([400., 800., 1600., 3200.])
noise  = 0.025

def rician(img, s):
    re = img + rng.normal(0, s, img.shape)
    im = rng.normal(0, s, img.shape)
    return np.sqrt(re**2 + im**2)

imgs = [rician(PDm * np.abs(1 - 2*np.exp(-ti / np.where(T1m > 0, T1m, 1))) * bmask, noise)
        for ti in TI_pts]

# Fast T1 estimation: template-matching LUT
T1_lut  = np.logspace(np.log10(100), np.log10(6000), 400)
tmpls   = np.abs(1 - 2*np.exp(-TI_pts[:, None] / T1_lut[None, :]))
tmpls_n = tmpls / np.maximum(tmpls.max(0), 1e-9)
S_flat  = np.stack(imgs).reshape(len(TI_pts), -1)
S_n     = S_flat / np.maximum(S_flat.max(0), 1e-9)
dp      = tmpls_n.T @ S_n     # (400, n²)
T1_fit  = T1_lut[dp.argmax(0)].reshape(n, n)
T1_fit  = np.where(bmask, T1_fit, np.nan)

fig, axes = plt.subplots(1, 5, figsize=(17, 3.8))
for ax, img, ti in zip(axes[:4], imgs, TI_pts):
    ax.imshow(img, cmap='gray', vmin=0, vmax=0.9, interpolation='bilinear')
    ax.set_title(f'TI = {int(ti)} ms', fontsize=9); ax.axis('off')
im = axes[4].imshow(T1_fit, cmap='magma', vmin=400, vmax=4500,
                    interpolation='bilinear')
axes[4].set_title('IR T1 map', fontsize=9); axes[4].axis('off')
div = make_axes_locatable(axes[4])
cax = div.append_axes('right', size='5%', pad=0.04)
plt.colorbar(im, cax=cax, label='T1 (ms)')

fig.suptitle('Inversion Recovery Brain Maps  (TR = 10 000 ms)', fontsize=11, y=1.02)
plt.tight_layout()
plt.show()
```

At TI = 400 ms the WM (T1 ≈ 840 ms) and GM (T1 ≈ 1300 ms) have both
partially recovered from the inversion and appear bright, while CSF
(T1 ≈ 4500 ms) is still largely negative and appears dark. At
TI = 800 ms WM is near its null crossing and appears darker than GM.
By TI = 3200 ms all parenchyma has recovered and only the ventricles
remain dark. The T1 map correctly assigns long T1 (bright, warm
colours) to the ventricles/CSF and progressively shorter T1 to GM and
WM.

## References

```{bibliography}
:filter: docname in docnames
```
