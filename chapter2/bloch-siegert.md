---
title: "2.3 Bloch-Siegert Shift"
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

# 2.3 Bloch-Siegert Shift

The Bloch-Siegert shift B1 mapping method, introduced by Sacolick et al. in
2010 {cite}`sacolick2010`, exploits a phase shift induced by an off-resonance
RF pulse to map B1. Unlike amplitude-based methods (DAM, AFI), the
Bloch-Siegert approach encodes B1 into the image *phase* rather than
magnitude, making it insensitive to T1, flip angle errors, and proton density.

The underlying physical phenomenon is the Bloch-Siegert shift: an
off-resonance RF field shifts the effective resonance frequency of nuclear
spins, causing them to accumulate an additional phase proportional to the
square of the B1 amplitude {cite}`bloch1940`.

## Pulse Sequence

A Bloch-Siegert B1 mapping sequence adds an off-resonance Fermi pulse (or
similar shaped RF pulse) to a standard imaging sequence (e.g., a GRE or SE
readout):

1. **Excitation** (standard 90° or $\alpha$ pulse).
2. **Off-resonance BS pulse** — a shaped RF pulse applied at offset $\pm\Delta\omega$
   from resonance. It is far enough off-resonance to not excite the spins, but
   close enough to induce a measurable phase shift.
3. **Readout** — standard gradient echo or spin echo.
4. **Second acquisition** with the BS pulse frequency at $-\Delta\omega$ (or
   at $+\Delta\omega$ in the alternating-sign version).

The B1 map is derived from the **phase difference** between the two
acquisitions, which doubles the BS phase shift and cancels other phase
contributions (B0 offsets, eddy currents).

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(4, 1, figsize=(11, 5.5), sharex=True)
fig.subplots_adjust(hspace=0.05)

t = np.linspace(0, 1, 2000)
TR = 1.0

def rect(t, t0, w, h=1.0):
    r = np.zeros_like(t)
    r[(t >= t0) & (t <= t0 + w)] = h
    return r

def sinc_pulse(t, t0, w, h=1.0):
    x = (t - (t0 + w/2)) / (w/6)
    v = h * np.sinc(x)
    v[(t < t0) | (t > t0 + w)] = 0
    return v

def fermi_pulse(t, t0, w, h=1.0, a=15.0):
    """Fermi pulse envelope."""
    tc = t0 + w/2
    v = h / (1 + np.exp(a * (np.abs(t - tc) - w/2.5)))
    v[(t < t0) | (t > t0 + w)] = 0
    return v

# RF: excitation + BS Fermi pulse
rf_exc = sinc_pulse(t, 0.02, 0.06, 0.60)
# BS pulse (off-resonance — drawn as envelope, different colour)
rf_bs = fermi_pulse(t, 0.14, 0.20, 0.90)

axes[0].fill_between(t, 0, rf_exc, color='royalblue', alpha=0.8, label='Excitation')
axes[0].fill_between(t, 0, rf_bs,  color='coral',     alpha=0.7, label='BS Fermi pulse (off-res.)')
axes[0].axhline(0, color='k', linewidth=0.8)
axes[0].annotate('90°', (0.05, 0.64), ha='center', fontsize=9)
axes[0].annotate('+Δω', (0.24, 0.95), ha='center', fontsize=9, color='coral')
axes[0].legend(loc='upper right', fontsize=8)
axes[0].set_ylabel('RF', fontsize=9)
axes[0].set_ylim(-0.15, 1.15)
axes[0].yaxis.set_tick_params(labelleft=False)

# Gz
gz = rect(t, 0.01, 0.08, 0.70) + rect(t, 0.09, 0.04, -0.35)
gz += rect(t, 0.35, 0.04, 0.50) + rect(t, 0.39, 0.03, -0.30)  # readout slice
axes[1].fill_between(t, 0, gz, where=gz > 0, color='forestgreen', alpha=0.7)
axes[1].fill_between(t, 0, gz, where=gz < 0, color='tomato', alpha=0.7)
axes[1].axhline(0, color='k', linewidth=0.8)
axes[1].set_ylabel('Gz', fontsize=9)
axes[1].yaxis.set_tick_params(labelleft=False)

# Gy
gy = rect(t, 0.37, 0.05, 0.55)
axes[2].fill_between(t, 0, gy, color='darkorange', alpha=0.7)
axes[2].axhline(0, color='k', linewidth=0.8)
axes[2].set_ylabel('Gy', fontsize=9)
axes[2].yaxis.set_tick_params(labelleft=False)

# Gx + ADC
TE_loc = 0.48
gx = rect(t, 0.37, 0.05, -0.45) + rect(t, TE_loc - 0.07, 0.14, 0.45)
axes[3].fill_between(t, 0, gx, where=gx > 0, color='purple', alpha=0.7)
axes[3].fill_between(t, 0, gx, where=gx < 0, color='plum', alpha=0.7)
axes[3].fill_between(t, -0.15, -0.05,
                     where=((t >= TE_loc - 0.07) & (t <= TE_loc + 0.07)),
                     color='gold', alpha=0.9, label='ADC')
axes[3].axhline(0, color='k', linewidth=0.8)
axes[3].set_ylabel('Gx / ADC', fontsize=9)
axes[3].yaxis.set_tick_params(labelleft=False)
axes[3].set_xlim(0, TR)
axes[3].set_xlabel('Time (one TR) — acquire twice with ±Δω BS pulse', fontsize=9)
axes[3].set_ylim(-0.30, 0.65)

# Bracket for BS pulse
axes[0].annotate('', xy=(0.34, -0.12), xytext=(0.14, -0.12),
                 arrowprops=dict(arrowstyle='<->', color='coral', lw=1.5))
axes[0].text(0.24, -0.14, 'BS pulse duration τ_p', ha='center',
             fontsize=8, color='coral')

axes[0].set_title('Bloch-Siegert Shift B1 Mapping Sequence', fontsize=11, pad=8)
plt.tight_layout()
plt.close()
```

:::{figure} bs_pulse_sequence.png
:name: fig-bs-sequence
:align: center
Bloch-Siegert B1 mapping pulse sequence. After the 90° excitation, an
off-resonance Fermi pulse (coral) at offset $+\Delta\omega$ induces a
phase shift proportional to B1². A second acquisition (not shown) applies
the pulse at $-\Delta\omega$; the phase difference yields the B1 map.
:::

## Mathematical Model

### Bloch-Siegert Phase Shift

For a shaped RF pulse with complex envelope $B_1^+(t) = B_1^0\,f(t)$
applied at off-resonance frequency $\Delta\omega \gg \gamma B_1^0$, the
accumulated phase shift in the rotating frame is {cite}`sacolick2010`:

$$
\phi_\text{BS} = \int_0^{\tau_p} \frac{\gamma^2 |B_1^+(t)|^2}{2\Delta\omega}\,dt
  = \frac{\gamma^2 (B_1^0)^2}{2\Delta\omega} \int_0^{\tau_p} |f(t)|^2\,dt
  = K_\text{BS} \cdot (B_1^0)^2
$$ (eq-bs-phase)

where $K_\text{BS}$ is a pulse-shape-dependent constant:

$$
K_\text{BS} = \frac{\gamma^2}{2\Delta\omega}\int_0^{\tau_p} |f(t)|^2\,dt
$$

### Differential Phase Measurement

Two acquisitions at $+\Delta\omega$ and $-\Delta\omega$ accumulate phase
shifts $+\phi_\text{BS}$ and $-\phi_\text{BS}$ respectively (plus common
background phase $\phi_0$ from B0, eddy currents, etc.). The differential
phase is:

$$
\Delta\phi = \phi_+  - \phi_- = 2\,K_\text{BS}\,(B_1^0)^2
$$ (eq-bs-delta-phi)

Solving for the B1 amplitude:

$$
B_1^0 = \sqrt{\frac{\Delta\phi}{2\,K_\text{BS}}}
$$ (eq-bs-b1)

The actual flip angle and B1 scaling factor follow directly from $B_1^0$ and
the nominal pulse calibration.

### $K_\text{BS}$ for a Fermi Pulse

For a Fermi pulse with peak amplitude $B_1^\text{max}$, duration $\tau_p$,
and transition width parameter $a$, $K_\text{BS}$ can be computed numerically:

$$
K_\text{BS} = \frac{\gamma^2}{2\Delta\omega}
  \int_0^{\tau_p} \left[\frac{1}{1+e^{a(|t - \tau_p/2| - \tau_p/2.5)}}\right]^2 dt
$$

In practice $K_\text{BS}$ is calibrated from the known pulse shape rather than
measured in the scanner.

## Simulations

```{code-cell} python
:tags: [hide-input]

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

gamma = 2 * np.pi * 42.577e6   # rad/(s·T)

# Fermi pulse shape
def fermi_shape(t_arr, tau_p, a=15.0):
    tc = tau_p / 2
    f = 1.0 / (1 + np.exp(a * (np.abs(t_arr - tc) - tau_p / 2.5)))
    return f

tau_p = 8e-3   # 8 ms pulse
dt = 1e-5
t_arr = np.arange(0, tau_p, dt)
f_t = fermi_shape(t_arr, tau_p)
integral_f2 = np.trapz(f_t**2, t_arr)

# K_BS vs off-resonance frequency
dw_arr = np.linspace(500, 4000, 300) * 2 * np.pi   # rad/s
K_BS_arr = gamma**2 * integral_f2 / (2 * dw_arr)

# Phase shift vs B1 for a fixed off-resonance
dw_fixed = 2000 * 2 * np.pi   # 2 kHz
K_BS_fixed = gamma**2 * integral_f2 / (2 * dw_fixed)

B1_arr = np.linspace(0, 25e-6, 300)   # T
phi_BS = K_BS_fixed * B1_arr**2       # radians

fig = make_subplots(rows=1, cols=2,
                    subplot_titles=[
                        'Fermi Pulse Shape and K_BS vs Δω/2π',
                        'Bloch-Siegert Phase vs B1 Amplitude'])

# Subplot 1a: Fermi pulse
t_ms = t_arr * 1e3
fig.add_trace(go.Scatter(x=t_ms, y=f_t,
                         name='Fermi envelope f(t)',
                         line=dict(color='coral', width=2.5)), row=1, col=1)

# Subplot 1b: K_BS vs Δω (secondary y — we overlay by scaling)
dw_Hz = dw_arr / (2 * np.pi)
K_scale = K_BS_arr * 1e9   # scale for visibility
fig.add_trace(go.Scatter(x=dw_Hz, y=K_scale,
                         name='K_BS × 10⁹ (rad/T², right axis)',
                         line=dict(color='royalblue', width=2),
                         xaxis='x2', yaxis='y3'), row=1, col=1)

# Subplot 2: phase vs B1
B1_uT = B1_arr * 1e6
fig.add_trace(go.Scatter(x=B1_uT, y=phi_BS,
                         name=f'Δω/2π = 2 kHz, τ_p = {tau_p*1e3:.0f} ms',
                         line=dict(color='royalblue', width=2.5)), row=1, col=2)
# 2π reference
fig.add_hline(y=2*np.pi, line_dash='dash', line_color='gray',
              annotation_text='2π rad', row=1, col=2)

# Noise floor (1° = π/180 rad, SNR-limited)
sigma_phi = np.radians(1.0)   # 1° phase noise
B1_noise = np.sqrt(sigma_phi / (2 * K_BS_fixed))
fig.add_vline(x=B1_noise * 1e6, line_dash='dot', line_color='tomato',
              annotation_text=f'Detection limit\n~{B1_noise*1e6:.1f} μT',
              row=1, col=2)

fig.update_xaxes(title_text='Time (ms)', row=1, col=1)
fig.update_xaxes(title_text='B1⁺ Amplitude (μT)', row=1, col=2)
fig.update_yaxes(title_text='Fermi envelope (a.u.)', row=1, col=1)
fig.update_yaxes(title_text='BS Phase Shift φ_BS (rad)', row=1, col=2)

fig.update_layout(
    height=420,
    title_text='Bloch-Siegert Shift Simulation',
    legend=dict(orientation='h', yanchor='bottom', y=-0.3,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The Fermi pulse shape (left panel, coral curve) determines $K_\text{BS}$ via
the integral of $|f(t)|^2$. The right panel shows the quadratic dependence of
the BS phase on B1 amplitude ([Eq. %s](eq-bs-phase)). The tomato vertical
line marks the detection limit for 1° phase noise; below this B1, the method
becomes unreliable at this off-resonance frequency.

## Data Fitting

```{code-cell} python
:tags: [hide-input]

import numpy as np
from scipy.ndimage import gaussian_filter
import plotly.graph_objects as go
from plotly.subplots import make_subplots

rng = np.random.default_rng(55)

# Simulation parameters
gamma = 2 * np.pi * 42.577e6
tau_p = 8e-3
dt = 1e-5
t_arr = np.arange(0, tau_p, dt)
a_fermi = 15.0
f_t = 1.0 / (1 + np.exp(a_fermi * (np.abs(t_arr - tau_p/2) - tau_p/2.5)))
integral_f2 = np.trapz(f_t**2, t_arr)
dw = 2000 * 2 * np.pi  # rad/s
K_BS = gamma**2 * integral_f2 / (2 * dw)

N = 150  # 1D profile pixels
x = np.linspace(0, np.pi, N)
B1_true = (10 + 6 * np.cos(x)**2 + 2 * np.sin(2*x)) * 1e-6  # T

# Phase images with noise (σ = 2°)
sigma_phi = np.radians(2.0)
phi_pos = K_BS * B1_true**2 + rng.normal(0, sigma_phi, N)
phi_neg = -K_BS * B1_true**2 + rng.normal(0, sigma_phi, N)

# B1 map reconstruction
delta_phi = phi_pos - phi_neg  # = 2 * K_BS * B1^2
B1_meas_raw = np.sqrt(np.abs(delta_phi) / (2 * K_BS))
B1_meas_smooth = gaussian_filter(B1_meas_raw, sigma=3)   # spatial smoothing

pixel = np.arange(N)
B1_uT = B1_true * 1e6
B1_raw_uT = B1_meas_raw * 1e6
B1_sm_uT = B1_meas_smooth * 1e6

rmse_raw = np.sqrt(np.mean((B1_raw_uT - B1_uT)**2))
rmse_sm = np.sqrt(np.mean((B1_sm_uT - B1_uT)**2))

fig = make_subplots(rows=1, cols=2,
                    subplot_titles=['Phase Images (±Δω)',
                                    'B1 Map Reconstruction'])

fig.add_trace(go.Scatter(x=pixel, y=np.degrees(phi_pos), name='φ_+ (noisy)',
                         line=dict(color='royalblue', width=1.5)), row=1, col=1)
fig.add_trace(go.Scatter(x=pixel, y=np.degrees(phi_neg), name='φ_− (noisy)',
                         line=dict(color='tomato', width=1.5)), row=1, col=1)
fig.add_trace(go.Scatter(x=pixel, y=np.degrees(K_BS * B1_true**2),
                         name='True φ_+', line=dict(color='royalblue',
                         dash='dot', width=1.2)), row=1, col=1)

fig.add_trace(go.Scatter(x=pixel, y=B1_uT, name='True B1 (μT)',
                         line=dict(color='gray', dash='dot', width=1.5)),
             row=1, col=2)
fig.add_trace(go.Scatter(x=pixel, y=B1_raw_uT,
                         name=f'Raw (RMSE={rmse_raw:.2f} μT)',
                         line=dict(color='coral', width=1.5, dash='dash')),
             row=1, col=2)
fig.add_trace(go.Scatter(x=pixel, y=B1_sm_uT,
                         name=f'Smoothed (RMSE={rmse_sm:.2f} μT)',
                         line=dict(color='royalblue', width=2)),
             row=1, col=2)

fig.update_xaxes(title_text='Pixel', row=1, col=1)
fig.update_xaxes(title_text='Pixel', row=1, col=2)
fig.update_yaxes(title_text='Phase (degrees)', row=1, col=1)
fig.update_yaxes(title_text='B1⁺ Amplitude (μT)', row=1, col=2)

fig.update_layout(
    height=400,
    title_text='Bloch-Siegert B1 Map from Phase Difference (σ_φ = 2°)',
    legend=dict(orientation='h', yanchor='bottom', y=-0.3,
                xanchor='center', x=0.5),
    template='plotly_white',
)
fig.show()
```

The left panel shows the noisy phase images at $\pm\Delta\omega$; their
difference (proportional to $B_1^{+2}$) is used in [Eq. %s](eq-bs-b1).
The right panel compares the raw and Gaussian-smoothed B1 maps against
the ground truth: smoothing exploits the known spatial smoothness of B1
to suppress noise without introducing significant bias.

## Example Brain Maps

A simulated Bloch-Siegert acquisition shows the anatomy-weighted
magnitude image, the BS phase difference map, and the recovered B1 map.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

rng = np.random.default_rng(2)

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
T2v={0:0,1:0,3:100, 4:75,  5:85,  6:80,  7:1500}
PDv={0:0,1:0,3:0.85,4:0.75,5:0.80,6:0.78,7:1.00}
tissue, bmask, B1m = brain_phantom()
T1m = np.vectorize(T1v.get)(tissue).astype(float)
T2m = np.vectorize(T2v.get)(tissue).astype(float)
PDm = np.vectorize(PDv.get)(tissue).astype(float)

KBS = 18.0   # normalised KBS constant
phi_plus  = KBS * B1m**2 + rng.normal(0, 0.04, B1m.shape)
phi_minus = -KBS * B1m**2 + rng.normal(0, 0.04, B1m.shape)
phi_diff  = (phi_plus - phi_minus) * bmask
B1_bs = np.where(bmask, np.sqrt(np.maximum(phi_diff, 0) / (2 * KBS)), np.nan)

img_mag = PDm * np.exp(-80. / np.where(T2m > 0, T2m, 1)) * bmask
re = img_mag + rng.normal(0, 0.02, img_mag.shape)
im_noise = rng.normal(0, 0.02, img_mag.shape)
img_mag_n = np.sqrt(re**2 + im_noise**2)

fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
axes[0].imshow(img_mag_n, cmap='gray', vmin=0, vmax=0.5, interpolation='bilinear')
axes[0].set_title('Magnitude image', fontsize=10); axes[0].axis('off')

pd_show = np.where(bmask, phi_diff, np.nan)
im2 = axes[1].imshow(pd_show, cmap='bwr', vmin=-2*KBS, vmax=2*KBS,
                     interpolation='bilinear')
axes[1].set_title('BS phase difference (rad)', fontsize=10); axes[1].axis('off')
div1 = make_axes_locatable(axes[1])
cax1 = div1.append_axes('right', size='5%', pad=0.04)
plt.colorbar(im2, cax=cax1, label='rad')

im3 = axes[2].imshow(B1_bs, cmap='RdBu_r', vmin=0.7, vmax=1.3,
                     interpolation='bilinear')
axes[2].set_title('BS B1 map', fontsize=10); axes[2].axis('off')
div2 = make_axes_locatable(axes[2])
cax2 = div2.append_axes('right', size='5%', pad=0.04)
plt.colorbar(im3, cax=cax2, label='B1 factor')

fig.suptitle(f'Bloch-Siegert Brain Maps', fontsize=11, y=1.02)
plt.tight_layout()
plt.show()
```

The phase difference map (centre panel) directly encodes B1²: the
centre of the brain, where B1 is highest (~1.22), shows the largest
phase accumulation. The B1 map (right) recovers the spatial field
pattern via $\sqrt{\Delta\Phi / (2K_\text{BS})}$, with excellent
agreement to the known phantom input.

## References

```{bibliography}
:filter: docname in docnames
```
