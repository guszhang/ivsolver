# Circuit and Poynting vector calculation

This note describes the frequency-domain features added in IVSolver 2.0: the coupled-circuit solution, and the evaluation of E, B and the complex Poynting vector at observation points. They implement the energy-flow method of the IPTVisual paper [1], which builds on the visualisation method of [2]. The original inductance (`-l`) and flux-density (`-b`) modes are unchanged; see [Compatibility](#compatibility).

## 1. Conventions

- Time-harmonic quantities are phasors with **peak** amplitude: $x(t) = \mathrm{Re}\{X e^{j\omega t}\}$, $\omega = 2\pi f$.
- Complex power is $\frac{1}{2}VI^*$, so its real part is the average power in watts.
- SI units throughout: metres, amperes, volts, henries, farads, hertz, teslas.

## 2. Circuit model

Each wire forms one loop. The loop has its self-inductance $L_{kk}$ from the Neumann integral, plus optional series elements:

| directive | element | default |
|---|---|---|
| `R r` | series resistance $R_k$ (Ω) | 0 |
| `C c` | series capacitance $C_k$ (F); `0` means no capacitor | none |
| `L l` | external series inductance $L_{\mathrm{ext},k}$ (H) | 0 |
| `T f` | choose $C_k = 1/((2\pi f)^2 (L_{kk}+L_{\mathrm{ext},k}))$; `T 0` uses the analysis frequency | off |
| `V re im` | voltage source (V) | 0 V |
| `i re im` | current source (A) | — |

A loop with neither `i` nor `V` is a passive closed loop (0 V source). The loops satisfy

$$
V_k = Z_k I_k + \sum_{m\neq k} j\omega M_{km} I_m,\qquad
Z_k = R_k + j\omega\,(L_{kk}+L_{\mathrm{ext},k}) + \frac{1}{j\omega C_k}.
$$

Rows of voltage-driven loops are solved for their currents by Gaussian elimination with partial pivoting. The source voltages of current-driven loops are then evaluated. The solution is written to `<output>.circuit`:

```text
# wire  driven_by  I_re I_im  V_re V_im  Z_re Z_im  C_series  P_source Q_source  P_resistive
```

where `Z` is the loop's own impedance $Z_k$, `C_series` is the capacitance actually used (after tuning), $P+jQ = \frac{1}{2}V_kI_k^*$, and $P_\mathrm{resistive}=\frac{1}{2}|I_k|^2R_k$. A passive loop is listed as `voltage` driven with a 0 V source.

The circuit is solved whenever `--freq` is given. In mode `-l` it is written alongside the `.lmat` file. In mode `-p` it supplies the currents for the field.

## 3. Fields

For wire segments $\mathrm{d}\mathbf{l}$ carrying current $I$, at an observation point separated from the segment by $\mathbf{R}$:

$$
\mathbf{B} = \frac{\mu_0}{4\pi}\sum I\,\frac{\mathrm{d}\mathbf{l}\times\mathbf{R}}{|\mathbf{R}|^3},\qquad
\mathbf{A} = \frac{\mu_0}{4\pi}\sum \frac{I\,\mathrm{d}\mathbf{l}}{|\mathbf{R}|},\qquad
\mathbf{E} = -j\omega\mathbf{A},\qquad \mathbf{H} = \mathbf{B}/\mu_0 .
$$

The complex Poynting vector and its time-domain form are

$$
\mathbf{S} = \tfrac{1}{2}\,\mathbf{E}\times\mathbf{H}^*,\qquad
\mathbf{s}(t) = \mathrm{Re}\{\mathbf{E}e^{j\omega t}\}\times\mathrm{Re}\{\mathbf{H}e^{j\omega t}\}
= \mathrm{Re}\{\mathbf{S}\} + \mathrm{Re}\{\tfrac{1}{2}\mathbf{E}\times\mathbf{H}\,e^{j2\omega t}\}.
$$

$\mathrm{Re}\{\mathbf{S}\}$ is the time-averaged power flow density (W/m²). $\mathrm{Im}\{\mathbf{S}\}$ is the reactive power density. $\mathbf{s}(t)$ oscillates at $2\omega$ about the average, and is what the animations in IPTVisual show.

**Scope.** $\mathbf{E}$ is the field induced by the time-varying currents. The electrostatic part $-\nabla\phi$, due to charges on compensation capacitors or the winding's distributed capacitance, is not modelled. This is consistent with the quasi-static, lumped-element circuit model and with [1]. Under this model, Poynting's theorem holds exactly for the real power exchanged between loops. The tests check this numerically (§5).

### Integration kernels

- `--kernel exact` (default for `-p`) integrates each straight segment of length $L$ in closed form. With the point at axial coordinate $s$ from the segment head and perpendicular distance $d$:

  $$
  \mathbf{B} = \frac{\mu_0 I}{4\pi}\,\frac{\hat{\mathbf{u}}\times\mathbf{r}_1}{d^2}\left[\frac{L-s}{R_2}+\frac{s}{R_1}\right],\qquad
  \mathbf{A} = \frac{\mu_0 I}{4\pi}\,\hat{\mathbf{u}}\left[\sinh^{-1}\frac{L-s}{d}+\sinh^{-1}\frac{s}{d}\right]
  $$

  The logarithms are evaluated in a cancellation-free form. Inside the conductor ($d<a$), $d$ is replaced by the wire radius $a$, which gives the field of a uniform current: $\mathbf{B}$ grows linearly with $d$ and $\mathbf{A}$ stays finite. Points on or inside a wire therefore give finite values.
- `--kernel midpoint` treats each segment as a current element at its midpoint. This is the original `-b` kernel, accurate away from the wire.

Each output row also contains the distance from the point to the nearest conductor axis. Use it to drop points inside or very close to windings.

## 4. Command line

```bash
ivsolver wire_def result -p obs_def --freq 85000 [--threads N] [--kernel exact|midpoint] [--circuit] [--json]
```

- The circuit is solved unless every wire has an `i` directive. `--circuit` forces it, for example to obtain the source voltages.
- `obs_def` is a grid (`x_min x_max x_n`, then the same for y and z, optionally preceded by `grid`) or a list (`points k` followed by k lines `x y z`).
- `result.poynting` has two `#` header lines, then one row of 22 columns per point:
  `x y z`, `Ex Ey Ez`, `Bx By Bz` and `Sx Sy Sz` (each as real and imaginary parts), and `min_distance`.
- `--json` writes a single document (to stdout with output `-`). It contains the wires and their trajectories, the inductance matrix, the circuit solution, the currents used, the field arrays (flattened `[xr, xi, yr, yi, zr, zi, …]` per point) and timings. The WebAssembly build returns the same document.

Read the text output in Python:

```python
import numpy as np
d = np.loadtxt("result.poynting")
xyz = d[:, 0:3]
E = d[:, 3:9:2] + 1j * d[:, 4:9:2]
B = d[:, 9:15:2] + 1j * d[:, 10:15:2]
S = d[:, 15:21:2] + 1j * d[:, 16:21:2]
S_avg = S.real                      # W/m^2
```

See [examples/poynting_two_coil](../examples/poynting_two_coil/README.md) for a complete example.

## 5. Validation

`make test` runs these checks against independent results:

| check | result |
|---|---|
| on-axis B of a circular loop vs. $\mu_0 I R^2 / 2(R^2+z^2)^{3/2}$ | rel. error ≈ 1e-5 |
| Faraday's law, $\oint\mathbf{E}\cdot\mathrm{d}\mathbf{l} = -j\omega\Phi$ (E from A, Φ from B) | rel. error ≈ 4e-7 |
| Poynting's theorem: $-\oint\mathrm{Re}\{\mathbf{S}\}\cdot\mathrm{d}\mathbf{A}$ around a resistive receiver $=\frac{1}{2}\vert I_2\vert^2R_2$ | rel. error ≈ 2e-4 |
| tuned receiver current $= -j\omega M I_1/R_2$; circuit power balance | machine precision |
| exact vs. midpoint kernel away from the wire | ≈ 1e-6 |
| results independent of the thread count | identical |
| committed examples (`result.lmat`, `wires.obs`) and the original command forms | identical |

## 6. Performance

The inductance matrix and the field evaluation run on all hardware threads (`--threads` to limit). Work is split into fixed chunks, so results do not depend on the thread count. For the two 11-turn helices of the example (about 7,000 segments each), the inductance matrix took 1.0 s with version 1, against 0.5 s single-threaded and 0.1 s on 12 threads with version 2.

## 7. Magnetic cores

Cores use a volume integral method. A body is voxelised into prisms with uniform magnetisation `M_v = chi_v H`
(`chi = mu_r - 1`). With `H0` the free-space field of a winding at 1 A and `G` the exact H-field tensor of a uniformly
magnetised prism (point dipole beyond 8 half-sizes), the cell magnetisations solve

    (1/chi_v) M_v - sum_u G(r_v - r_u) M_u = H0(r_v)

which is symmetric positive definite for `mu_r >= 1` (conjugate gradients, Jacobi preconditioner, tolerance 1e-9).
One solve per winding; everything else is superposition:

    B(r)   = B_air(r) + mu0 (H_M(r) + M(r) inside a cell)
    dL_ij  = sum_v V_v M_j(v) . B_i(v)      (symmetrised; added to the Neumann matrix)

Checks: one cubic cell reproduces the exact cube result (demagnetising factor 1/3); a small sphere in a large loop
matches `dL = V 3(mu_r-1)/(mu_r+2) B0^2/mu0`; `tests/run_tests.py` covers both. Known limits: staircased curved
surfaces converge slowly at high `mu_r`, `A` (hence `E`) is a point-dipole approximation, no saturation or loss.

## Compatibility

- `ivsolver wire_def out`, `ivsolver wire_def out -l x` and `ivsolver wire_def out -b obs_def` behave as before. They produce the same `.lmat`, `.trj` and `.obs` files with the same layout and values, and `-b` still uses the midpoint kernel.
- Existing input files are read unchanged. New in 2.0: `#` comments, upper-case circuit directives, the `points` observation format, and descriptive errors for malformed input. Previously, malformed input was silently misread.

## References

1. C. Zhang et al., "IPTVisual: Visualisation of the Spatial Energy Flows in Inductive Power Transfer Systems with Arbitrary Winding Shapes", *World Electric Vehicle Journal*, 13(4):63, 2022. [doi:10.3390/wevj13040063](https://doi.org/10.3390/wevj13040063)
2. H. Wang, C. Zhang and S. Y. R. Hui, "Visualization of Energy Flow in Wireless Power Transfer Systems," *2019 IEEE Wireless Power Transfer Conference (WPTC)*, London, UK, 2019, pp. 536–541. [doi:10.1109/WPTC45513.2019.9055544](https://doi.org/10.1109/WPTC45513.2019.9055544). Interactive demonstration: <https://guszhang.github.io/IPTVisual/WPW2019Demo>
