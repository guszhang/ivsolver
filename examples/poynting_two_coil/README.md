# Two-coil energy flow (complex Poynting vector)

This is the two-coil system of the IPTVisual paper and the WPW 2019 demonstration (case A). It has two 11-turn helical coils of radius 0.15 m, wound with 0.75 mm radius wire and placed 0.3 m apart.

- **Transmitter:** driven by a 1 V source, with 1 Ω series resistance and a 1 nF series capacitor.
- **Receiver:** a passive loop with a 10 Ω load and a 1 nF series capacitor.

The analysis frequency is 562.698 kHz (ω = 3.5355×10⁶ rad/s), close to the 566 kHz resonance of 79 µH with 1 nF.

Files:

- `wire_def`: the windings and circuit directives (`V`, `R`, `C`);
- `obs_def`: a 9 × 1 × 9 grid in the xz plane through the coil axis.

Run from this folder:

```bash
../../build/ivsolver wire_def result -p obs_def --freq 562698
```

This writes:

- `result.lmat` with the inductance matrix (L = 79.09 µH, M = 2.58 µH, k = 0.033);
- `result.circuit` with the solved currents and powers. The transmitter current is 116 mA, and 57.98 mW is delivered. Of that, 51.2 mW is dissipated in the receiver load, so the efficiency is 88 %;
- `result.poynting` with E, B and S at the 81 points;
- `result_0.trj` and `result_1.trj` with the winding trajectories.

Plot the time-averaged energy flow in the xz plane:

```python
import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("result.poynting")
x, z = d[:, 0], d[:, 2]
Sx, Sz = d[:, 15], d[:, 19]          # Re{Sx}, Re{Sz}
keep = d[:, 21] > 0.01               # drop points within 1 cm of a conductor
plt.quiver(x[keep], z[keep], Sx[keep], Sz[keep], np.log10(np.hypot(Sx, Sz)[keep]))
plt.gca().set_aspect("equal")
plt.xlabel("x (m)"); plt.ylabel("z (m)"); plt.title("Re{S} (colour: log10 |S| in W/m²)")
plt.show()
```

Energy leaves the transmitter (z ≈ 0) and flows into the receiver (z ≈ 0.3 m). The same project opens in [IPTVisual](https://github.com/guszhang/iptvisual) as the example “Two-coil system (WPW 2019 case A)”, where the instantaneous flow can be animated.
