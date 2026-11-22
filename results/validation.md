# Validation

## 1. SWS example run (100 kg, A*Cd = 1 m^2, 300 -> 180 km, F10.7 = 70, Ap = 0)

| | Lifetime [days] |
|---|---|
| SWS note, section 4 | 46 |
| This code (RK4, dt = 60 s) | 46.9 |
| Difference | +2.0 % |

## 2. Time-step convergence (time-varying area, F10.7 = 150, start 500 km)

| dt [s] | Steps per orbit | Lifetime [days] | Change vs. dt = 15 s |
|---|---|---|---|
| 15 | 377 | 561.48 | +0.000 % |
| 60 | 94 | 561.48 | -0.000 % |
| 300 | 19 | 561.48 | +0.000 % |
| 900 | 6 | 561.48 | +0.000 % |
| 1800 | 3 | 561.60 | +0.022 % |
| 2700 | 2 | 561.25 | -0.041 % |
