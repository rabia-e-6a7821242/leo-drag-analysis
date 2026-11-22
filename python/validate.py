"""
validate.py
Two simple checks of the Fortran orbit_decay program:

1. Reproduce the example run printed in the SWS/IPS note
   "Satellite Orbital Decay Calculations" (section 4):
   100 kg, effective area A*Cd = 1 m^2, start 300 km, F10.7 = 70, Ap = 0,
   re-entry at 180 km  ->  "Re-entry after 46 days".
   The SWS program uses a simple Euler step of 0.1 day on the orbital period,
   our program uses RK4 on the semi-major axis, so small differences are expected.

2. Time-step convergence: run the same case with different RK4 steps and
   check that the lifetime does not depend on dt.

Run from the repository root (after `make`):
    python3 python/validate.py
"""

import subprocess
from pathlib import Path

from run_scenarios import EXE, make_case, write_namelist, parse_summary

ROOT = Path(__file__).resolve().parent.parent
VAL_DIR = ROOT / "runs" / "validation"


def run(case):
    VAL_DIR.mkdir(parents=True, exist_ok=True)
    nml = VAL_DIR / f"{case['name']}.nml"
    out = VAL_DIR / f"{case['name']}.csv"
    write_namelist(case, nml)
    res = subprocess.run([str(EXE), str(nml), str(out)],
                         capture_output=True, text=True, check=True)
    return float(parse_summary(res.stdout)["LIFETIME_DAYS"])


def main():
    lines = ["# Validation", ""]

    # ---- 1) SWS example ----
    sws_case = make_case("sws_example", "validation",
                         mass_kg=100.0, cd=1.0, area_body=1.0, area_panel=0.0,
                         area_mode="fixed_max", h0_km=300.0, h_end_km=180.0,
                         f107_const=70.0, ap_quiet=0.0, out_dt_s=3600.0)
    ours = run(sws_case)
    ref = 46.0
    lines += ["## 1. SWS example run (100 kg, A*Cd = 1 m^2, 300 -> 180 km, F10.7 = 70, Ap = 0)",
              "",
              "| | Lifetime [days] |",
              "|---|---|",
              f"| SWS note, section 4 | {ref:.0f} |",
              f"| This code (RK4, dt = 60 s) | {ours:.1f} |",
              f"| Difference | {100 * (ours - ref) / ref:+.1f} % |", ""]
    print(f"SWS example: reference {ref:.0f} d, this code {ours:.1f} d")

    # ---- 2) time-step convergence ----
    lines += ["## 2. Time-step convergence (time-varying area, F10.7 = 150, start 500 km)",
              "",
              "| dt [s] | Steps per orbit | Lifetime [days] | Change vs. dt = 15 s |",
              "|---|---|---|---|"]
    results = {}
    for dt in (15, 60, 300, 900, 1800, 2700):
        case = make_case(f"dt_{dt}", "validation", area_mode="varying",
                         f107_const=150.0, dt_s=float(dt))
        results[dt] = run(case)
    for dt, life in results.items():
        change = 100 * (life - results[15]) / results[15]
        lines.append(f"| {dt} | {5650 / dt:.0f} | {life:.2f} | {change:+.3f} % |")
        print(f"dt = {dt:4d} s  lifetime = {life:.2f} d")

    (ROOT / "results" / "validation.md").write_text("\n".join(lines) + "\n")
    print("saved", ROOT / "results" / "validation.md")


if __name__ == "__main__":
    main()
