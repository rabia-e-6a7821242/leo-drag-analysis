"""
compare_msis.py
Compares the simple SWS density model (Fortran) with the empirical
NRLMSIS 2.1 model (Emmert et al., 2021), using the `pymsis` package.

NRLMSIS depends on latitude, local time and season, while the SWS model gives
one value per altitude. To compare them fairly, the MSIS density is averaged
over a latitude/longitude grid (area-weighted with cos(latitude)) on two dates
(March equinox and June solstice, 12 UT). The same F10.7 (daily and 81-day
mean) and Ap values as in the Fortran profiles are given to MSIS by hand, so
no space-weather files have to be downloaded.

Run from the repository root (after density_profile has been run):
    python3 python/compare_msis.py
"""

import sys
from pathlib import Path

import numpy as np

try:
    import pymsis
except ImportError:
    sys.exit("pymsis is not installed (pip install pymsis), skipping MSIS comparison")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_results import C, GREY, INK2, load, save  # same style as other figures

ROOT = Path(__file__).resolve().parent.parent

DATES = [np.datetime64("2019-03-20T12:00"), np.datetime64("2019-06-21T12:00")]
LATS = np.arange(-87.5, 90, 5.0)
LONS = np.arange(0, 360, 10.0)
CASES = [(70.0, "f107_70", "F10.7 = 70"),
         (150.0, "f107_150", "F10.7 = 150"),
         (230.0, "f107_230", "F10.7 = 230")]
AP = 15.0


def msis_global_mean(alts, f107, ap):
    """Area-weighted global mean density [kg/m^3], averaged over DATES."""
    w = np.cos(np.radians(LATS))[None, :, None]
    means = []
    for date in DATES:
        out = pymsis.calculate(date, LONS, LATS, alts, f107s=[f107], f107as=[f107],
                               aps=[[ap] * 7], version=2.1)
        rho = np.squeeze(out[..., pymsis.Variable.MASS_DENSITY])   # (lon, lat, alt)
        means.append((rho * w).sum(axis=(0, 1)) / (w.sum() * len(LONS)))
    return np.mean(means, axis=0)


def main():
    prof = load("density_profiles")
    keep = prof["alt_km"] <= 600
    alts = prof["alt_km"][keep]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    rows = ["| F10.7 | Altitude [km] | SWS [kg/m^3] | NRLMSIS 2.1 [kg/m^3] | SWS / MSIS |",
            "|---|---|---|---|---|"]
    for i, (f107, key, lab) in enumerate(CASES):
        msis = msis_global_mean(alts, f107, AP)
        sws = prof[key][keep]
        ax1.plot(sws, alts, color=C[i], label=f"SWS, {lab}")
        ax1.plot(msis, alts, color=C[i], ls="--", label=f"MSIS, {lab}")
        ax2.plot(sws / msis, alts, color=C[i], label=lab)
        for h in (300.0, 400.0, 500.0):
            j = int(np.argmin(abs(alts - h)))
            rows.append(f"| {f107:.0f} | {h:.0f} | {sws[j]:.3e} | {msis[j]:.3e} | "
                        f"{sws[j] / msis[j]:.2f} |")

    ax1.set_xscale("log")
    ax1.set_xlabel("Neutral density [kg m$^{-3}$]")
    ax1.set_ylabel("Altitude [km]")
    ax1.set_title("a) SWS (solid) vs. NRLMSIS 2.1 (dashed)")
    ax1.legend(fontsize=7, loc="upper right", ncol=1)

    ax2.axvline(1, color=GREY, lw=1)
    ax2.axvspan(0.8, 1.25, color=GREY, alpha=0.12, lw=0)
    ax2.axhspan(500, 600, color=GREY, alpha=0.08, lw=0)
    ax2.text(0.98, 550, "outside SWS\nvalid range", transform=ax2.get_yaxis_transform(),
             ha="right", va="center", fontsize=8, color=INK2)
    ax2.axhline(180, color=GREY, lw=0.8, ls=":")
    ax2.text(0.98, 170, "static table below 180 km", transform=ax2.get_yaxis_transform(),
             ha="right", va="center", fontsize=8, color=INK2)
    ax2.set_xlabel("Density ratio SWS / MSIS (grey band: 0.8–1.25)")
    ax2.set_title("b) Ratio SWS / MSIS")
    ax2.legend(fontsize=8, loc="upper left")
    save(fig, "fig7_sws_vs_msis.png")

    (ROOT / "results" / "msis_comparison.md").write_text("\n".join(rows) + "\n")
    print("saved", ROOT / "results" / "msis_comparison.md")


if __name__ == "__main__":
    main()
