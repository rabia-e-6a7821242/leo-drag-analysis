"""
plot_results.py
Makes all figures in results/figures/ and a markdown summary table from
the CSV files written by the Fortran programs.

Run from the repository root (after run_scenarios.py):
    python3 python/plot_results.py
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "results" / "data"
FIGS = ROOT / "results" / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

# colour-blind checked categorical palette (fixed order) + neutral grey
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
GREY = "#8a8985"
INK = "#0b0b0b"
INK2 = "#52514e"

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 150, "font.size": 10,
    "axes.edgecolor": "#c9c8c2", "axes.labelcolor": INK2,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": INK2, "ytick.color": INK2, "lines.linewidth": 2,
    "legend.frameon": False,
})


def load(name):
    """Read a CSV (skipping '#' comment lines) into a dict of numpy arrays."""
    arr = np.genfromtxt(DATA / f"{name}.csv", delimiter=",", names=True,
                        comments="#")
    return {k: arr[k] for k in arr.dtype.names}


def label_end(ax, x, y, text, color=INK2, dy=0, dx=4):
    ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points",
                color=color, fontsize=9, va="center",
                ha="left" if dx >= 0 else "right")


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIGS / name)
    plt.close(fig)
    print("saved", FIGS / name)


# ---------------------------------------------------------------------------
def fig_density_profiles():
    d = load("density_profiles")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)

    series = [("f107_70", "F10.7 = 70 (solar min)", C[0]),
              ("f107_150", "F10.7 = 150", C[1]),
              ("f107_230", "F10.7 = 230 (solar max)", C[2]),
              ("f107_150_storm", "F10.7 = 150, Ap = 200 (storm)", C[3])]
    ax1.plot(d["static"], d["alt_km"], color=GREY, ls="--",
             label="Static exponential (Vallado)")
    for key, lab, col in series:
        ax1.plot(d[key], d["alt_km"], color=col, label=lab)
    ax1.set_xscale("log")
    ax1.set_xlabel("Neutral density [kg m$^{-3}$]")
    ax1.set_ylabel("Altitude [km]")
    ax1.set_title("a) Thermospheric density profiles")
    ax1.legend(fontsize=8, loc="upper right")

    for key, lab, col in series[1:]:
        ax2.plot(d[key] / d["f107_70"], d["alt_km"], color=col, label=lab)
    ax2.axvline(1, color=GREY, lw=1)
    ax2.set_xlabel("Density ratio to solar minimum (F10.7 = 70)")
    ax2.set_title("b) Solar/geomagnetic amplification")
    ax2.legend(fontsize=8, loc="lower right")
    save(fig, "fig1_density_profiles.png")


def fig_solar_decay():
    fig, ax = plt.subplots(figsize=(8, 4.8))
    cases = [("solar_low", "Low (F10.7 = 70)", C[0]),
             ("solar_moderate", "Moderate (F10.7 = 150)", C[1]),
             ("solar_high", "High (F10.7 = 230)", C[2])]
    s = load("solar_static_model")
    ax.plot(s["time_day"] / 365.25, s["alt_km"], color=GREY, ls="--",
            label="Static model (no solar dependence)")
    for name, lab, col in cases:
        d = load(name)
        ax.plot(d["time_day"] / 365.25, d["alt_km"], color=col, label=lab)
        label_end(ax, d["time_day"][-1] / 365.25, 160,
                  f"{d['time_day'][-1]:.0f} d", dy=6,
                  dx=-4 if name == "solar_moderate" else 4)
    ax.set_xlabel("Time since start [years]")
    ax.set_ylabel("Altitude [km]")
    ax.set_title("Altitude decay from 500 km under different solar activity")
    ax.set_ylim(140, 510)
    ax.legend(loc="lower center", bbox_to_anchor=(0.64, 0.04))
    save(fig, "fig2_decay_solar_activity.png")


def fig_solar_cycle():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True,
                                   gridspec_kw={"height_ratios": [1, 2]})
    for name, lab, col in [("cycle_launch_at_min", "Launch at solar minimum", C[0]),
                           ("cycle_launch_at_max", "Launch at solar maximum", C[1])]:
        d = load(name)
        t = d["time_day"] / 365.25
        ax1.plot(t, d["f107"], color=col, label=lab)
        ax2.plot(t, d["alt_km"], color=col, label=lab)
        label_end(ax2, t[-1], 160, f"{d['time_day'][-1]:.0f} d", dy=6)
    ax1.set_ylabel("F10.7 [sfu]")
    ax1.set_title("a) Idealised 11-year solar cycle")
    ax1.legend(loc="lower right", fontsize=8)
    ax2.set_ylabel("Altitude [km]")
    ax2.set_xlabel("Time since launch [years]")
    ax2.set_title("b) Orbital decay")
    ax2.set_ylim(140, 510)
    save(fig, "fig3_solar_cycle.png")


def fig_storm():
    q, s = load("storm_none"), load("storm_3day")
    # compare only while the storm case is still above 250 km; near re-entry
    # the difference blows up and hides the actual storm effect
    n = int(np.argmax(s["alt_km"] < 250.0))
    t = s["time_day"][:n]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax1.plot(q["time_day"], q["density_kgm3"], color=GREY, label="Quiet (Ap = 15)")
    ax1.plot(s["time_day"], s["density_kgm3"], color=C[1],
             label="3-day storm (Ap = 200)")
    ax1.axvspan(20, 23, color=C[3], alpha=0.15, lw=0)
    ax1.set_yscale("log")
    ax1.set_ylabel("Density at satellite [kg m$^{-3}$]")
    ax1.set_title("a) Density encountered by the satellite (start 400 km)")
    ax1.legend(loc="upper left", fontsize=8)

    ax2.plot(t, s["alt_km"][:n] - q["alt_km"][:n], color=C[1])
    ax2.axvspan(20, 23, color=C[3], alpha=0.15, lw=0)
    ax2.axhline(0, color=GREY, lw=1)
    ax2.set_ylabel("Altitude difference\nstorm - quiet [km]")
    ax2.set_xlabel("Time [days]")
    ax2.set_title("b) Extra altitude loss caused by the storm (until 250 km)")
    loss = q["alt_km"][int(np.searchsorted(q["time_day"], 23.0))] - \
        s["alt_km"][int(np.searchsorted(s["time_day"], 23.0))]
    ax2.annotate(f"{loss:.1f} km lost during the storm", (23, -loss),
                 xytext=(30, -loss / 2), color=INK2, fontsize=9,
                 arrowprops={"arrowstyle": "-", "color": GREY})
    save(fig, "fig4_geomagnetic_storm.png")


def fig_area():
    # a) instantaneous frontal area along two orbits (same formula as Fortran)
    a_body, a_panel = 0.6, 2.4
    u = np.linspace(0, 4 * np.pi, 600)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
    ax1.axhline(a_body + a_panel, color=C[0], label="Fixed, maximum")
    ax1.axhline(a_body + 2 / np.pi * a_panel, color=C[1], label="Fixed, orbit average")
    ax1.plot(u / (2 * np.pi), a_body + a_panel * np.abs(np.sin(u)), color=C[2],
             label=r"Time-varying, $\beta$ = 0$^\circ$")
    ax1.plot(u / (2 * np.pi), a_body + a_panel * np.abs(np.cos(np.radians(60)) * np.sin(u)),
             color=C[3], label=r"Time-varying, $\beta$ = 60$^\circ$")
    ax1.set_xlabel("Orbits")
    ax1.set_ylabel("Frontal area [m$^2$]")
    ax1.set_title("a) Ram-facing area of a sun-tracking satellite")
    ax1.set_ylim(0, 3.6)
    ax1.legend(fontsize=8, loc="lower right", ncol=2)

    # b) decay at moderate solar activity
    cases = [("area_fixed_max_moderate", "Fixed, maximum", C[0]),
             ("area_fixed_avg_moderate", "Fixed, orbit average", C[1]),
             ("area_varying_moderate", r"Time-varying, $\beta$ = 0$^\circ$", C[2]),
             ("area_varying_beta60_moderate", r"Time-varying, $\beta$ up to 60$^\circ$", C[3])]
    for name, lab, col in cases:
        d = load(name)
        ls = ":" if "varying_moderate" in name else "-"   # overlaps fixed avg
        ax2.plot(d["time_day"], d["alt_km"], color=col, ls=ls, label=lab)
    ax2.set_xlabel("Time [days]")
    ax2.set_ylabel("Altitude [km]")
    ax2.set_title("b) Altitude decay, F10.7 = 150")
    ax2.set_ylim(140, 510)
    ax2.legend(fontsize=8, loc="lower left")
    save(fig, "fig5_area_assumptions.png")


def read_summary():
    with open(ROOT / "results" / "summary.csv") as f:
        return {r["case"]: r for r in csv.DictReader(f)}


def fig_lifetime_bars(summary):
    levels = ["low", "moderate", "high"]
    modes = [("fixed_max", "Fixed, maximum"), ("fixed_avg", "Fixed, orbit average"),
             ("varying", r"Time-varying, $\beta$ = 0$^\circ$"),
             ("varying_beta60", r"Time-varying, $\beta$ up to 60$^\circ$")]
    x = np.arange(len(levels))
    w = 0.2
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for i, (mode, lab) in enumerate(modes):
        vals = [float(summary[f"area_{mode}_{lv}"]["lifetime_days"]) for lv in levels]
        bars = ax.bar(x + (i - 1.5) * w, vals, width=w - 0.02, color=C[i], label=lab)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.0f}", ha="center",
                    va="bottom", fontsize=7, color=INK2)
    ax.set_xticks(x, ["Low\n(F10.7 = 70)", "Moderate\n(F10.7 = 150)", "High\n(F10.7 = 230)"])
    ax.set_ylabel("Orbital lifetime from 500 km [days]")
    ax.set_title("Predicted lifetime: area assumption vs. solar activity")
    ax.grid(axis="x", visible=False)
    ax.legend(fontsize=8)
    save(fig, "fig6_lifetime_comparison.png")


def write_markdown(summary):
    lines = ["| Case | Start alt. [km] | Density model | F10.7 | Area model | Lifetime [days] |",
             "|---|---|---|---|---|---|"]
    for r in summary.values():
        area = r["area_mode"] + (f" (beta {float(r['beta_amp_deg']):.0f} deg)"
                                 if float(r["beta_amp_deg"]) > 0 else "")
        lines.append(f"| {r['case']} | {float(r['h0_km']):.0f} | {r['density_model']} | "
                     f"{r['f107']} | {area} | {float(r['lifetime_days']):.1f} |")
    (ROOT / "results" / "summary.md").write_text("\n".join(lines) + "\n")
    print("saved", ROOT / "results" / "summary.md")


if __name__ == "__main__":
    summary = read_summary()
    fig_density_profiles()
    fig_solar_decay()
    fig_solar_cycle()
    fig_storm()
    fig_area()
    fig_lifetime_bars(summary)
    write_markdown(summary)
