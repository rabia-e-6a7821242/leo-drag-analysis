#!/usr/bin/env bash
# Build the Fortran code, run every scenario and make the figures.
# Usage (from anywhere):  bash scripts/run_all.sh
# Uses .venv/bin/python if a virtual environment exists, otherwise python3.
set -euo pipefail

cd "$(dirname "$0")/.."

if [[ -x .venv/bin/python ]]; then PY=.venv/bin/python; else PY=python3; fi

echo "== Building Fortran programs =="
make

echo "== Density profiles =="
./bin/density_profile results/data/density_profiles.csv

echo "== Orbit decay scenarios =="
$PY python/run_scenarios.py

echo "== Figures =="
$PY python/plot_results.py

echo "== Validation =="
$PY python/validate.py

echo "== Comparison with NRLMSIS 2.1 (needs pymsis) =="
$PY python/compare_msis.py || echo "skipped"

echo "Done. Figures are in results/figures/"
