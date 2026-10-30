#!/usr/bin/env python3
"""V4 conditional Nd/Dy citrate model anchored to prior PHREEQC no-citrate output.

This is a hybrid SCREENING calculation, NOT a full PHREEQC coupled equilibrium,
not a surface complexation model, not DFT/MD and not a membrane flux prediction.
Inputs are explicitly versioned in configs/pilot_v4.json.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

from citrate_speciation_v3 import acid_fractions, gamma

R_KJMOL_K = 0.008314462618
REPO = Path(__file__).resolve().parents[1]

def config(path=None):
    path = Path(path) if path else REPO / "configs/pilot_v4.json"
    with path.open(encoding="utf-8") as f:
        return json.load(f)

def baseline_table(path=None):
    path = Path(path) if path else REPO / "results/phreeqc_pilot_grid.csv"
    lookup = {}
    with path.open(encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            key = (float(row["pH"]), float(row["NaNO3_background_mmolkgw"]),
                   float(row["each_Ln_mmolkgw"]))
            lookup[key] = {
                "Nd": float(row["Nd_free_fraction"]),
                "Dy": float(row["Dy_free_fraction"]),
                "SI_DyOH3_amorphous": float(row["SI_DyOH3_amorphous"]),
                "SI_DyOH3_crystalline": float(row["SI_DyOH3_crystalline"]),
            }
    assert len(lookup) == 72, "Missing V1 PHREEQC grid rows"
    return lookup

def solve_condition(pH, salt_mmolkgw, citrate_mM, base, config_data,
                    beta_shift_Nd=0.0, beta_shift_Dy=0.0):
    """Ltotal = Lfree + 2 sum Ln(Cit)(HCit); baseline odds from PHREEQC.
    'Other aqueous' is a LUMPED calibrated no-citrate nonfree species pool.
    Citrate reaction equilibrium uses source log(beta0), V3 Davies activity
    approximation. Ionic strength I=salt/1000 is held fixed (not self-consistent).
    All mol/kgw -> mol/L substitutions are dilute-solution approximations.
    """
    if citrate_mM < 0:
        raise ValueError("negative citrate")
    ion_total = config_data["each_Ln_total_mM"] / 1000.0
    I = salt_mmolkgw / 1000.0
    if I <= 0 or I > 0.5:
        raise ValueError("I out of Davies pilot bounds")
    fractions = acid_fractions(pH, I, pka=tuple(config_data["free_citrate_pKa"]))
    g3 = gamma(3, I)
    shifts = {"Nd": beta_shift_Nd, "Dy": beta_shift_Dy}
    Ltotal = citrate_mM / 1000.0

    def get_at(Lfree):
        hc = Lfree * fractions[2]
        cit = Lfree * fractions[3]
        metal = {}
        for el in ("Nd", "Dy"):
            f0 = base[el]
            if not 0.0 < f0 <= 1.0:
                raise ValueError("Invalid PHREEQC baseline fraction")
            noncitrate_odds = 1.0 / f0 - 1.0
            beta = config_data["literature_logbeta_mixed_citrate"][el] + shifts[el]
            mixed_odds = 10.0 ** beta * (g3 ** 2) * hc * cit
            den = 1.0 + noncitrate_odds + mixed_odds
            metal[el] = {
                "free": ion_total / den,
                "mixed": ion_total * mixed_odds / den,
                "other": ion_total * noncitrate_odds / den,
            }
        bound = 2 * sum(metal[el]["mixed"] for el in ("Nd", "Dy"))
        return metal, bound

    low, high = 0.0, Ltotal
    for _ in range(160):
        mid = (low + high) / 2
        _, bound = get_at(mid)
        if mid + bound > Ltotal:
            high = mid
        else:
            low = mid
    free_ligand = (low + high) / 2
    metal, bound = get_at(free_ligand)
    ferr = free_ligand + bound - Ltotal
    if abs(ferr) > 1e-13:
        raise ArithmeticError("Citrate mass balance failed")
    for el in ("Nd", "Dy"):
        if abs(sum(metal[el].values()) - ion_total) > 1e-13:
            raise ArithmeticError(el + " mass balance failed")

    def frac(el, group):
        return metal[el][group] / ion_total

    return dict(
        pH=pH, NaNO3_background_mmolkgw=salt_mmolkgw,
        each_Ln_mmolkgw=config_data["each_Ln_total_mM"],
        citrate_total_mM=citrate_mM,
        logbeta_shift_Nd_scenario=beta_shift_Nd,
        logbeta_shift_Dy_scenario=beta_shift_Dy,
        Nd_free_fraction=frac("Nd", "free"),
        Dy_free_fraction=frac("Dy", "free"),
        Nd_mixed_fraction=frac("Nd", "mixed"),
        Dy_mixed_fraction=frac("Dy", "mixed"),
        Nd_other_aq_fraction=frac("Nd", "other"),
        Dy_other_aq_fraction=frac("Dy", "other"),
        citrate_free_mM=free_ligand * 1000,
        citrate_bound_mM=bound * 1000,
        citrate_balance_error_M=ferr,
        free_Nd_over_free_Dy_if_equal_totals=metal["Nd"]["free"] / metal["Dy"]["free"],
        model="conditional_partial_v4",
    )

def calculate(config_data=None, baselines=None):
    c = config_data or config()
    b = baselines or baseline_table()
    scenarios = []
    for salt in c["NaNO3_background_mmolkgw"]:
        for ph in c["pH"]:
            key = (float(ph), float(salt), float(c["each_Ln_total_mM"]))
            if key not in b:
                raise ValueError(f"V1 baseline unavailable for {key}")
            for cm in c["citrate_total_mM"]:
                for sn in c["logbeta_shift_scenarios"]:
                    for sd in c["logbeta_shift_scenarios"]:
                        scenarios.append(solve_condition(ph, salt, cm, b[key], c, sn, sd))
    base = [row for row in scenarios
            if row["logbeta_shift_Nd_scenario"] == 0
            and row["logbeta_shift_Dy_scenario"] == 0]
    assert len(base) == 36 and len(scenarios) == 324
    transport = []
    for row in base:
        for energy in c["assumed_entry_GDy_minus_GNd_kJmol"]:
            for diffusivity in c["assumed_mobile_DNd_over_DDy"]:
                transport.append({
                    "pH": row["pH"],
                    "NaNO3_background_mmolkgw": row["NaNO3_background_mmolkgw"],
                    "citrate_total_mM": row["citrate_total_mM"],
                    "assumed_entry_GDy_minus_GNd_kJmol": energy,
                    "assumed_mobile_DNd_over_DDy": diffusivity,
                    "free_ion_feed_Nd_over_Dy": row["free_Nd_over_free_Dy_if_equal_totals"],
                    "apparent_JNd_over_JDy_if_only_free_Ln_transports":
                    row["free_Nd_over_free_Dy_if_equal_totals"]
                    * math.exp(energy / (R_KJMOL_K * c["temperature_K"]))
                    * diffusivity,
                })
    assert len(transport) == 324
    return base, scenarios, transport

def write_rows(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def assert_matching_csv(reference_path, computed, relative_tolerance=3e-9):
    with Path(reference_path).open(encoding="utf-8", newline="") as f:
        old = list(csv.DictReader(f))
    if len(old) != len(computed):
        raise AssertionError(f"Row count mismatch {len(old)} != {len(computed)}")
    for index, (expected, got) in enumerate(zip(old, computed), start=1):
        for field, val in got.items():
            if field not in expected:
                raise AssertionError(f"Missing column {field}")
            if isinstance(val, str):
                if expected[field] != val:
                    raise AssertionError((index, field, val, expected[field]))
            elif not math.isclose(float(expected[field]), val, rel_tol=relative_tolerance, abs_tol=3e-13):
                raise AssertionError((index, field, val, expected[field]))
    return True

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", default=str(REPO / "configs/pilot_v4.json"))
    p.add_argument("--baseline", default=str(REPO / "results/phreeqc_pilot_grid.csv"))
    p.add_argument("--out-dir", default="results/v4_regenerated")
    p.add_argument("--check-committed", action="store_true")
    args = p.parse_args()
    cfg = config(args.config)
    base, scenarios, transport = calculate(cfg, baseline_table(args.baseline))
    out = Path(args.out_dir)
    write_rows(out / "v4_citrate_baseline_36.csv", base)
    write_rows(out / "v4_beta_scenarios_324.csv", scenarios)
    write_rows(out / "v4_conditional_transport_324.csv", transport)
    if args.check_committed:
        for path, rows in (
            ("v4_citrate_baseline_36.csv", base),
            ("v4_beta_scenarios_324.csv", scenarios),
            ("v4_conditional_transport_324.csv", transport),
        ):
            assert_matching_csv(REPO / "results" / path, rows)
    print(f"V4 conditional: {len(base)} baseline, {len(scenarios)} beta scenarios, "
          f"{len(transport)} conditional-transport rows.")
    example = next(r for r in base if r["pH"] == 5.5
                   and r["NaNO3_background_mmolkgw"] == 10
                   and r["citrate_total_mM"] == 0.1)
    print("pH5.5 / NaNO3 10mM / citrate 0.1mM",
          "Nd free", f'{example["Nd_free_fraction"]:.6f}',
          "Dy free", f'{example["Dy_free_fraction"]:.6f}')
    print("WARNING: partial conditional hybrid; no comprehensive speciation, "
          "no adsorption/DFT/MD or measured membrane selectivity.")

if __name__ == "__main__":
    main()
