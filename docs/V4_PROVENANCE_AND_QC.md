# Pilot V4 — Nd/Dy–citrate–nitrate conditional equilibrium, reproducibility and audit

Date: 2026-10-31. This is a **numerical hybrid screening model**, not a thermodynamically complete aqueous-speciation solver and **not DFT / MD**. It reuses *existing PHREEQC output*; PHREEQC was not rerun for V4.

## File inventory and provenance

| Stage | Source code | Committed results | Status |
|---|---|---|---|
| V1 | \`scripts/phreeqc_pilot.py\` | \`results/phreeqc_pilot_grid.csv\`, 72 rows | previously computed with PHREEQC 3.8.6 + llnl.dat; no citrate |
| V2 | \`scripts/literature_pilot_v2.py\` | \`results/v2_*_reproduced.csv\` | rebuilt from identical published equations, toy parameters **not fitted** |
| V3 | \`scripts/citrate_speciation_v3.py\` | \`results/v3_citrate_aqueous_partial_reproduced.csv\`, 18 rows | re-evaluated V3 equations, checked against V3 successful Actions result |
| V4 | \`scripts/integrated_speciation_v4.py\` | 36 baseline, 324 beta scenarios, 324 conditional flux scenarios | computed with mass-balance bisection; V4 Python and CI match CSV |

V1 earlier single-metal/mixture 18-case baseline plots and prior ZIP are **not present in this Git repository**; no unverified path/file is presented as uploaded. V1 72-case grid is the available quantitative anchor. V2 archived curves are re-generated formula data, not a preserved original file hash. V3 archived CSV is re-computed from the tracked source rather than downloaded from Actions, and should only be called a reproduction after CI comparison.

CP2K: \`simulations/cp2k_hydration/generate.py\` generates **starting geometries and inputs only**. GitHub run [41115528530](https://github.com/lleeoo0611-stack/MXene-for-Nd-Dy/actions/runs/41115528530) was canceled after the runner received a shutdown signal during the SCF step. No independently verified converged DFT energy exists. The expensive CP2K job is now opt-in via manual workflow dispatch. A canceled runner is not evidence of SCF failure.

## Fixed literature sources

- **Wen, Liu, Wang, Zhao & Tang (2024)**, *ACS Earth and Space Chemistry* **8** 1048–1061, DOI [10.1021/acsearthspacechem.4c00009](https://doi.org/10.1021/acsearthspacechem.4c00009). Table 1 lists thermodynamic \`log β°\` for soluble **Ln(Cit)(HCit)²⁻**: Nd=**13.05**, Dy=**13.28**, hydroxides NdOH²⁺=5.82, DyOH²⁺=6.41. This is the source used by V3 and inherited by V4. The separate **He et al. (2024)**, DOI [10.1021/acsearthspacechem.4c00032](https://doi.org/10.1021/acsearthspacechem.4c00032), establishes goethite surface chemistry and is **not** the source of V3 mixed-citrate constants.
- Uncomplexed citric acid pKa: **3.13/4.76/6.40** from V2 references, not surface-bound citrate protonation constants.
- The old \`llnl.dat\` generated non-citrate free-ion fractions; constants from that database have **not** been independently audited against another model database.

## Hybrid equations (mass balances)

For fixed \`pH\`, \`I≈[NaNO3]\`, \`each total Ln=0.1 mM\`, and citrate total \`L_T\`:

1. Read \`f0_i(pH,salt)\` from the V1 no-citrate PHREEQC grid.
2. Empirical **nitrate+hydroxo+other** odds in no-citrate fluid: \`B_i=1/f0_i−1\`; these are **not** calibrated physical nitrate formation constants. The odds are held fixed on adding citrate.
3. Using the same Davies treatment and soluble Ln(Cit)(HCit)²⁻ formation constant as V3, calculate \`C_i(Lfree)=10^(logβ_i) × γ3² × [HCit²⁻] × [Cit³⁻]\`.
4. \`[Ln³⁺]_i = T_i/(1+B_i+C_i)\`, \`[Ln-mixed]_i=T_i×C_i/(1+B_i+C_i)\`, \`[Ln-other]_i=T_i×B_i/(1+B_i+C_i)\`.
5. Solve \`L_T = L_free + 2×([Nd-mixed]+[Dy-mixed])\` by monotonic bisection. Each row saves residual mass balance.
6. **Hypothetical permeation illustration** only if *free Ln³⁺ alone* enters channel:
  \`J_Nd/J_Dy = ([Nd³⁺]_feed/[Dy³⁺]_feed) × exp[(G_entry,Dy−G_entry,Nd)/(RT)] × D_Nd/D_Dy\`,
  taking equal total feed Nd=Dy. Parameters \`ΔG_entry=−2/0/+2 kJ mol−1\`, \`D ratio=0.5/1/2\` are **scenarios**, not measured or atomistically calculated.

Three citrate shifts \`−0.5,0,+0.5\` log units per metal independently are arbitrary **uncertainty/sensitivity stress tests**, **not** empirical uncertainty intervals.

V4 reuses the V3 citrate speciation mechanism but anchors its noncitrate fluid at V1 PHREEQC and should not be mistaken for *joint PHREEQC equilibrium*. Mobile nitrate ion activities, ligand effects on activity coefficients, OH solids, carbonate and other 1:1/1:2 citrate species are **not** recalculated. In particular, at citrate=0.5 mM the apparent ~100% mixed complex is conditional on excluding the missing complexes.

## Computational parameters

Single source of truth: \`configs/pilot_v4.json\`.
- **pH** 4.5, 5.5, 6.5; **NaNO3** 10 and 100 mmol/kgw; **Nd/Dy** each 0.1 mmol/kgw (treated approximately as mmol/L).
- **Citrate total** 0, 0.05, 0.1, 0.2, 0.5, 1 mM.
- **V4 grid** 3×2×6 = 36; 3×2×6×3²=324 beta perturbations; 36×3×3=324 hypothetical feed+membrane scenarios.
- Constant I approximated by background \`NaNO3\`; Davies is limited for Ln³⁺ at I=0.1 M.

## Run, verify, compare

No third-party dependencies for V3/V4:
\`\`\`bash
python scripts/citrate_speciation_v3.py --output build/v3.csv
python scripts/integrated_speciation_v4.py --out-dir build/v4 --check-committed
python -m unittest discover -s tests -p 'test_integrated_v4.py' -v
\`\`\`

GitHub workflow \`.github/workflows/nddy-v4.yml\` regenerates V4 and checks all numerical columns versus committed reference CSV within floating-point tolerance. Do not claim completed CI before its job result is inspected.

## Experimental gate

Select pH 5.5 and 6.5, 10/100 mM NaNO3, no-MXene blank, untreated MXene, MXene/PSS negative charge control, and CA-MXene, with and without soluble citrate leachate. Measure pH, citrate/TOC, ICP mass balance and precipitate-free filtrate. Any selectivity claim must separate soluble citrate speciation from surface immobilization, equilibrium adsorption, transient breakthrough and true steady-state permeation.
