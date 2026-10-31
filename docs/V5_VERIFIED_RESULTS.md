# V5 verified atomistic pilot results (2026-10-31)

## 1. REAL CP2K calculations

**Both ions succeeded locally and on GitHub Actions**. See [official-CP2K CI log artifacts](https://github.com/lleeoo0611-stack/MXene-for-Nd-Dy/actions/runs/41161762479).

| Geometry | CI PBE/LnPP2 250 Ry, 12 Å, fixed geometry E (Ha) | SCF |
| --- | ---: | --- |
| Nd(H₂O)³⁺, initial Nd–O 2.50 Å | −50.560063222794035 | converged |
| Dy(H₂O)³⁺, initial Dy–O 2.38 Å | −54.334334295052095 | converged |

These are one-water **gas-phase, unoptimized, single-point electronic total energies**, using distinct pseudopotentials. **Do not directly subtract Nd total E from Dy total E** or infer selectivity, hydration free energy or their preferred coordination numbers.

Local static binary low-cutoff probes with 180 Ry (not rigorously grid-converged) yielded Nd −50.558868654924169 Ha, Dy −54.332836524869862 Ha. Additional local *within-element* fixed-water bond-length scans: [CSV](../results/v5_cp2k_local_bond_distance_scan.csv). Both energies become lower on shortening the initial Ln–O distance by 0.10 Å, but a minimum is **not bracketed or optimized**, and numerical sensitivity is still significant. **No optimized coordination preference is established.**

## 2. REAL OpenMM 12-6-4 hydration MD

Both Nd and Dy actual MD jobs **succeeded** in [workflow 41161762479](https://github.com/lleeoo0611-stack/MXene-for-Nd-Dy/actions/runs/41161762479). Uses Li, Song & Merz (*J. Phys. Chem. B*, DOI 10.1021/jp505875v) TIP3P-specific Ln³⁺ 12-6-4 parameters. Trajectories/summary processed on GitHub Actions and original full 100-bin CSVs retained in workflow artifacts.

| Statistic | Nd | Dy |
| --- | ---: | ---: |
| Water count | 500 | 500 |
| Box side | 2.5 nm | 2.5 nm |
| Temperature target | 298.15 K | 298.15 K |
| NVT equilibration | 4 ps | 4 ps |
| NVT sampling | 8 ps | 8 ps |
| Frames | 160 | 160 |
| First Ln–O RDF peak | 2.45 Å | 2.35 Å |
| Mean count water O ≤0.35 nm | 9.05 | 9.00 |

- RDF first peak difference **0.10 Å** is consistent in sign with lanthanide contraction, but based on **one short trajectory**, cannot support a quantitative equilibrium CN comparison.
- Averaged RDF [50-bin CSV](../results/v5_openmm_12-6-4_rdf_comparison_0p02nm.csv) and [plotted SVG](../figures/v5_md_Nd_Dy_RDF_0p02nm.svg) are computed from the original GH artifacts' pair histogram counts. The original artifact CSVs have 100 bins, Δr=0.01 nm.
- [Nd full MD artifact](https://github.com/lleeoo0611-stack/MXene-for-Nd-Dy/actions/runs/41161762479/artifacts/13026884841); [Dy full MD artifact](https://github.com/lleeoo0611-stack/MXene-for-Nd-Dy/actions/runs/41161762479/artifacts/13027531697). Both accessible to authorized GitHub accounts. Local exported copies can also be supplied as ZIP.

## 3. Limits and next step

- **No citrate, MXene or actual membrane simulations**; no permeation flux, PMF, surface attachment energies or water exchange kinetics established.
- No explicit counterions; PME neutralizing plasma +3, finite box, 8 ps dynamics; note published validation studies used orders of magnitude longer runs.
- Next DFT: verify basis/grid/cell convergence and relax geometry, then Nd/Dy isolated 8-/9-water structures; compute *balanced within-ion* water addition energies and explicit solvation.
- Next MD: three independent seeds, ≥500 ps equilibration + ~2 ns production, density/box-size tests, CN from RDF first minimum not arbitrary radius and standard errors; then add validated citrate–MXene cross-interactions, membranes and PMF.
- Current DFT and MD establish **functional computational tools and initial hydration structure trends**; they do not establish CA-MXene Nd/Dy selectivity.
