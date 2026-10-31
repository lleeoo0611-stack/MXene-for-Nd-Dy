# V7 verified Dy³⁺ CN8/CN9 initial-shell MD results (2026-10-31)

## Data provenance

GitHub Actions [rerun #41211990302](https://github.com/lleeoo0611-stack/MXene-for-Nd-Dy/actions/runs/41211990302), attempt 2: all **four actual OpenMM simulations completed**, using separate constructed starting shells and two seed-dependent bulk water placements. Raw per-trajectory RDF/CN/initial XYZ/metadata are retained in the run's artifacts. Extracted numeric table: [CSV](../results/v7_Dy_CN8_CN9_initial_120ps_replicas.csv).

| Starting CN | Seed | CN immediately after minimization | Production mean CN (RDF first minimum) | Sampled CN8 frames | Sampled CN9 frames | First Ln–O RDF peak |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 20261031 | 8 | 9.000 | 0 | 400 | 2.35 Å |
| 8 | 20261101 | 8 | 9.000 | 0 | 400 | 2.35 Å |
| 9 | 20261031 | 9 | 8.995 | 2 | 398 | 2.35 Å |
| 9 | 20261101 | 9 | 9.000 | 0 | 400 | 2.35 Å |

**Observed result**: Both initially 8-coordinate Dy³⁺ clusters became nine-coordinate by the first sampled frame following 40 ps equilibration, for this implementation and field parameterization. The CN9-start trajectories mostly remained nine-coordinate. **No unbiased equilibrium probability or water-exchange time constant** can be deduced from four 80 ps correlated trajectories.

## Scientific interpretation

This direct initial-structure comparison weakens the hypothesis that the older V5/V6 Dy CN9 result arose solely from initial nine-coordinate setup. It raises the possibility that (i) TIP3P 12-6-4 (Li/Song/Merz 2015) favors CN9 with the chosen parameters; (ii) water model + PME charged-background + finite box biases coordination; or (iii) prepared CN8 shell relaxes rapidly but rare reverse transitions need longer enhanced sampling. **Cannot distinguish without longer dynamics and force-field cross-checks.**

Reference: Li, Song & Merz, *Journal of Physical Chemistry B* **119** (2015), DOI [10.1021/jp505875v](https://doi.org/10.1021/jp505875v); comparison of Nd and Dy EXAFS/MD: [10.1021/acs.inorgchem.0c03438](https://doi.org/10.1021/acs.inorgchem.0c03438). This model has no citrate, no MXene and no explicit counterion; all water simulation box is fixed 2.5 nm, 298.15 K, PME +3 charge compensated by uniform plasma.

## QA correction: actual water count

The previous V7 version computed RDF normalization and `n_water` summary using the pre-shell solvent list size (`len(waters)`) rather than the actual regenerated shell size (`len(water_positions)`). Raw trajectory/coordinate counts, per-frame CN and RDF peak **position** are unaffected, but the reported RDF peak **height** and water density can be biased by about 1–3%. `run.py` has now been patched for all future jobs to normalize by the actual water count and explicitly report `water_density_molecules_per_nm3`. **Do not reuse pre-patch RDF peak heights without correcting the normalization.**

## Next computational gate

- Compare production trajectories at >=0.5 ns equilibration + 2 ns production × >=3 initial packings, both CN8 and CN9, with block error / decorrelation.
- Use at least one independent Ln-specific water force field, with water model matched to published parameters. Do not simply hand-tune C4 to enforce CN8.
- If apparent CN remains stuck or transitions are rare, model `G(CN_{\text{water}})` using umbrella sampling/MBAR or enhanced sampling.
- V7 CP2K Dy8 geometric checkpoint status must be checked separately and never misread as convergence.
