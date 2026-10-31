# V5 atomistic hydration shakedown — exact scope and scientific limits

## DFT: CP2K 2026.2

Uploaded static executable CP2K 2026.2 (GCC 15.3, OpenMP) executed in a Linux x86-64 container. **Actual successful initial local tests**:

| species | INPUT assumption | cutoff | cell | EPS_SCF | total electronic E (hartree) | status |
| --- | --- | --- | --- | --- | ---: | --- |
| Nd³⁺–H₂O | Nd–O = 2.50 Å; unoptimized gas phase | 180 Ry | 12 Å | 5×10⁻⁵ | -50.558868654924169 | SCF converged |
| Dy³⁺–H₂O | Dy–O = 2.38 Å; unoptimized gas phase | 180 Ry | 12 Å | 5×10⁻⁵ | -54.332836524869862 | SCF converged |

- PBE, GPW, 4f-in-core LnPP2 q11, DZVP-MOLOPT-GTH, nonperiodic Martyna–Tuckerman Poisson solver, initial 1 water. Source 4f-in-core pseudopotential/basis: Lu et al., *JCTC* 2023, DOI 10.1021/acs.jctc.2c00922.
- These are **different elements with different pseudopotentials**, so subtracting absolute total energies is **physically invalid** as a comparison of hydration preference. Neither is an adsorption energy, optimized bonding distance, bulk hydration free energy, or real water coordination prediction.
- 180 Ry cutoff and 5e-5 SCF tolerance are intentionally **low accuracy smoke tests**. Use numerical convergence (cutoff 250/350/450 Ry; cell 12/16/20 Å; EPS_SCF ≤1e-6), geometry optimizations and balanced reactions before scientific interpretation.
- Local minimal basis/potential subset was manually prepared from official CP2K v2026.2 LnPP2_POTENTIALS, BASIS_MOLOPT_LnPP2, GTH_POTENTIALS, BASIS_MOLOPT; use **unmodified CP2K distributions** for production. The GitHub CI uses conda-forge CP2K data, not that manually truncated local subset.
- Local unoptimized positions and full output logs are available in the associated exported archive when separately supplied; the repo records an easily audited energy summary, not a claim of completed production DFT.

## MD: OpenMM as free 12-6-4 platform

- Real atomistic MD implementation: [../simulations/openmm_hydration/run.py](../simulations/openmm_hydration/run.py).
- Parameters **Li, Song & Merz 2015**, *J. Phys. Chem. B* 119, 883–895 DOI [10.1021/jp505875v](https://doi.org/10.1021/jp505875v) — TIP3P column: Nd³⁺ Rmin/2=1.724 Å, ε=0.15486311 kcal/mol, C4=213 kcal·Å⁴/mol; Dy³⁺ Rmin/2=1.637 Å, ε=0.09900804 kcal/mol, C4=207 kcal·Å⁴/mol.
- TIP3P water: O charge −0.834e, H +0.417e; O sigma 0.315075 nm; epsilon 0.636386 kJ/mol.
- 12-6-4 = standard PME Coulomb plus 12-6 LJ plus explicit \`−C4/r⁴\` on Ln-water oxygen pairs only. \`C4[kJ nm^4/mol]=C4[kcal Å^4/mol]×4.184×10⁻⁴\`.
- One metal ion and periodic TIP3P water with net +3; PME implicitly neutralizes charge. No explicit anions/surface/ligands; the small box causes finite-size effects. Thus **no hydration free energy, solvation equilibrium, or actual selectivity** can be inferred.
- Quick runs are deliberately short (default 4 ps relaxation + 8 ps sampling), much shorter than published 500 ps equilibration + 2 ns sampling. RDF and r<0.35 nm CN are diagnostic only. No water exchange, radial dynamical residence times, citrate or MXene membrane simulation.
- **Duvail et al. JCTC 2024 10.1021/acs.jctc.3c01162** is a *different* 12-6-4 parameter set. The current code uses the separately verified older Li–Merz TIP3P values because full appropriate pair parameters are publicly accessible.
- Use ≥3 independent trajectories, larger boxes + finite-size checks, ~0.5 ns equilibration, ~2 ns production and experimental Ln–O distances/CN comparisons for a defensible baseline.

## Next research gates

1. CP2K more converged Nd/Dy monohydrate **within-ion** energy vs distance, then optimize 8/9 H₂O and compare balanced water-addition reactions or statistical solution free energy.
2. Validate MD Nd–O/Dy–O RDF and CN vs lanthanide literature; test force-field sensitivity and explicitly neutralized aqueous backgrounds.
3. Only then construct citrate surface sites and Ti₃C₂Tₓ MXene slit geometry (surface chemistry, charges, wet interlayer gap).
4. Run PMF over ion axial position and first-shell water CN; obtain error bars and independent replicas. Cannot conclude membrane alpha from a short hydration benchmark.

**No physical inference that Nd diffuses/absorbs more strongly than Dy until controls are run.**
