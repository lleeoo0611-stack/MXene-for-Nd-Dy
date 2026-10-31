# V7 repair: CP2K checkpoints + Dy CN8/CN9 distinct initial states

## Ground truth from V6 (2026-10-31)

- GitHub run 41163885518: Nd(H2O)8^3+ converged (E = -171.69516239247577 Hartree; mean Nd–O = 2.58583227677 Å; preliminary PBE/LnPP2 220 Ry / 18 Å).
- Nd(H2O)9^3+, Dy(H2O)8^3+, Dy(H2O)9^3+ exited 124 after 1550 seconds. This is a time limit, not a geometry convergence result.
- V6 artifact output included last-structure XYZ, but its upload mask omitted genuine .restart / .wfn. Previous runs therefore only support reconstruction of positions; not Hessian, history or electronic-density reuse.
- Official CP2K docs: https://manual.cp2k.org/trunk/methods/restarting.html and https://manual.cp2k.org/trunk/methods/optimization/geometry.html .

## CP2K recovery path

1. Read last complete XYZ with exact atom count/species order using simulations/cp2k_hydration/checkpoint_v7.py.
2. Build fresh GEO_OPT with same geometry box, basis, pseudo, XC and numerical settings, but initialize from that last XYZ.
3. Configure MOTION/PRINT/RESTART on each GEO_OPT and DFT/SCF/PRINT/RESTART to emit geometric .restart and wavefunction .wfn; save those artifacts with the complete log and trajectory.
4. Separate GitHub job downloads the stage-1 artifact. Prefer real .restart positions and compatible .wfn via EXT_RESTART / WFN_RESTART_FILE_NAME and SCF_GUESS RESTART. If absent, fall back to last complete XYZ, explicitly marked positions-only.
5. CP2K exit 124 is a checkpointable interruption, NEVER a converged geometry. The run must pass force, displacement, energy and SCF validation before results are scientific.

Demonstration workflow: .github/workflows/nddy-checkpoint-v7.yml. Original checkpoint source: GitHub Actions 41163885518 artifact cp2k-Dy-H2O8-GEO_OPT-status-and-log. Both 240-second segments are a workflow **capability demonstration**, not guaranteed completion of a large DFT optimization.

## Distinct Dy CN8 vs CN9 initialization

- Previous V6 runs used identical initial bulk packing with independent velocity seeds, and essentially remained Dy CN=9 at 40-ps equilibration / 80-ps sampling.
- simulations/openmm_hydration/cn_initial_v7.py creates 8 cube-corner or 9 Fibonacci-sphere oxygen positions, at r(Dy–O)=0.245 nm. Water H atoms are positioned radially outward from the +3 center, preserving TIP3P OH and HH geometry.
- Each seed translates the entire bulk solvent configuration independently; waters within 0.52 nm are deleted and a new first-shell water cluster is inserted.
- The initialization is highly non-equilibrium. Records include requested CN, observed initial CN, minimized CN, initial XYZ, water count/density, CN time series, RDF histogram and per-frame first-minimum CN.
- Test workflow: .github/workflows/nddy-dy-cn8-cn9-v7.yml, 2 initial CN values × 2 independent solvent seeds, each 40 ps equilibration + 80 ps production. This is a **short diagnostic only**.
- scripts/analyze_dy_cn_initial_v7.py collects downloaded, actual summary and CN time series. Frame-wise fractions are correlated and MUST NOT be described as thermodynamic probabilities.

## What would count as a repair?

- CP2K test: second segment must load a genuine positional checkpoint and, if compatible wavefunction exists, SCF uses it. Artifacts and logs must demonstrate the mechanism. Convergence of Dy8/9 and Nd9 still requires further computational work.
- Dy test: both starts are demonstrated by saved XYZ and initial CN; whether each reaches a stable/equilibrated CN distribution needs extended sampling and at least one other validated lanthanide hydration force field.
- Production MD target: >= 0.5 ns equilibration + >= 2 ns sampling for at least 3 independent initial packings, starting from both CN8 and CN9; more sampling may be needed for rare exchange.

## Execution status (most recent check)

- The code, workflow definitions and regression tests have been committed to the V6 branch in PR #3.
- GitHub Actions runs 41211407165 and 41211428990 (and contemporaneous previously working V5/V6 workflows) failed within seconds, before readable job steps were available. Consequently **new V7 CP2K continuation and Dy 8/9 molecular dynamics results have not been validated**. Do not mistake committed inputs for completed calculations.
- This is a current execution blocker; use GitHub Actions run UI to inspect runner allocation, workflow permissions, billing and job setup. Do not assert a specific root cause from the unavailable logs.

Reference for the TIP3P 12-6-4 Lennard-Jones and C4 interactions: Li, Song & Merz, Journal of Physical Chemistry B (2015), DOI https://doi.org/10.1021/jp505875v .
