# CP2K: first-principles Nd³⁺/Dy³⁺ hydration starter

Status: **INPUT GENERATOR ONLY. No CP2K electronic energy or geometry optimization has yet been independently verified**.

## Why these basis and pseudopotentials?

- Lu, J.-B., Jiang, X.-L., Hu, H.-S., & Li, J. (2023), *Journal of Chemical Theory and Computation* **19**, 82–96, DOI [10.1021/acs.jctc.2c00922](https://doi.org/10.1021/acs.jctc.2c00922): specifically designed and validated trivalent lanthanide **4f-in-core LnPP2 GTH PBE** family (Nd, Dy: **GTH-PBE-q11**).
- CP2K official test `tests/QS/regtest-LnPP2/dycl3.inp`: `BASIS_MOLOPT_LnPP2`, `LnPP2_POTENTIALS`, `DZVP-MOLOPT-GTH-q11`.
- CP2K docs: [Run a First Calculation](https://manual.cp2k.org/trunk/getting-started/first-calculation.html), [Basis sets](https://manual.cp2k.org/trunk/methods/dft/basis_sets.html).

The original Zhang et al. (2026) *Science* SI used **VASP** for its NH₄-MXene DFT and **GROMACS** for MD. Our Ln hydration calculation is a **new independent starting model**, not a reproduction of the article.

## Generate 4 initial geometries

```bash
python simulations/cp2k_hydration/generate.py --run-type ENERGY
python -m unittest discover -s tests -v
```

The script generates `Nd(H₂O)₈³⁺`, `Nd(H₂O)₉³⁺`, `Dy(H₂O)₈³⁺`, `Dy(H₂O)₉³⁺`, using a cubic/approximately spherical water arrangement. **Starting Ln–O radii (Nd 2.50 Å, Dy 2.38 Å) are initialization guesses, not measured equilibrium values.** H directions initially point outward.

## Data setup

CP2K distribution should include `BASIS_MOLOPT_LnPP2`, `BASIS_MOLOPT`, `LnPP2_POTENTIALS` and `GTH_POTENTIALS`. Since both LnPP2 and ordinary O/H potentials are needed, make a local concatenated potential file:

```bash
cat /path/cp2k/data/LnPP2_POTENTIALS /path/cp2k/data/GTH_POTENTIALS > POTENTIALS_COMBINED
cp /path/cp2k/data/{BASIS_MOLOPT_LnPP2,BASIS_MOLOPT} .
cp simulations/cp2k_hydration/generated/nd_w8.inp .
OMP_NUM_THREADS=2 cp2k.ssmp -i nd_w8.inp -o nd_w8.out
```

This copies the normal CP2K distribution data for local computation only, rather than checking in a third-party dataset.

## Interpretation limits / next quality gates

1. First confirm CP2K startup, compatible data, successful **SCF convergence**, no numerical instability.
2. Initial `RUN_TYPE ENERGY` is just a single point on unoptimized, gas-phase geometry: does **not** determine preferred coordination.
3. Then run `--run-type GEO_OPT` and repeat from independent initial orientations; compare Ln–O radial distribution within each optimized structure. The bare +3 cluster in vacuum is strongly idealized.
4. Bulk hydration coordination probabilities and water exchange rates require explicit liquid water (validated force field or AIMD), not isolated clusters.
5. 4f-in-core is appropriate to many non-redox trivalent Ln structural questions, but cannot resolve 4f electronic fine structure or unusual redox.
6. Do **not** subtract absolute total energies between Nd and Dy or different stoichiometries to infer selectivity. To get hydration incremental free energies, use balanced exchange reactions and solvation/entropy corrections.
7. A future GROMACS hydration MD model needs validated Ln³⁺–water 12-6-4 or polarizable interactions and appropriate citrate–MXene cross terms; using generic LJ+charge without parameter validation is insufficient.

No unverified computations are labeled 'completed'.
