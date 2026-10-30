# CP2K starting-input snapshots (NOT computed DFT outputs)

These four `*.inp` and corresponding `*.xyz` files are reproducible initial configurations from `../generate.py --run-type ENERGY --cutoff 450`:

- Nd(H₂O)₈³⁺ / Nd(H₂O)₉³⁺
- Dy(H₂O)₈³⁺ / Dy(H₂O)₉³⁺

4f-in-core LnPP2 q11, PBE, DZVP-MOLOPT, cutoff 450 Ry, REL_CUTOFF 55 Ry, 26 Å nonperiodic box, net charge +3, 298 K **not prescribed** (single-point).

**No relaxation, SCF-converged energy, AIMD, DFT binding energy or transport PMF is contained here.**
The CP2K potentials and basis libraries are not vendored; provide those from a proper CP2K installation. GitHub job 41115528530 canceled due runner shutdown; the CP2K single-point must be manually opted into via workflow_dispatch, and output checked for convergence.
