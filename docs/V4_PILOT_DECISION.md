# V4 pilot decision summary (model-conditional, no actual membrane data)

The 36-case V4 conditional baseline combines nitrate + hydroxo + other **non-free Ln** odds from prior PHREEQC (no citrate) and one literature aqueous mixed citrate complex from Wen et al. 2024. Assumes each Ln total 0.10 mM, pH fixed and ionic-strength proxy from NaNO3.

At pH 5.5:

| NaNO₃ | Citrate total | Free Nd³⁺ | Free Dy³⁺ | Free Nd / free Dy |
|---|---:|---:|---:|---:|
| 10 mM | 0 mM | 95.963% | 98.907% | 0.970 |
| 10 mM | 0.1 mM | 77.181% | 69.363% | 1.113 |
| 10 mM | 0.5 mM | 0.0488% | 0.0288% | 1.698 |
| 100 mM | 0 mM | 80.371% | 96.569% | 0.832 |
| 100 mM | 0.1 mM | 66.096% | 67.029% | 0.986 |
| 100 mM | 0.5 mM | 0.3916% | 0.2312% | 1.694 |

**Interpretation:** citrate may reverse the free-ion feed ordering under the *selected partial-speciation assumptions*, particularly for 10 mM nitrate. It is **not** evidence for intrinsic CA-MXene uptake or Nd/Dy membrane permeability.

For pH 5.5 / 10 mM NaNO₃ / 0.1 mM citrate, ±0.5 log-unit *independent, uncalibrated* perturbations to each Ln's soluble mixed citrate formation constant yield free-Nd/free-Dy ratio spanning **0.656–1.688**. This is a **scenario range**, not a 95% confidence interval.

If citrate leaching occurs, apparent transport discrimination may change solely through feed speciation. Prioritize a citrate-free reference and a soluble-citrate-spiked blank before attributing effects to fixed CA-MXene.

**Figure:** [V4 conditional free trivalent fractions](../figures/v4_free_ions_pH5_5.svg), with all plotted values from [36-case CSV](../results/v4_citrate_baseline_36.csv).

This does not calculate stable adsorbed site populations or nanoconfinement PMF; these require controlled adsorption/transport experiments and optionally a validated ion force field, CP2K/QE DFT, and GROMACS/PLUMED.

The previous 18-case PHREEQC baseline *raw* output and older PNGs are not confirmed available in this repository; do not confuse this new reconstructed visualization with those prior originals.
