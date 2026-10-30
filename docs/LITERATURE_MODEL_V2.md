# Literature-led pilot V2 (2026-10-31)

## Publications (peer-reviewed, established publishers)

| Paper | DOI | Evidence transferred | Not transferable |
| --- | --- | --- | --- |
| Xiong, Liu & Schatz (2024), *International Journal of Smart and Nano Materials*, carboxylated MoS₂ lanthanide transport | [10.1080/19475411.2024.2387926](https://doi.org/10.1080/19475411.2024.2387926) | Hydration/dehydration + binding + channel spacing jointly affect transport; Sm maximum means no monotonic assumption | Acetate MoS₂ is not citrate MXene; no fixed Nd/Dy barrier numbers |
| He et al. (2024), *ACS Earth and Space Chemistry*, Nd/Dy/Yb citrate–goethite | [10.1021/acsearthspacechem.4c00032](https://doi.org/10.1021/acsearthspacechem.4c00032) | Citrate may promote surface sorption at pH<5 and inhibit it at higher pH by soluble complexation; SCM framework | Goethite constants do not transfer to MXene |
| Duvail et al. (2024), *JCTC*, lanthanide 12-6-4 force field | [10.1021/acs.jctc.3c01162](https://doi.org/10.1021/acs.jctc.3c01162) | Future atomistic MD validation against Ln–O distance, hydration free energy and CN | Does not give citrate–MXene PMF |
| Violet et al. (2024), *Nature Water*, selective binding site design | [10.1038/s44221-024-00279-6](https://doi.org/10.1038/s44221-024-00279-6) | Site identity and spacing for ultraselective membranes | Qualitative design principles, not calibrated citrate K |
| Liao et al. (2025), *Nature Communications*, 2D channel transport | [10.1038/s41467-025-61735-9](https://doi.org/10.1038/s41467-025-61735-9) | Hydration-shell distortions / friction alter transport | Studied mostly mono-/divalent ions, not Nd/Dy |
| Ding & Azimi (2024), *Scientific Reports*, Dy/Nd/Pr electrodialysis | [10.1038/s41598-024-62885-4](https://doi.org/10.1038/s41598-024-62885-4) | Extended Nernst–Planck transport after feed equilibrium | Uses EDTA and electrodialysis, not MXene |
| Perri et al. (2026), *ACS Omega*, citrate-cellulose Nd/Y | [10.1021/acsomega.5c07380](https://doi.org/10.1021/acsomega.5c07380) | pH/metal concentration and binding stoichiometry matter | Nd/Y not Nd/Dy; cellulose not MXene |
| Zhang et al. (2026), *Science*, rotatable NH4-MXene | [10.1126/science.aed8319](https://doi.org/10.1126/science.aed8319) | Controlling functionality mobility; NH4 residual confounding | Water transport only; no Nd/Dy evidence |
| Bulk citric acid pKa reference (ACS Omega 2019) | [10.1021/acsomega.9b00270](https://doi.org/10.1021/acsomega.9b00270) | pKa 3.13/4.76/6.40 **free citrate only** | Not grafted citrate; not soluble Ln–citrate formation constants |

## Physics-based substitute for MD/DFT: what this code actually solves

1. **Acid/base mass action** for free uncomplexed H3Cit/H2Cit-/HCit2-/Cit3- at 25°C. No metal complexation.
2. **Competitive Langmuir** occupancy for two ions with unknown surface K and pKa; requires independent fitting.
3. **Exact linear reversible sorption–diffusion breakthrough**: B≈S K C, retardation R=1+S K, ∂C/∂t=(D0/R)∂²C/∂x². Boundary mobile C(0,t)=C0, C(L,t)=0, C(x,0)=0.
4. **Steady permeation scenario**: P_i=Kentry_i Dmobile_i, PNd/PDy=exp[(Gentry,Dy-Gentry,Nd)/RT]*(DNd/DDy), input factors uncalibrated.

The analytical breakthrough function is
`J(t)/Jsteady = 1+2 Σ (-1)^n exp[-n²π² D0t/(R L²)]`;
`t50*D0/(R L²)≈0.138785`.

**Fundamental correction:** immobilized reversible binding delays breakthrough but does **not automatically decrease steady permeation** under these limits. The earlier 'toy hopping' conditional first-passage times are not the experimental steady permeability. Adsorption ≠ uptake ≠ lag time ≠ steady permeability.

## Parameter status & must-have measurements

- **Literature-anchored**: bulk acid pKa (3.13,4.76,6.40); model *form* for transport and binding.
- **Unknown, scanned only**: grafted site apparent pKa=5.5; K_Nd=10^4 M^-1; K_Dy/K_Nd=10^(−0.5..+0.5); max S*K_Nd=20; entry ΔG=−4..+4 kJ/mol; D_Nd/D_Dy=0.5..2.
- **Must measure**: site density, wet d-spacing and membrane thickness, citrate leachate, free metal activity, adsorption distribution, time-resolved breakthrough and eventual steady permeation; background nitrate & pH controls, OH precipitation checks, NH4 residual from HF2 etching.
- **Must retrieve**: Nd/Dy soluble citrate complex formation constants for all relevant protonation/stoichiometry species (e.g., 2023 Guseva et al. 1:2 citrate complexes, DOI 10.31857/S0044460X23120156) before PHREEQC soluble citrate extension.

Run:
```bash
pip install numpy pandas matplotlib scipy
python scripts/literature_pilot_v2.py
```

Full expanded model, generated CSVs/PNGs and 6 passing unit tests are included in the downloadable ZIP. No DFT, AIMD, MD, or PMF was run.