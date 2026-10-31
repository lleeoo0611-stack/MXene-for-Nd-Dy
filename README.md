# Nd/Dy–citrate–MXene | Wet-lab 前 pilot 模擬

**計算狀態（2026-10-31）**

本專案研究是否能藉由 Ti₃C₂Tₓ MXene 奈米侷限與 citrate 配位差異，放大 Nd³⁺/Dy³⁺ 的水合層重組及傳輸選擇性。

## 證據等級（不可混用）

1. **實際執行：** PHREEQC 3.8.6-17100 + `llnl.dat`，18 組單／雙金屬基準，以及 **72 組雙金屬 pH × NaNO₃ × 濃度** 網格。計算**不含** carbonate、citrate、MXene、固相平衡或動力學。
2. **假說計算：** Python 活化能／自由能能井敏感度，以及條件式連續時間 Markov 跳躍模型。所有能量參數都是掃描假設，**不是** DFT、AIMD、PMF 或實測選擇性。
3. **尚未執行：** Nd/Dy–citrate–MXene 真正原子尺度 MD、DFT、PMF 或濕實驗。此 repo **不能**判定 Nd 或 Dy 真正較快穿透。

## 核心結果

PHREEQC 混合溶液，各離子 0.1 mmol/kgw，pH 6.5：
- NaNO₃ 10 mmol/kgw：free Nd³⁺ = 94.680%，free Dy³⁺ = 96.492%；Dy(OH)₃(am) SI = −2.4152。
- NaNO₃ 100 mmol/kgw：free Nd³⁺ = 79.700%，free Dy³⁺ = 95.386%；Dy(OH)₃(am) SI = −2.9672。
- pH 7.0 的兩種鹽濃度下，**crystalline** Dy(OH)₃ 的 SI 已 >0；而 **amorphous** 參照相在 pH 7.5 左右開始有超飽和風險。SI>0 **不代表實際沉澱一定發生**。

不同 nitrate 絡合常數是資料庫模型假設，必須以其他資料庫／實驗驗證。第一輪 wet lab 優先選 pH 4.5–6.5。

## 重現

需要 Python 3.10+，以及 **若重新解熱力學平衡**才需要 PHREEQC 3.8.6 executable 和相同的 `llnl.dat`：

```bash
pip install -r requirements.txt
python scripts/phreeqc_pilot.py --phreeqc /path/to/phreeqc --database /path/to/llnl.dat
python scripts/energy_sensitivity.py
python scripts/toy_transport.py
python -m unittest discover -s tests -v
```

預先算好的 `results/phreeqc_pilot_grid.csv` 可直接使用，不必安裝 PHREEQC。完整原始輸入、既有 18 組基準輸出、PNG 圖表、研究限制與 72 組資料也包含在同次交付的 ZIP 封存。

- 實驗安排：`docs/WET_LAB_PILOT_MATRIX.md`
- 科學限制：`docs/MODEL_LIMITATIONS.md`
- 數據：`results/phreeqc_pilot_grid.csv`

## 文獻定位

Zhang et al., *Science* **394**, 434–439 (2026), DOI [10.1126/science.aed8319](https://doi.org/10.1126/science.aed8319) 與 SI 研究的是 **NH₄⁺ rotatable crosslinker–water transport**；尚無 Nd/Dy–citrate 的直接結果。不可把該文的 NH₄⁺旋轉機制直接當成稀土水合離子的機制。

## Pilot V2: literature-grounded analytical models (2026-10-31)

Peer-reviewed source audit and equations: [docs/LITERATURE_MODEL_V2.md](docs/LITERATURE_MODEL_V2.md).

```bash
pip install -r requirements.txt
python scripts/literature_pilot_v2.py
python -m unittest discover -s tests -v
```

New V2 outputs under `results/literature_v2/`: citrate uncomplexed speciation; reversible binding breakthrough; site-affinity sensitivity; steady partition/mobility sensitivity. **All surface K, site pKa, binding/entry energies and mobility differences are scenarios, not DFT/MD or empirically fitted values.**

**Key correction:** fixed reversible sorption increases transient lag but may have no effect on steady permeability when mobile partition and diffusion are unchanged. Neither uptake nor conditional first-passage time equals membrane steady flux.

## CP2K Nd/Dy hydration starter

A reproducible LnPP2 4f-in-core trivalent lanthanide CP2K hydration cluster input generator and tests are in [`simulations/cp2k_hydration/`](simulations/cp2k_hydration/README.md). These are **initial geometry / input tests**, not evidence of hydrated-ion free energies or Nd/Dy membrane selectivity. The associated GitHub Actions CI checks input generation and runs a bounded single-point SCF smoke test.

## Pilot V3：Nd/Dy–citrate 文獻常數的部分配位平衡

- [完整方法、數值與重大限制](docs/CITRATE_SPECIATION_V3.md)。
- [可執行的純 Python 模型](scripts/citrate_speciation_v3.py) 與 [回歸測試](tests/test_citrate_speciation_v3.py)。
- 文獻 Table 1 取得混合 `Ln(Cit)(HCit)²−` 的 `logβ°`：Nd **13.05**，Dy **13.28**（ACS Earth and Space Chemistry, DOI [10.1021/acsearthspacechem.4c00009](https://doi.org/10.1021/acsearthspacechem.4c00009)）。
- pH 5.5、I=0.01 M、Nd/Dy 各 0.1 mM、游離 citrate 0.1 mM：模型計算 Nd 指定配合物 **19.9%**、Dy **29.6%**。
- **僅部分水溶液物種**（缺其他 LnCit/水解/硝酸根等配合物），也**不代表固定 citrate 的配位常數**；禁止將數字當完整物種預測或 Nd/Dy 穿透率。
- GitHub Actions 自動產生 18 組條件的 CSV 並檢查金屬／citrate 質量平衡。純 Python，無需 DFT/MD 商用授權。

```bash
python scripts/citrate_speciation_v3.py
python -m unittest discover -s tests -p test_citrate_speciation_v3.py -v
```

## V4：保留 V1–V3 並整合來源、設定與產出（2026-10-31）

- [V4 溶液物種數值模型與核查](docs/V4_PROVENANCE_AND_QC.md)；[機器可讀設定](configs/pilot_v4.json)。
- [V1 PHREEQC 無 citrate 基準](results/phreeqc_pilot_grid.csv) 72 組；[V2 重現曲線資料](results/v2_breakthrough_reproduced.csv)；[V3 重算結果](results/v3_citrate_aqueous_partial_reproduced.csv) 18 組。
- [V4 36 組基準](results/v4_citrate_baseline_36.csv)、[324 組 citrate β 敏感度](results/v4_beta_scenarios_324.csv)、[324 組條件式自由離子／膜進入假設](results/v4_conditional_transport_324.csv)。
- **V4 是 PHREEQC 無 citrate 校正分率 + 文獻單一 Ln(Cit)(HCit)²⁻ 配合物的混合近似；不是完整水溶液平衡、表面錯合模型、DFT、MD 或真實膜通量。**
- 前一次 CP2K SCF job 因 GitHub runner shutdown 被取消，**沒有已驗證的 DFT 能量**；已將耗資源的 CP2K 工作切換成手動觸發。

`python scripts/integrated_speciation_v4.py --out-dir build/v4 --check-committed` 可重算並核對 V4 數據。

- [V4 計算結果與實驗判讀](docs/V4_PILOT_DECISION.md)；[自由 Nd/Dy 的 SVG 圖](figures/v4_free_ions_pH5_5.svg)。
- [CP2K 四組未最佳化的 .inp / .xyz 結構快照](simulations/cp2k_hydration/snapshots/)，並由 generator regression test 驗證；**不是 DFT 能量結果**。

## V5: real atomistic DFT / MD hydration baseline (2026-10-31)

- [DFT CP2K Nd/Dy monohydrate input generator](simulations/cp2k_hydration/generate_dimer.py): one water, PBE/LnPP2 4f-in-core; [locally converged low-precision single-points CSV](results/v5_cp2k_monohydrate_local_smoke.csv). Absolute total energies across **different Ln species must not be subtracted**.
- [12-6-4 TIP3P OpenMM hydration MD](simulations/openmm_hydration/run.py) with published Li/Song/Merz parameters. CI workflow attempts **real MD trajectories and RDF**, and CP2K with official library data.
- [Methods, limitations, provenance and next scientific gates](docs/ATOMISTIC_V5.md).
- [Automated verification](.github/workflows/nddy-atomistic-v5.yml) artifacts retain raw logs/MD summaries. These are short **proof-of-execution** calculations, not MXene adsorption or Nd/Dy membrane selectivity predictions.

**V5 executions verified:** GitHub Actions [run #41161762479](https://github.com/lleeoo0611-stack/MXene-for-Nd-Dy/actions/runs/41161762479) completed five jobs successfully (input tests, Nd/Dy CP2K PBE/LnPP2 singlepoints, Nd/Dy OpenMM hydration short MD). Actual figures, RDF and metadata: [verified results](docs/V5_VERIFIED_RESULTS.md); [RDF comparison](figures/v5_md_Nd_Dy_RDF_0p02nm.svg). These are unoptimized DFT singlepoints and very short hydration MD, not citrate-MXene selective transport.
