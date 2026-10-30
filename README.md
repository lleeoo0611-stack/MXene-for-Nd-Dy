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
