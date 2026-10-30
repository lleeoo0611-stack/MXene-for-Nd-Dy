# Wet-lab 前 pilot 與模型限制

## 1. 實際 PHREEQC 計算（72 組）

- 資料庫：USGS PHREEQC 3.8.6-17100 bundled `llnl.dat`；25°C，pH 4.0–8.0，0.5 間距。
- 雙金屬：Nd = Dy = 0.05 / 0.10 / 0.50 / 1.00 mmol/kg 水。
- 背景鹽：NaNO₃ = 10 / 100 mmol/kg 水；分析濃度含 Nd(NO₃)₃ 與 Dy(NO₃)₃ 硝酸根當量 `N(5)=NaNO3+3*Nd+3*Dy`。
- 假設 **無 CO₂／碳酸鹽、無 citrate、無 MXene**；pH 是指定輸入；未建立 `EQUILIBRIUM_PHASES`。
- 輸出：水溶液自由三價離子比例、硝酸根絡合、羥基絡合、Nd(OH)₃ 與 Dy(OH)₃ 的 crystalline/amorphous saturation indices。
- **SI>0 僅代表對指定固相的超飽和，不是沉澱量／沉澱速率／純化效率。** 非晶與晶態參照固相的 SI 可有很大差異。
- `llnl.dat` 的 NdNO₃²⁺ / DyNO₃²⁺ 形成常數差距明顯。背景鹽導致的 free-ion fraction 差異需交叉驗證，不宜直接引用為可靠的實驗平衡定量值。

## 2. 假說模型（非 atomistic）

- Energy sensitivity：`k_Nd/k_Dy = exp((ΔG‡_Dy−ΔG‡_Nd)/(RT))`，預設 298.15K、兩者相同前因子、單一有效活化步驟。
- Toy hopping：對給定離散位點自由能與過渡態能障，使用 detailed-balance-compatible rates 求 **已入膜、左側反射、右側吸收** 的 conditional MFPT。
- 其輸出不是 PMF，不包含擴散/分配係數、離子競爭、電中性、Donnan 效應、膜的實際透過通量或 citrate 的真實配位能。

## 3. 真實機制驗證前的門檻

- 測量 wet-XRD 層間距、膜厚與缺陷、表面官能基及其固定 citrate 的實際連接方式。
- 以 XPS N 1s／洗出液驗證 NH₄HF₂ 製程殘留 NH₄⁺是否構成混淆因素。
- 取得經交叉驗證的 Nd/Dy–citrate 水溶液物種穩定常數；**游離 citrate 的形成常數不能直接用於固定表面 citrate**。
- Nd/Dy³⁺ 原子尺度力場需再現配位數、Ln–O RDF、hydration dynamics（可能須極化或 12-6-4 模型）；再針對 2D 坐標 `G(z,CN_water)` 計算 PMF。

## 4. 文獻與原始材料

- 使用者提供：PHREEQC 3.8.6-17100 原始碼與 `llnl.dat`；未重新散布編譯執行檔或資料庫。
- Zhang et al., *Science* 394, 434–439 (2026) DOI 10.1126/science.aed8319，以及 Supplementary Materials，尤其 DFT/AIMD/MD 的 Materials and Methods（SI pp. 6–8）。
