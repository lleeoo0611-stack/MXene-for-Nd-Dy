# Proposed wet-lab pilot matrix (post-simulation)

## 起手式：12 個化學條件

各 Nd/Dy 0.10 mmol/kg 水，混合 1:1，在 **pH 5.5、6.5 × NaNO₃ 10、100 mM × 三種材料**：baseline MXene、負電惰性對照、citrate-MXene，共 **12 組**。每組至少三個獨立樣本，另需無材料空白；n=3 只適合初步可行性，不是最終 power analysis。

模型篩選：在 pH 6.5 及更低條件，研究所用的參照氫氧化物 SI 尚未顯示與 pH 7.5–8.0 一樣強的超飽和風險。pH≥7 需特別控制 Dy(OH)₃ 可能沉澱／空白損失。

## 先做 Aim 1 adsorption，再決定 Aim 2 膜

1. **空白／金屬回收**：沒有 MXene 的金屬溶液，在各 pH/鹽濃度追蹤 Nd、Dy 總量與 pH；分辨容器吸附和沉澱。
2. **自由 citrate vs immobilized**：檢查 citrate 溶出，以游離 citrate 補加空白；確保觀察的效應來自材料，而非溶液配合物。
3. **NH₄⁺ confound**：有 NH₄HF₂ 製程需檢查 N 1s、NH₄⁺洗出液；特別參考 Zhang et al. (2026) 的 NH₄⁺ crosslinking。
4. **ICP-OES QC**：雙譜線、基質匹配、加標回收、LOD/LOQ；同時測總金屬和物種相關的溶液條件。
5. **吸附與選擇性**：報告質量平衡及 `Δlnα = lnα_CA − lnα_control`，跨批次比較而不是只看最高 α。
6. **進入膜階段**：優先挑選 Aim 1 較有差異、無沉澱混淆的 1–2 個化學條件；測 wet XRD 的 effective interlayer gap，分開報告 partition K、permeability P、apparent D。

## 明確不能由此 pilot 推論

此次 PHREEQC 沒有 citrate，沒有膜或 transport；MFPT 與活化能掃描參數為假設值。Nd 比 Dy 快或相反，都尚未被模型證實。
