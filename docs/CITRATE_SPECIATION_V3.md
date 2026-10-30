# V3：Nd/Dy–citrate 的文獻約束部分平衡模型

## 原始數據

文獻：*ACS Earth and Space Chemistry* (2024), DOI [10.1021/acsearthspacechem.4c00009](https://doi.org/10.1021/acsearthspacechem.4c00009), Table 1。

| 反應 | Nd logβ° | Dy logβ° |
| --- | ---: | ---: |
| Ln³⁺ + Cit³⁻ + HCit²⁻ ⇌ Ln(Cit)(HCit)²⁻ | 13.05 | 13.28 |
| Ln³⁺ + OH⁻ ⇌ LnOH²⁺ | 5.82 | 6.41 |

**該文的 Nd/Dy 單 citrate（1:1）常數為空白**：沒有從相鄰稀土補造數值。

檸檬酸的**游離態** pKa = 3.13、4.76、6.40，參考 *ACS Omega* (2019), DOI [10.1021/acsomega.9b00270](https://doi.org/10.1021/acsomega.9b00270)。

## 這套程式已經求解的內容

- 25°C，pH 定義為 `−log10(aH+)`，外加 `I=0.01/0.10 M` 作為固定背景離子強度。
- 以 Davies approximation 估計活度係數。對三價離子尤其在 0.1 M 時精度有限。
- 酸鹼平衡分配 `H3Cit/H2Cit-/HCit2-/Cit3-`。以共同未配位 citrate 濃度計算配合物。
- 金屬質量平衡：`Ln_total = Ln3+ + LnOH2+ + Ln(Cit)(HCit)2-`。
- citrate 質量平衡：`Cit_total = Cit_free + 2*(Nd_mixed + Dy_mixed)`。
- 單調二分法解算 citrate 自由濃度，輸出 Nd、Dy 自由／水解／混合 citrate 形態的分率。

### 條件式結果（請勿當成完整水溶液物種分布）

| 條件 | Nd mixed | Dy mixed |
| --- | ---: | ---: |
| pH 5.5、I 0.01 M、[Nd]=[Dy]=0.1 mM、citrate=0.1 mM | 約 19.9% | 約 29.6% |
| 同條件、citrate=0.5 mM | 約 99.95% | 約 99.97% |

這表示文獻允許溶液中 Dy 形成上述**指定**配合物的優勢；**沒有**證實 citrate–MXene 對 Dy 的固定表面選擇性，更沒有推導出 Nd/Dy 的穿膜通量。

## 嚴重限制

1. 未加入 `LnCit`、`LnHCit`、`LnCit2`、`LnCitOH` 等其他 1:1／1:2 citrate 配合物。特別是條件式 99.9% 不可解釋為真實體系的百分比。
2. 未加入 nitrate/chloride/carbonate 及其絡合、反離子、競爭表面官能基、固相沉澱。
3. 未加入 MXene /固定 citrate，不能將 soluble logβ 當作固定表面配位常數。
4. 不曾執行 DFT、AIMD、MD、PMF；不可估計 hydration shell 的動態或真實膜選擇性。
5. 所有常數必須在一致的熱力學標準態、溫度與活度模型下使用；目前的 Davies 外推只是近似。

## 可重現執行

```bash
python scripts/citrate_speciation_v3.py
python -m unittest discover -s tests -p test_citrate_speciation_v3.py -v
```

需要 Python 3.10+，**不用安裝任何第三方 Python 套件**。

### 下一關：全物種與表面模型

若取得 Nd/Dy 其他 citrate 1:1／1:2 配位常數、質子化常數與來源，可擴充至 PHREEQC `SOLUTION_SPECIES` 或通用 mass action solver。接著由 ICP adsorption experiments 擬合固定 CA-MXene 的 site density 與 conditional binding constants；在此之前不要宣稱「citrate 固定化提升 Nd/Dy 的傳輸選擇性」。對後續 MD，可參考 Duvail et al. *JCTC* (2024), DOI [10.1021/acs.jctc.3c01162](https://doi.org/10.1021/acs.jctc.3c01162) 的 12-6-4 Ln 水合力場作驗證基準。
