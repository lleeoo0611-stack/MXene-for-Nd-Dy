# V6 — CP2K 八／九配位最佳化與 Nd/Dy MD 水合基準驗證

## 主要目的與證據等級

**Goal 1: DFT.** 用相同 PBE、4f-in-core LnPP2 (Lu et al., JCTC 2023, DOI 10.1021/acs.jctc.2c00922)、基組與收斂閾值最佳化 Nd(H₂O)₈³⁺、Nd(H₂O)₉³⁺、Dy(H₂O)₈³⁺、Dy(H₂O)₉³⁺。此時只比較**孤立氣相局部最低能構型**，**不得**判定水溶液八／九配位占比、傳輸選擇性或水合自由能。

**Goal 2: MD.** 使用 V5 已實測可執行的 OpenMM TIP3P 12-6-4 Li–Song–Merz (2015) 參數，獨立三個 random seeds 重複各 Nd/Dy 短時間 MD。保存逐幀 CN、RDF 原始 bin 計數、RDF 第一谷界線、六時間區塊均值，並檢查重複間變異；每個 replica 是一次獨立 Langevin 軌跡。

### V6 自動工作之設計（不是已完成的聲明）

- DFT：共 4 個獨立工作，\`GEO_OPT\`，18 Å 無週期盒，PBE/LnPP2 q11，220 Ry 網格 cutoff，REL_CUTOFF 45 Ry，EPS_SCF=1e-5、BFGS-like \`LBFGS\`（所用 generator 原設置為 BFGS? 實際取自 V5：\`OPTIMIZER BFGS\`），最多 55 結構步；收斂條件 MAX_FORCE 4.5e-4、RMS_FORCE 3e-4 Hartree/Bohr，MAX_DR 3.0e-3、RMS_DR 1.5e-3 Bohr。**只有 CP2K 輸出包含 GEOMETRY OPTIMIZATION COMPLETED、所有強制條件通過且沒有 SCF failure，才算完成。** 這是 preliminary screening resolution；最後要用 ≥300/450 Ry、較大盒與 EPS_SCF ≤1e-6 對最佳化後構型重算並比較力、能量、Ln–O。
- MD：Nd/Dy × 3 seeds (20261031,20261101,20261102)，NVT 298.15K，2.5nm 邊長 ~500 TIP3P waters，12-6-4 PMF parametrization，**40 ps equil + 80 ps sample** 每個 replica。輸出每個 Ln 的 Ln–O RDF、逐幀 0.35 nm CN、估計 RDF 首個谷值及 CN、seed、短區塊平均、軟體版本。
- 40+80 ps 屬於**較 V5 進步的複本篩選**，不是充分平衡；用於發現數值不穩定、檢查力場是否可能違背文獻 CN 趨勢。要宣稱「水合基準驗證完成」仍需 **≥0.5 ns equil + ≥2 ns sample × ≥3 replica**，有限大小與力場交叉檢驗，RDF first-min 積分、CN8/9 轉換、CN block error、收斂趨勢與文獻差距，並避免在未觀察事件時報告 water exchange time。
- 帶總電荷 +3 的 PME 使用 uniform neutralizing plasma，沒有 nitrate、citrate、MXene 或 explicit counterion；僅適合初步水合診斷。

## 文獻比較（有來源的參考，不能移植為模型參數）

1. Neutron scattering: Nd³⁺ first-shell hydration number **8.5±0.2**, Dy³⁺ **7.4±0.5**：*Journal of Non-Crystalline Solids* 1986, [ScienceDirect](https://www.sciencedirect.com/science/article/pii/0022508886904169). 這是所列實驗體系，依方法/濃度可能與現今 EXAFS 有差異。
2. Shiery et al. (2021), *Inorganic Chemistry* **60**, 3117–3130, DOI [10.1021/acs.inorgchem.0c03438](https://doi.org/10.1021/acs.inorgchem.0c03438)：Nd/Dy 水合離子 AIMD + EXAFS，顯示動態且不對稱的第一配位殼。
3. Alanís-Manzano et al. (2024), *ACS Omega* **9**, 50978–50991, DOI [10.1021/acsomega.4c04947](https://doi.org/10.1021/acsomega.4c04947)：r²SCAN-3c explicit-water cluster 模型顯示 Nd（水相）主要 CN9，重鑭系較偏 CN8，並需考慮 cluster 多構型與振動頻率才可推判穩定結構。
4. Li, Song & Merz (2015), *J. Phys. Chem. B* **119**, 883–895, DOI [10.1021/jp505875v](https://doi.org/10.1021/jp505875v)：12-6-4 離子/水 TIP3P 參數。
5. CP2K 官方 [GEO_OPT 收斂判準](https://manual.cp2k.org/trunk/methods/optimization/geometry.html)：同時檢查 maximum/RMS force 與 displacement，**MAX_ITER reached ≠ converged**。

## 正確的 8/9 配位比較

單獨比較 \`E[Ln(H₂O)9³⁺]-E[Ln(H₂O)8³⁺]\` 不平衡。可定義：

\`ΔE_add,gas(Ln) = E[Ln(H₂O)9³⁺] - E[Ln(H₂O)8³⁺] - E[H₂O]\`

僅是氣相局部勢能，不是溶液八/九配位熱力學分配比。必須納入水的化學位能、振動/轉動/平動自由能、外部水層與溶劑環境，且基組誤差要可控。跨元素的裸總能差更不能比較。

## Reproduction

\`\`\`bash
python -m unittest discover -s tests -p 'test_hydration_v6.py' -v
python simulations/cp2k_hydration/generate_v6.py --out build/v6_dft
pip install openmm==8.4.0 numpy
python simulations/openmm_hydration/run.py --ion Nd --seed 20261031 --equil-steps 20000 --sample-steps 40000 --interval 100 --out build/v6_md
\`\`\`

真正分子動力學數據與 CP2K 結構、日誌由 GitHub Actions artifacts 保存，程式本身與輸入/模型參數保存在 Repo。**該 CI 工作若逾時或 fail，必須如實記錄未完成。**

## 2026-10-31 computation correction

Six 40+80 ps MD jobs on [Actions run 41163249828](https://github.com/lleeoo0611-stack/MXene-for-Nd-Dy/actions/runs/41163249828) succeeded; first-min CN across independent velocity seeds Nd=9.035 ± 0.015 (sample SD), Dy=8.999 ± 0.004. All six started from identical solvent packing and Dy fails to show expected dominant eight coordination, so **hydration baseline NOT validated**. Full original records in GH artifacts; [PR discussion](https://github.com/lleeoo0611-stack/MXene-for-Nd-Dy/pull/3) has each seed and source. Nd9 GEO_OPT runner failed due MPICH/UCX `ud_verbs/mana_0` initialization, not demonstrated nonconvergence. Separate CP2K Nd–water **UCX_TCP smoke test passed** [run 41163608195](https://github.com/lleeoo0611-stack/MXene-for-Nd-Dy/actions/runs/41163608195); GEO_OPT now specifies `UCX_TLS=tcp,self`, `UCX_NET_DEVICES=all`. Expensive V5 and V6 MD repeat jobs are now opt-in via `workflow_dispatch` to avoid duplicate GitHub Actions usage after establishing six baseline traces.
