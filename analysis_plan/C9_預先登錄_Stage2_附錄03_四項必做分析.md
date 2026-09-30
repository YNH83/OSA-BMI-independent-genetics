# C9 預先登錄 Stage 2 附錄 03：四項必做分析（A6 到 A9）

**登錄時間** 見同夾 `SHA256_登錄時間戳.txt`。本檔在 A6 到 A9 任何分析執行之前寫成。
**使用者裁示** 2026-09-27：文獻調查（`09_圖表設計參考/C9_可套用分析_統整建議.md`）列出的四項必做全部執行。
**登錄時的狀態（須揭露）**
- v2 多變量 GWAS（F_NonBMI）仍在執行，尚未合併、尚未看過任何全基因體結果；先前只看過 v1 的 2,000 SNP 測試輸出。
- Pan-UKB 三個檔案下載中，內容未開啟。僅查證了 Pan-UKB 官方文件的共變項：年齡、性別、年齡乘性別、年齡平方、年齡平方乘性別、前 10 個 PC，**不含 BMI**（pan.ukbb.broadinstitute.org/docs/qc，2026-09-27 讀取）。
- 已知結果：A1 PASS（FinnGen 高血壓）、A2 FAIL（Evangelou 血壓，校正 BMI）、A3 反向效應旗標成立。A7 是針對 A1 已知型態的外部複製，照實揭露。
**固定不變** Stage 1 與 Stage 2 的資料、SNP 集合、LDSC 設定（sep_weights = TRUE、MHC 25 到 35 Mb 排除）、Cholesky 模型。

## A9 F_NonBMI GWAS 後處理（其餘三項的前提，先做）

- **QC**：剔除 userGWAS 回報 lavaan 警告或錯誤的 SNP；若輸出有平滑後 Z，剔除 |Z_Estimate − Z_smooth| > 1.96 的 SNP（Park 2024 Nat Genet；Foote 2025 Nat Genet）。
- **有效樣本數**：N_hat = 平均值 of 1 / (2p(1−p)SE²)，只用 MAF 0.1 到 0.4 的 SNP（Mallard 2022 公式）。LAVA、MAGMA、MR 一律用此單一 N_hat。
- **P3 判定不變**：依 Stage 1 與補充說明 01，閾值 5 × 10⁻⁸、500 kb 距離法。另報 2.5 × 10⁻⁸（兩因子多重比較）的位點數，只作敏感度，不改 P3 等級。
- **位點分類**：保留（FinnGen OSA 位點 500 kb 內有 F_NonBMI 全基因體顯著 SNP）、消失（無）、新增（F_NonBMI 位點 500 kb 內無 OSA 全基因體顯著 SNP）。

## A6 LAVA 局部遺傳相關（回應弱點 W3，並在局部層級扣除 BMI）

- **區塊**：LAVA 官方 2,495 個 LD 區塊（blocks_s2500_m25_f1_w200，GRCh37）；排除與 MHC（chr6 25 到 35 Mb）重疊的區塊。LD 參考為 1000G EUR（g1000_eur）。
- **樣本重疊**：以 LDSC 跨性狀截距矩陣輸入。
- **性狀**：OSA (FinnGen)、F_NonBMI (N_hat)、BMI (GIANT)、T2D (Xue)、HTN (FinnGen)、HTN (Pan-UKB, A7 資料)、打鼾、失眠、憂鬱、FEV1/FVC、睡眠時長、晝夜型。
- **步驟**：單變量局部 h2 P < 1 × 10⁻⁴ 的區塊才進雙變量（Werme 2022 慣例）；雙變量局部 rg 以每一性狀對內 FDR q < 0.05 為顯著，另報 Bonferroni。
- **A6a（描述）**：F_NonBMI 與 HTN 顯著局部 rg 的區塊數，對比 F_NonBMI 與 T2D。
- **A6b（主要）**：在 OSA、BMI、外部性狀三者皆有局部 h2 的區塊，以 LAVA partial correlation 估計 rg(OSA, X | BMI)，X 為 HTN 與 T2D。
  - **支持「高血壓在局部層級仍有不經 BMI 的共享」**：至少 1 個區塊的 partial rg(OSA, HTN | BMI) 達 FDR q < 0.05 且為正。
  - 同時報 T2D 的對應數字；不預設 T2D 必須為 0。
- **A6c（抵銷檢查）**：F_NonBMI 對全基因體 rg 接近 0 的性狀（FEV1/FVC、睡眠時長、晝夜型）。
  - **抵銷旗標**：同一性狀對同時出現顯著正、負局部 rg 各至少 1 個區塊。旗標成立時，主張 3 必須改寫為「全基因體淨相關為 0，但有局部方向相反的共享」，不得寫「無關」。

## A7 以 Pan-UKB 複製高血壓軸（回應弱點 W1、W2）

- **資料**：Pan-UKB 歐洲血統；收縮壓 4080（irnt）、舒張壓 4079（irnt）、高血壓 phecode 401（病例與對照數取自 Pan-UKB 官方 manifest 的 EUR 欄位）。來源不經 FinnGen，GWAS 未校正 BMI。
- **模型**：與 A1 相同（`07_twoaxis_test.R`），代謝軸 = T2D (Xue)；第二軸 = HTN (Pan-UKB) 為主要，SBP、DBP (Pan-UKB) 為次要。
- **判讀（與 A1 相同規則）**
  - **PASS（高血壓軸在 FinnGen 之外複製）**：cor(F_NonBMI, HTN_PanUKB) > cor(F_NonBMI, T2D) 且卡方差異 P < 0.05。
  - **FAIL**：A1 的結果可能來自與 OSA 共用的 FinnGen 登錄，主軸須改寫，高血壓不得作為主要發現；改以上呼吸道（打鼾）為主軸並回報使用者。
  - SBP、DBP 規則同上，屬次要；主要 PASS 而 SBP、DBP 皆 FAIL 時，寫「高血壓診斷層級複製，連續血壓未複製」。
- **限制（預先寫明）**：Pan-UKB 血壓未校正降壓藥；UKB 與 GIANT BMI、Xue T2D 有樣本重疊，由 LDSC 截距處理。

## A8 多變量 MR 與穩健性（把兩軸改以因果估計表達）

- **暴露**：OSA (FinnGen R12) 與 BMI (GIANT)。工具變數為兩者 P < 5 × 10⁻⁸ SNP 的聯集，以兩者較小 P 排序 clumping（r² < 0.001、10,000 kb、1000G EUR），去除回文 SNP。
- **結局**：HTN (Pan-UKB) 為主要（與 FinnGen OSA 樣本不重疊）；SBP (Pan-UKB)、T2D (Xue)、CAD (Aragam) 為次要；HTN (FinnGen) 只作敏感度（與 OSA 同庫重疊）。
- **方法**：MVMR-IVW 為主，報條件 F 值；OSA 條件 F < 10 標為弱工具。單變量 OSA → 結局另做 MR-PRESSO（離群校正）、Steiger、CAUSE（相關多效性）；BMI → OSA 與 BMI → 結局（GIANT 含 UKB，重疊）以 MRlap 校正樣本重疊。
- **判讀**
  - **支持「OSA 對高血壓有不經 BMI 的直接效應」**：MVMR 中 OSA 對 HTN (Pan-UKB) 的直接效應為正且 P < 0.05。
  - **支持兩軸對比**：上述成立，且 OSA 對 T2D 的直接效應 95% CI 含 0。只有前者成立時寫「OSA 對高血壓有直接效應，T2D 未區分」。
  - CAUSE 偏好共享模型（相關多效性）時，因果解讀降為「遺傳共享」，並在結果中揭露。
- 只作方向與因果順序證據，不換算臨床效果大小。

## 預寫的可能失敗陳述

1. A7 主要 FAIL：高血壓軸無法在 FinnGen 之外複製，主軸改為上呼吸道軸，高血壓降為次要描述並揭露。
2. A6b 無任何區塊達顯著：局部層級未見不經 BMI 的高血壓共享，高血壓軸只在全基因體平均層級成立。
3. A6c 抵銷旗標成立：主張 3 改寫為「無淨相關但有局部相反方向共享」。
4. A8 OSA 對高血壓直接效應不顯著或 CAUSE 偏好共享模型：不得使用因果語言。
5. A9 N_hat 過低或 QC 剔除比例 > 10%：F_NonBMI GWAS 的下游分析（MAGMA、LAVA、MR）加註功效限制。
