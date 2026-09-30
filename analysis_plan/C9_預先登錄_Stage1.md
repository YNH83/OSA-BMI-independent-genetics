# C9 預先登錄 Stage 1：OSA 的非肥胖遺傳成分

**登錄時間** 2026-09-26（在任何 LDSC、Genomic SEM 結果產生之前寫成；本檔 SHA256 記錄於同夾 `SHA256_登錄時間戳.txt`）
**門檻來源** `02_題目評估/09262026 開放數據GWAS題目評估與執行建議_v2.docx` 第五節草案，原樣照錄。
**使用者裁示** 2026-09-26 使用者指示「開始進行 C9 和 C1b」。三個主要門檻未經使用者逐條修改；若使用者在停損測試結果出來前要求調整，以 Stage 1 修正版另存，舊版保留。

## 一、研究問題

扣除 BMI 的遺傳成分後，OSA 還剩多少遺傳力？剩下的成分是否仍與心代謝疾病共享遺傳？
參考文章（Liao CC 等，DMS 2026）的共同因子由 BMI 主導（808 個領頭變異中 697 個已在 BMI 顯著），本題以 GWAS-by-subtraction 直接拆開。

## 二、資料（已下載並驗證，見 `04_公開資料庫原始檔/`）

| 角色 | 資料 | 來源 |
|---|---|---|
| 主要 OSA | FinnGen R12 `G6_SLEEPAPNO`，56,885 例 / 441,137 對照 | FinnGen R12 |
| 主要 BMI | GIANT Yengo 2018，約 70 萬人 | GWAS Catalog GCST006900 |
| 心代謝外部相關 | T2D Xue 2018 (GCST006867)；HTN FinnGen `I9_HYPTENSESS`；CAD Aragam 2022 (GCST90132314) | |
| 虛無對照 | 痛風 FinnGen `GOUT_STRICT` | |
| 三方比對 | BMI 校正 OSA (Campos 2023, GCST90267363)；BMI < 30 OSA (Strausz 2025, GCST90624797) | |
| 敏感度 | OSA `G6_SLEEPAPNO_INCLAVO`；BMI UKB `ukb-b-19953` (OpenGWAS)；WHRadjBMI (GCST008999) | |
| 獨立歐洲複製 | Campos 2023 OSA (GCST90267361)，是否含 FinnGen 樣本須讀原文確認 | |

## 三、固定的技術設定

- SNP：HapMap3 (w_hm3.snplist)，排除 MHC (chr6 25 到 35 Mb)，MAF ≥ 0.01，INFO 以來源檔可得者為準。
- LD 參考：1000G Phase3 EUR LD 分數（本機 Microbiome-Epi `raw data/ldsc_ref`）。
- 二元性狀 liability 轉換：樣本盛行率當族群盛行率（OSA 0.114、T2D、HTN 依各自樣本）。三個主要門檻都是比例或遺傳相關，不受此假設影響。
- 有效樣本數：二元性狀用 4 / (1/Ncase + 1/Ncontrol)。
- 模型：Genomic SEM Cholesky 分解。BMI 因子負荷在 BMI 與 OSA；NonBMI 因子只負荷在 OSA；兩因子正交。估計法 DWLS。
- 座標：FinnGen 為 hg38，GIANT 為 hg19；LDSC 與 Genomic SEM 以 rsID 合併，不做 liftover。

## 四、停損測試（5 分鐘，最先跑）

- OSA liability SNP 遺傳力 z < 5 → 退場。
- rg(OSA, BMI) > 0.85 → 退場（非 BMI 成分太小）。
- 兩者皆未觸發才往下跑。

## 五、主要門檻（三個，判讀先寫死）

| 代號 | 指標 | PASS | PARTIAL | FAIL |
|---|---|---|---|---|
| P1 | NonBMI 成分佔 OSA SNP 遺傳變異的比例 | ≥ 0.60 | 0.40 到 0.60 | < 0.40 |
| P2 | rg(NonBMI-OSA, T2D) 相對 rg(OSA, T2D) 的衰減 | ≥ 50% | 25% 到 50% | < 25% |
| P3 | FinnGen OSA 單一 GWAS 的全基因體顯著位點中，在 NonBMI GWAS 仍 P < 5 × 10⁻⁸ 的比例 | ≥ 30% | 15% 到 30% | < 15% |

P1 理論值為 1 − rg²（rg 約 0.50 時約 0.75）。P2 的方向：衰減 ≥ 50% 表示 OSA 與 T2D 的共享遺傳大多經由肥胖；衰減 < 25% 表示存在不經肥胖的 OSA 心代謝路徑。兩種結果都照實寫，任何一種都不改框架。

## 六、競爭假說

- H1 肥胖主導：P1 < 0.40，心代謝相關幾乎消失。
- H2 顱顏與神經肌肉路徑：NonBMI 成分富集在顱顏骨骼、腦與呼吸控制。
- H3 脂肪分布而非總量：rg(NonBMI-OSA, WHRadjBMI) > 0.20。

## 七、虛無對照與敏感度

- 虛無：(1) 痛風扣 BMI 做同樣 subtraction，殘餘成分不應富集於顱顏與腦；(2) LDSC 截距與比值；(3) 報告領頭變異中已在任一輸入 GWAS 顯著的比例。
- 敏感度：(1) OSA 改 INCLAVO 定義；(2) BMI 改 UKB ukb-b-19953；(3) HapMap3 與全部 SNP 兩版（全部 SNP 版視運算時間決定）。

## 八、預寫的可能失敗陳述

1. 停損觸發：OSA 遺傳力不足或與 BMI 幾乎等同，本題退場，不改題硬做。
2. P1 FAIL：OSA 的 SNP 遺傳力大多經由 BMI，參考文章的共同因子本質上是肥胖因子；這本身可寫成方法學註記。
3. P3 FAIL 但 P1 PASS：非肥胖成分存在但高度多基因，缺乏大效應位點。
4. 與 Campos 2023 方向不一致：不得宣稱複製。
5. 痛風虛無對照也出現相同富集：富集結果為 subtraction 流程產物，不得解讀。

## 九、主要分析以外不做的事

不做 PRS 預測效度宣稱（無個體層資料）；MR 結果只談方向，不換算臨床效果大小。
