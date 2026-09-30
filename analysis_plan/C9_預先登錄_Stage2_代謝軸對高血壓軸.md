# C9 預先登錄 Stage 2：代謝軸對高血壓軸（A1 到 A5）

**登錄時間** 見同夾 `SHA256_登錄時間戳.txt`。本檔在 A1 到 A5 任何分析執行、任何新資料下載之前寫成。
**使用者裁示** 2026-09-26：以代謝軸對高血壓軸為主軸，A1 到 A5 全部進行；B1 (All of Us) 未啟動。
**登錄時已知的結果**（須揭露）：v2 的 P1 0.750、P2 0.691；次要描述 rg(OSA,HTN) 0.452 → rg(F_NonBMI,HTN) 0.288，rg(OSA,T2D) 0.385 → 0.119。
A1 正是對這個已知的描述性差異做正式檢定，所以 A1 屬「確認一個已看到的型態」，不是盲測；判讀時照實揭露。A2 到 A5 的資料在登錄時尚未下載。
**固定不變** Stage 1 的資料、SNP 集合、LDSC 設定（sep_weights = TRUE、MHC 25 到 35 Mb 排除）、Cholesky 模型。

## A1 主要：兩條軸的正式檢定

- **模型**：F_BMI =~ BMI + OSA；F_NonBMI =~ OSA（同 Stage 1）。T2D 與 HTN 各設為單指標潛在因子 F_T2D、F_HTN，變異數固定 1、殘差固定 0，
  因此所有因子間共變數即為遺傳相關。F_T2D、F_HTN 與 F_BMI、F_NonBMI 相關皆自由估計；F_T2D ~~ F_HTN 自由。
- **資料**：T2D = Xue 2018 (GCST006867)；HTN = FinnGen R12 I9_HYPTENSESS。
- **檢定**：限制 cor(F_NonBMI, F_T2D) = cor(F_NonBMI, F_HTN) 的模型，對未限制模型做卡方差異檢定 (df 1)。
- **判讀**
  - **PASS（兩條軸成立）**：cor(F_NonBMI, F_HTN) > cor(F_NonBMI, F_T2D) 且卡方差異 P < 0.05。
  - **FAIL**：P ≥ 0.05 或方向相反。FAIL 時全文不得使用「兩條軸」或「分離的路徑」，改寫為「程度不同的經肥胖比例」並揭露檢定不顯著。
- **A1b（次要，同規則，不影響主軸判讀）**：CAD (Aragam 2022) 對 T2D。

## A2 連續血壓複製

- **資料**：收縮壓 GCST006624、舒張壓 GCST006630（Evangelou 2018，PMID 30224653，757,601 人，連續性狀）。
- **做法**：A1 模型中以收縮壓取代 HTN；舒張壓為次要。
- **判讀**：**PASS** = cor(F_NonBMI, SBP) > cor(F_NonBMI, T2D) 且卡方差異 P < 0.05；**FAIL** 則「高血壓軸」只能以診斷型 HTN 描述，並揭露連續血壓未複製。

## A3 BMI 與 OSA 雙向孟德爾隨機化（支撐 Cholesky 的 BMI 在前）

- **工具變數**：P < 5 × 10⁻⁸，LD clumping r² < 0.001、視窗 10,000 kb，1000G Phase3 EUR 參考（本機 plink 對 g1000_eur）。
- **方向**：BMI (GIANT) → OSA (FinnGen R12)；OSA (FinnGen R12) → BMI (GIANT)。兩資料集樣本不重疊。
- **方法**：IVW 為主；MR-Egger、加權中位數、加權眾數；Steiger 方向過濾；Cochran Q；Egger 截距。
- **判讀**
  - **Cholesky 順序受支持**：BMI → OSA 的 IVW 與加權中位數皆為正且 P < 0.05。
  - **反向效應旗標**：若 OSA → BMI 的 IVW 與加權中位數皆 P < 0.05 且同向，則 F_BMI 可能部分包含 OSA 對 BMI 的反向效應，P1 因此可能低估不經 BMI 的比例；須在限制中量化，且不得宣稱 BMI 單向在前。
- 只作方向性與因果順序的支撐證據，不換算臨床效果大小。

## A4 偵測與社經偏差檢查

- **資料**：家戶收入 GCST009523（Hill 2019，PMID 31844048，286,301 人，非 MTAG）為主；認知表現 GCST006572（Lee 2018）為次要。
- **做法**：以 A1 同一 Cholesky，外部相關變項為收入。
- **判讀**：|cor(F_NonBMI, 收入)| ≥ 0.20 且 95% CI 不含 0 → **旗標**：非 BMI 成分可能含就醫或社經偏差，須在結果與限制明寫，且 A1 的解讀須附此旗標。未達者寫「未見明顯社經相關」。

## A5 機轉面板（描述性，無單項門檻）

- **性狀**：打鼾 GCST009761（Campos 2020）、睡眠時長 GCST007561（Dashti 2019）、晝夜型 GCST007576（Jones 2019）、
  FEV1/FVC GCST007431（Shrine 2019）、憂鬱 GCST007342（Howard 2019）、失眠 FinnGen F5_INSOMNIA、氣喘 FinnGen J10_ASTHMA_EXMORE。
- **報告**：每一性狀的 rg(OSA,X)、cor(F_BMI,X)、cor(F_NonBMI,X)；多重比較以 Bonferroni（7 性狀，P < 0.0071）標註。
- **H2 支持的操作定義**：至少 2 個非肥胖路徑性狀（打鼾、失眠、FEV1/FVC、憂鬱、晝夜型）同時滿足
  cor(F_NonBMI,X) 過 Bonferroni 且 |cor(F_NonBMI,X)| > |cor(F_BMI,X)|。未達者寫「面板未指向單一非肥胖路徑」。

## 預寫的可能失敗陳述

1. A1 FAIL：代謝軸與高血壓軸的差異不顯著，主軸改為程度描述，標題不得使用 "separates"。
2. A2 FAIL：高血壓軸限於診斷型 HTN，可能反映診斷或治療行為，不能推到血壓本身。
3. A3 反向效應成立：Cholesky 的 BMI 在前假設受質疑，P1 須加註上下界解讀。
4. A4 旗標：非 BMI 成分與社經相關，就醫偏差不能排除，結論降為假設生成。
5. A5 未指向路徑：機轉段落只寫「待個體層與功能研究」。
