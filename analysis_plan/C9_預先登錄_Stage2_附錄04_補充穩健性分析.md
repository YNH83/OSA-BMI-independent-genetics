# C9 預先登錄 Stage 2 附錄 04：補充穩健性分析（執行前登錄）

**登錄時間** 見 `SHA256_登錄時間戳.txt`。本附錄寫於任何一項下列分析執行之前；門檻與判讀規則在此寫死，看到結果後不再更改。
**來源** `09_圖表設計參考/C9_可套用分析_統整建議.md` 第二節「建議」四項。

## 不執行的兩項與理由（解析結果，非事後取捨）

1. **放寬 F_BMI 與 F_NonBMI 的正交假設（相關 0.1、0.2）**：F_BMI 只有 BMI 一個指標，標準化負荷為 1，因此 F_BMI 等同 BMI 的遺傳成分。無論兩因子相關設為多少，OSA 遺傳變異中可由 BMI 遺傳解釋的比例恆為 rg(OSA, BMI)²，不經 BMI 的比例恆為 1 − rg²。放寬正交只會把同一個量在兩個因子之間重新分配名稱，不改變結論，也無法作為穩健性證據。此點與已報告的「P1 約等於 1 − rg²、痛風虛無比例更高、P1 本身不具專一性」一致，將寫入方法與限制。
2. **逐位點的模型比較（只經 BMI 對另有 F_NonBMI 路徑）**：在每個 SNP 自由度為 0 的 Cholesky 模型中，兩模型的差異就是 F_NonBMI 路徑上的 SNP 效應是否為 0，其檢定等同 userGWAS 已報告的 F_NonBMI SNP 的 Wald P。重做不提供新資訊。

## S1 奇偶染色體分半（P1 與 P2 的穩定度）

- **做法**：與主分析相同的三性狀 LDSC（OSA_FG、BMI_GIANT、T2D_XUE；sep_weights = TRUE；MHC 排除），以 GenomicSEM `ldsc(select = "ODD")` 與 `select = "EVEN"` 各估一次，套用主分析的 Cholesky 模型，算 P1 與 P2。
- **判讀**：兩半的 P1 都 ≥ 0.60 且兩半 95% CI 重疊 → P1 分半穩定；任一半 < 0.60 → 報告為「P1 依賴特定染色體區段」。P2 同理，門檻 0.50。只作穩健性描述，不改變 Stage 1 判定。

## S2 Genomic SEM 多元迴歸（打鼾是否獨立於憂鬱與失眠）

- **做法**：五性狀 LDSC（OSA_FG、BMI_GIANT、SNORING_CAMPOS、INSOMNIA_FG、DEPRESSION_FG）。主分析 Cholesky 結構不變，加上 F_NonBMI ~ SNORING_CAMPOS + INSOMNIA_FG + DEPRESSION_FG，三個外部性狀之間的遺傳共變異自由估計。報各自的標準化迴歸係數與 95% CI。
- **敏感度**：以 BMI 校正版打鼾（SNORING_CAMPOS_ADJBMI）取代。
- **判讀**：打鼾的標準化係數為正、95% CI 不含 0，且點估計大於失眠與憂鬱 → 寫「上呼吸道訊號不能由憂鬱與失眠重複解釋」；打鼾係數 CI 含 0 → 寫「三者共享，無法區分獨立貢獻」。失眠、憂鬱係數只描述。模型不收斂或出現負變異數時如實報告，不改模型。

## S3 位點共定位（coloc）

- **位點**：Table 4 的 5 個 F_NonBMI 位點（KCNK3、CDKN2B、METTL15、MSRB3、FTO），lead ± 500 kb。
- **配對**：F_NonBMI 對 BMI (GIANT)、F_NonBMI 對打鼾 (Campos 2020)、原始 OSA 對 BMI。
- **方法**：coloc.abf（coloc R 套件），預設先驗 p1 = p2 = 1e-4、p12 = 1e-5；連續性狀以 β 與 varβ 輸入，F_NonBMI 的 sdY = 1、N = N_hat；二元性狀 type = "cc"、以有效樣本數與病例比例輸入。只用兩檔共有的 HapMap3 SNP。
- **判讀（每一對）**：PP.H4 ≥ 0.80 → 共享因果變異；PP.H3 ≥ 0.80 → 不同因果變異；其餘 → 無法判定。
- **預測**：若相減有效，非 FTO 的四個位點與 BMI 的 PP.H4 應 < 0.50；FTO 與 BMI 的 PP.H4 預期高（與 P3 的既定限制一致，不得作為「不經肥胖」的證據）。與打鼾共享者只作支持性描述。
- **限制（預先寫明）**：摘要統計只含 HapMap3 SNP，SNP 密度較低時 coloc 傾向 H3 或無法判定；單一因果變異假設。

## 預寫的可能失敗陳述

1. S1 任一半 P1 < 0.60：P1 的 PASS 須附「分半不穩」說明。
2. S2 打鼾係數 CI 含 0：上呼吸道主張降為「與憂鬱、失眠共享的非肥胖成分」，不得寫上呼吸道特異。
3. S3 非 FTO 位點與 BMI 的 PP.H4 ≥ 0.80：該位點不得在正文列為不經 BMI 的位點，Figure 5 與 Table 4 加註。
