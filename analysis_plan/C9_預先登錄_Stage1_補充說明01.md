# C9 預先登錄 Stage 1 補充說明 01（操作定義，在多變量 GWAS 執行前補登）

**登錄時間** 見同夾 `SHA256_登錄時間戳.txt` 第二列。登錄時已有的結果只有停損測試與 P1（皆與本補充說明的定義無關）；P3 所需的多變量 GWAS 尚未執行。

1. **位點的操作定義（P3、Table 3 用）**：本機無 plink LD clumping，改用距離法。
   全基因體顯著 SNP (P < 5 × 10⁻⁸) 依 P 值由小到大，取最小者為領頭 SNP，其前後 500 kb 內的 SNP 併入同一位點並移除，重複至無剩餘。
   座標一律用 1000G Phase3 bim 的 hg19 位置（以 rsID 對應），兩個 GWAS 用同一套座標。
2. **P3 的比對規則**：FinnGen OSA 單一 GWAS 的每個位點，若其領頭 SNP 前後 500 kb 內任一 SNP 在 NonBMI-OSA GWAS 達 P < 5 × 10⁻⁸，記為「保留」。P3 = 保留位點數 / FinnGen OSA 位點數。
3. **SNP 範圍**：主分析為 HapMap3 SNP（與 LDSC 相同），排除 MHC (hg19 chr6 25 到 35 Mb)。
4. **多變量 GWAS 模型**：Demange 2021 GWAS-by-subtraction 設定，SNP 同時預測 F_BMI 與 F_NonBMI；報告 F_NonBMI ~ SNP 的效應與 P。
5. **FUMA 定位**（使用者上傳）屬補充分析，不改變 P3 的判讀。
