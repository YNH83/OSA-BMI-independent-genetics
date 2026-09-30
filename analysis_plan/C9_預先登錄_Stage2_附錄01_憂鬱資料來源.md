# C9 Stage 2 附錄 01：A5 憂鬱資料來源更換（A5 任何資料下載與分析之前）

- Stage 2 登錄的憂鬱資料為 GCST007342（Howard 2019，PMID 30718901）。
- 查證 GWAS Catalog 目錄後發現該研究只公開前 10,000 個 SNP（檔名 PGC_UKB_23andMe_depression_10000.txt，因含 23andMe），無法做 LDSC 遺傳相關。
- 改用 FinnGen R12 `F5_DEPRESSIO`（Depression，59,333 例 / 434,831 對照），與 OSA 同庫，樣本重疊由 LDSC 截距處理。
- 其餘 A5 性狀、判讀規則與 H2 操作定義不變。
