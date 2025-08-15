# monocyte-redox-scRNA

Data and outputs for the manuscript on redox-associated programs in human peripheral blood mononuclear cells (PBMCs).

- **Dataset:** 10x Genomics **“PBMCs from Citrate-Treated Cell Preparation Tubes (3′ v3.1 Chemistry)”**  
  File used: `3p_Citrate_CPT_molecule_info.h5` (gene-expression only; no ADT).  
- **Cells used:** **3,940 pass-filter cells**  
- **HVGs:** top **2,000** highly variable genes  
- **PCA:** **25 components**, cumulative explained variance **0.668**  
- **Normalization:** library-size scaling to **10,000 counts/cell** → **natural log** transform (`log1p`)  
- **PCA method:** sparse SVD via **`irlba`** on the log-normalized matrix

---

## Repository layout

