###############################################################################
#  Skin Marker Reference  →  scMRMA format (gene | Level1 | Level2 | Level3 | Level4)
#  ───────────────────────────────────────────────────────────────────────────
#  1.  編輯  marker_sets_KGD_ChatGPT      ← 想增刪任何細胞類型或基因，只改這裡
#  2.  編輯  level_info       ← 新增細胞類型時，在這裡定義 4 個階層名稱
#  3.  執行 build_scMRMA()   ← 自動產生 scMRMA 參考表
###############################################################################

## ---------------------------------------------------------------------------
##  1)  marker_sets_KGD_ChatGPT  (每一個 cell type 對應一組 marker gene 向量)
## ---------------------------------------------------------------------------
marker_sets_KGD_ChatGPT <- list(
  "Epithelial cells"        = c("EPCAM","CDH1","MUC1","KRT5","KRT7","KRT8",
                                "KRT14","KRT18","KRT19","DSP","SDC1","CLDN1",
                                "CLDN4","CD24","ITGB4"),
  "Endothelial cells"       = c("PECAM1","CDH5","ENG","CD34","ICAM2","TIE1",
                                "TEK","CLDN5","ESAM","FLT1","KDR","EMCN",
                                "ERG","NOS3","ID3"),
  "Granulocytes cells"      = c("ITGAM","ITGAX","ITGB2","ANPEP","CD33","FUT4",
                                "CEACAM8","FCGR2A","CD63","CR1","C5AR1","CEBPE",
                                "PTPRC","CSF2RB","CSF3R"),
  "Keratinocytes"           = c("KRT5","KRT14","KRT1","KRT10","KRT6A","KRT16",
                                "KRT17","KRT15","KRT19","IVL","LOR","FLG",
                                "DSP","AQP3","SFN"),
  "Sweat gland cells"       = c("MUCL1","PIP","AQP5","KRT7","KRT8","KRT18",
                                "CA2","DCD","LYZ","SCGB1D2","CRNN","KRT19",
                                "KRT14","ACTA2","CALML5"),
  "Sebaceous gland cells"   = c("MUC1","KRT7","IHH","HRH1","SOX9","PPARG",
                                "MC5R","DGAT2","FABP5","NPC2","CYP11A1","ACSL5",
                                "KRT5","KRT14","MUC4"),
  "Smooth muscle cells"     = c("ACTA2","MYH11","TAGLN","CNN1","DES","LMOD1",
                                "SMTN","MYL9","MYLK","MCAM","ITGA8","GJA4",
                                "KCNMB1","RGS5","PDGFRB"),
  "Arrector pili muscle cells" = c("ACTG2","ITGA8","TAGLN","MYH11","CNN1",
                                   "DES","NOTCH3","JUNB","MCAM","PDGFRB",
                                   "LBH","MEF2C","FOXO1","PRPH","NGFR"),
  "Pericytes"               = c("PDGFRB","MCAM","CSPG4","RGS5","KCNJ8","ABCC9",
                                "ACTA2","DES","NOTCH3","CD248","ANPEP","ZIC1",
                                "VIM","ANGPT1","P2RY14"),
  "Fibroblasts"             = c("COL1A1","COL1A2","PDGFRA","DCN","LUM",
                                "COL3A1","COL5A1","COL6A3","POSTN","FAP",
                                "THY1","PDPN","CD34","SFRP2","ASPN"),
  "Vascular endothelial cells" = c("PECAM1","VWF","CD34","KDR","CLDN5","CDH5",
                                   "ENG","ESM1","ICAM2","PLVAP","NOS3","TIE1",
                                   "TEK","VCAM1","ACKR1"),
  "Lymphatic endothelial cells"= c("PROX1","LYVE1","PDPN","FLT4","NRP2","CCL21",
                                   "ACKR4","VCAM1","ITGA9","EMILIN1","MADCAM1",
                                   "RELN","CCL21","PECAM1","ENG"),
  "T cells"                 = c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B",
                                "CCR7","IL7R","GZMK","GNLY","NKG7","IFNG",
                                "CXCR3","ITGAE","FOXP3"),
  "NK cells"                = c("NCAM1","FCGR3A","NKG7","GNLY","PRF1","GZMB",
                                "GZMH","KLRD1","KLRF1","NCR1","XCL1","FGFBP2",
                                "IL2RB","KLRC1","KLRC2"),
  "B cells"                 = c("MS4A1","CD19","CD79A","CD79B","IGHD","IGHM",
                                "TCL1A","CR2","FCRL5","PAX5","SELL","CCR7",
                                "SEC11C","MME","BANK1"),
  "Plasma cells"            = c("MZB1","JCHAIN","XBP1","SDC1","SLAMF7",
                                "TNFRSF17","PRDM1","IGHG1","IGHA1","IGKC",
                                "FKBP11","DERL3","SRPX2","HSP90B1","TNFRSF13B"),
  "Monocytes"               = c("CD14","FCGR3A","S100A8","S100A9","LYZ","VCAN",
                                "CCR2","IL1B","CSF1R","ITGAM","ITGAL","HLA-DRA",
                                "FCGR2A","NR4A1","APOBEC3A"),
  "Macrophages"             = c("CD68","CD163","MRC1","MARCO","FCGR3A","HLA-DRA",
                                "APOE","C1QA","C1QB","C1QC","IL10RA","CCR2",
                                "FCER1G","NLRP3","SPI1"),
  "Dendritic cells"         = c("HLA-DRA","CD1C","CLEC9A","XCR1","CCR7","ITGAX",
                                "ITGAM","CD207","CD1A","IRF7","IL3RA","LAMP3",
                                "SIGLEC6","FCER1A","CCR6"),
  "Neutrophils"             = c("S100A8","S100A9","MPO","ELANE","PRTN3","CEACAM8",
                                "FCGR3B","CSF3R","CXCR2","ITGAM","FUT4","OLFM4",
                                "RETN","TREM1","IL1R2"),
  "Eosinophils"             = c("RNASE2","RNASE3","PRG2","EPX","IL5RA","CCR3",
                                "SIGLEC8","HPGDS","GATA1","CCR1","CEACAM8",
                                "ITGAX","ITGAM","ANPEP","PTGDR2"),
  "Basophils"               = c("FCER1A","IL3RA","ENPP3","CCR3","CD9","CD36",
                                "CD38","CD40LG","HDC","CPA3","HPGDS","SIGLEC8",
                                "KIT-","HLADR-","MS4A2"),
  "Mast cells"              = c("KIT","TPSAB1","TPSB2","CPA3","HPGDS","HDC",
                                "FCER1A","MS4A2","VEGFA","IL1RL1","VWA5A",
                                "SIGLEC6","MRGPRX2","TNF","CXCL8"),
  "Melanocytes"             = c("TYR","DCT","TYRP1","PMEL","MLANA","MITF",
                                "SOX10","GPR143","RAB27A","OCA2","SLC45A2",
                                "KIT","MC1R","ATP7A","SLC24A5"),
  "Neuronal cells"          = c("UCHL1","TUBB3","PRPH","SCN9A","TRPV1","CALCA",
                                "NTRK1","NGFR","MAP2","ELAVL4","STMN2","GAP43",
                                "CHRNA3","SUBP","NEFH"),
  "Schwann cells"           = c("MPZ","MBP","PMP22","S100B","GFAP","NCAM1",
                                "NGFR","ERBB3","SOX10","CCN3","PTN","L1CAM",
                                "PRX","CTSB","CDH19"),
  "Adipocytes"              = c("ADIPOQ","PLIN1","FABP4","LEP","GPD1","CIDEC",
                                "LPL","LIPE","CIDEA","PPARG","AKT2","SLC2A4",
                                "CEBPA","DNASE1L3","FGF21"),
  
  "Basal keratinocyte" = c(
    "KRT5","KRT14","KRT15","KRT17","COL17A1","ITGA3","ITGB4",
    "LAMB3","LAMC2","TP63","CDH3","CXCL14","IRF6","S100A14","EGFR"
  ),
  
  "Spinous keratinocyte" = c(
    "KRT1","KRT10","DSG1","IVL","PKP1","PKP3","DMKN","CALML5",
    "SBSN","KRT2","KRT16","TRIM29","DSC3","CSTA","SFN"
  ),
  
  "Granular keratinocyte" = c(
    "LOR","FLG","TGM1","TGM3","SPINK5","CASP14",
    "CRNN","IVL","CDSN","SPRR2A","CLDN1","CST6","PPL","EVPL","KRTDAP"
  ),
  
  "Mitotic keratinocyte" = c(
    "MKI67","TOP2A","AURKB","UBE2C","BIRC5","CCNB1","CCNB2",
    "CDC20","MCM2","MCM4","PCNA","CENPF","CDC6","TYMS","GMNN"
  ),
  
  "Secretory papillary fibroblast" = c(
    "APCDD1","WIF1","ID1","PTGDS","COL18A1","AXIN2","COLEC12",
    "COL6A1","ENPP2","ITM2A","DIO2","HSPB3","PDPN","COL6A3","CD34"
  ),
  
  "Secretory reticular fibroblast" = c(
    "WISP2","SLPI","TSPAN8","SFRP2","ACTA2","CNN1","COL11A1","FMO1",
    "DPP4","MFAP5","ELN","LOX","MMP2","MMP14","PDGFRA"
  ),
  
  "Pro-inflammatory fibroblast" = c(
    "CCL19","CXCL2","CXCL3","CXCL1","CXCL8","CXCL10",
    "IL6","IL11","MMP1","MMP3","CCL2","CXCL5","IL32","PTGS2","TNF"
  ),
  
  "Mesenchymal fibroblast" = c(
    "ASPN","POSTN","GPC3","TNN","SFRP1","ADAM12","COL1A1","COL3A1",
    "FN1","DCN","PDGFRA","THY1","PDGFRB","LRRC15","TNC"
  )
)

## ---------------------------------------------------------------------------
##  2)  level_info  (為每個 cell-type 指定 Level1~Level4 名稱)
##      ➜ 若新增 / 刪除 cell-type，僅需在此 list 增減對應項
## ---------------------------------------------------------------------------
level_info <- list(
  #── Epithelial 系
  "Epithelial cells"         = c("Epithelial cell",
                                 "Non-keratinizing epithelial cell",
                                 "Skin epithelial cell",
                                 "Epithelial cells"),
  "Keratinocytes"            = c("Epithelial cell",
                                 "Keratinizing epithelial cell",
                                 "Skin epidermal keratinocyte",
                                 "Keratinocytes"),
  "Sweat gland cells"        = c("Epithelial cell",
                                 "Secretory epithelial cell",
                                 "Skin eccrine sweat gland cell",
                                 "Sweat gland cells"),
  "Sebaceous gland cells"    = c("Epithelial cell",
                                 "Secretory epithelial cell",
                                 "Skin sebaceous gland cell",
                                 "Sebaceous gland cells"),
  
  #── Endothelial 系
  "Endothelial cells"        = c("Endothelial cell",
                                 "Endothelial cell",
                                 "Skin endothelial cell",
                                 "Endothelial cells"),
  "Vascular endothelial cells" = c("Endothelial cell",
                                   "Blood endothelial cell",
                                   "Skin vascular endothelial cell",
                                   "Vascular endothelial cells"),
  "Lymphatic endothelial cells"= c("Endothelial cell",
                                   "Lymphatic endothelial cell",
                                   "Skin lymphatic endothelial cell",
                                   "Lymphatic endothelial cells"),
  
  #── Mesenchymal 系
  "Smooth muscle cells"      = c("Mesenchymal cell",
                                 "Smooth muscle cell",
                                 "Skin vascular smooth muscle cell",
                                 "Smooth muscle cells"),
  "Arrector pili muscle cells" = c("Mesenchymal cell",
                                   "Smooth muscle cell",
                                   "Skin arrector pili muscle cell",
                                   "Arrector pili muscle cells"),
  "Pericytes"                = c("Mesenchymal cell",
                                 "Perivascular mural cell",
                                 "Skin pericyte",
                                 "Pericytes"),
  "Fibroblasts"              = c("Mesenchymal cell",
                                 "Fibroblast",
                                 "Skin dermal fibroblast",
                                 "Fibroblasts"),
  "Adipocytes"               = c("Mesenchymal cell",
                                 "Adipocyte",
                                 "Subcutaneous adipocyte",
                                 "Adipocytes"),
  
  #── Immune – Granulocyte & Myeloid
  "Granulocytes cells"       = c("Immune cell",
                                 "Myeloid granulocyte",
                                 "Granulocyte (undefined)",
                                 "Granulocytes"),
  "Neutrophils"              = c("Immune cell",
                                 "Myeloid granulocyte",
                                 "Neutrophil",
                                 "Neutrophils"),
  "Eosinophils"              = c("Immune cell",
                                 "Myeloid granulocyte",
                                 "Eosinophil",
                                 "Eosinophils"),
  "Basophils"                = c("Immune cell",
                                 "Myeloid granulocyte",
                                 "Basophil",
                                 "Basophils"),
  
  #── Immune – Mononuclear phagocyte
  "Monocytes"                = c("Immune cell",
                                 "Myeloid mononuclear phagocyte",
                                 "Monocyte",
                                 "Monocytes"),
  "Macrophages"              = c("Immune cell",
                                 "Myeloid mononuclear phagocyte",
                                 "Macrophage",
                                 "Macrophages"),
  
  #── Immune – Others
  "Dendritic cells"          = c("Immune cell",
                                 "Dendritic cell",
                                 "Skin dendritic cell",
                                 "Dendritic cells"),
  "T cells"                  = c("Immune cell",
                                 "Lymphocyte",
                                 "T cell",
                                 "T cells"),
  "B cells"                  = c("Immune cell",
                                 "Lymphocyte",
                                 "B cell",
                                 "B cells"),
  "Plasma cells"             = c("Immune cell",
                                 "Lymphocyte",
                                 "Plasma cell",
                                 "Plasma cells"),
  "NK cells"                 = c("Immune cell",
                                 "Lymphocyte",
                                 "Natural killer cell",
                                 "NK cells"),
  "Mast cells"               = c("Immune cell",
                                 "Mast cell / Basophil lineage",
                                 "Mast cell",
                                 "Mast cells"),
  
  #── Neural / Neural-crest
  "Melanocytes"              = c("Neural-crest-derived cell",
                                 "Melanocyte",
                                 "Skin melanocyte",
                                 "Melanocytes"),
  "Neuronal cells"           = c("Neural cell",
                                 "Peripheral sensory neuron",
                                 "Cutaneous neuron",
                                 "Neuronal cells"),
  "Schwann cells"            = c("Neural-crest-derived cell",
                                 "Peripheral glial cell",
                                 "Schwann cell",
                                 "Schwann cells"),
  
  #──
  "Basal keratinocyte" = c(
    "Epithelial cell",
    "Keratinizing epithelial cell",
    "Basal keratinocyte",
    "Basal keratinocyte"
  ),
  
  "Spinous keratinocyte" = c(
    "Epithelial cell",
    "Keratinizing epithelial cell",
    "Spinous keratinocyte",
    "Spinous keratinocyte"
  ),
  
  "Granular keratinocyte" = c(
    "Epithelial cell",
    "Keratinizing epithelial cell",
    "Granular keratinocyte",
    "Granular keratinocyte"
  ),
  
  "Mitotic keratinocyte" = c(
    "Epithelial cell",
    "Keratinizing epithelial cell",
    "Mitotic keratinocyte",
    "Mitotic keratinocyte"
  ),
  
  "Secretory papillary fibroblast" = c(
    "Mesenchymal cell",
    "Fibroblast",
    "Secretory papillary fibroblast",
    "Secretory papillary fibroblast"
  ),
  
  "Secretory reticular fibroblast" = c(
    "Mesenchymal cell",
    "Fibroblast",
    "Secretory reticular fibroblast",
    "Secretory reticular fibroblast"
  ),
  
  "Pro-inflammatory fibroblast" = c(
    "Mesenchymal cell",
    "Fibroblast",
    "Pro-inflammatory fibroblast",
    "Pro-inflammatory fibroblast"
  ),
  
  "Mesenchymal fibroblast" = c(
    "Mesenchymal cell",
    "Fibroblast",
    "Mesenchymal fibroblast",
    "Mesenchymal fibroblast"
  )
  
)

## ---------------------------------------------------------------------------
##  3)  通用函式：marker_sets → scMRMA 資料框
## ---------------------------------------------------------------------------
build_scMRMA <- function(marker_sets, level_info) {
  
  # 3-1. 檢查：有沒有 cell-type 沒在 level_info 定義
  missing <- setdiff(names(marker_sets), names(level_info))
  if (length(missing)) {
    warning("【警告】以下 cell type 尚無階層定義，已跳過：\n",
            paste(missing, collapse = ", "))
    marker_sets <- marker_sets[setdiff(names(marker_sets), missing)]
  }
  
  # 3-2. 組裝為長格式
  do.call(
    rbind,
    lapply(names(marker_sets), function(ct) {
      genes <- marker_sets[[ct]]
      lv    <- level_info[[ct]]  # 取得對應 4 個階層名稱
      
      data.frame(
        gene   = genes,
        Level1 = lv[1],
        Level2 = lv[2],
        Level3 = lv[3],
        Level4 = lv[4],
        stringsAsFactors = FALSE
      )
    })
  )
}

## ---------------------------------------------------------------------------
##  4)  產生 scMRMA 參考表
## ---------------------------------------------------------------------------
marker_df_scMRMA_KGD_ChatGPT <- build_scMRMA(marker_sets_KGD_ChatGPT, level_info)

## ---------------------------------------------------------------------------
##  5)  選擇性：匯出 CSV
## ---------------------------------------------------------------------------
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 10) # Generate a unique time-based ID
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))

Set_note <- "scMRMA_Test"
Name_Export <- paste0(Name_FileID)
Name_ExportFolder_KGD_CTAnnot <- paste0("Export_",Set_note)
# Create export folder if it does not exist
if (!dir.exists(Name_ExportFolder_KGD_CTAnnot)){dir.create(Name_ExportFolder_KGD_CTAnnot)}

write.csv(marker_df_scMRMA_KGD_ChatGPT,
          file = paste0(Name_ExportFolder_KGD_CTAnnot, "/", Name_Export, "Skin_Marker_scMRMA_reference.csv"),
          row.names = FALSE)

## (檢視前幾列)
head(marker_df_scMRMA_KGD_ChatGPT, 10)




################################################################################
# Set_Run_DotPlots <- TRUE

if (isTRUE(get0("Set_Run_DotPlots", ifnotfound = FALSE, inherits = FALSE))) {
  #### Visualization: DotPlots ####
  marker_sets <- marker_sets_KGD_ChatGPT 
  
  Set_Idents <- "seurat_clusters"
  Set_Marker_DotPlot_ncol <- 4
  Set_Marker_DotPlot_width <- 30
  Set_Marker_DotPlot_height <- 60
  
  source("RUNPlot_Marker_DotPlots.R")
}




