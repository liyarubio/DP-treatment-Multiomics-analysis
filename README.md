# DP-treatment-Multi-omics-analysis
Code for "Decitabine-induced T cell remodeling facilitates a high antitumor response to PD-1 blockade therapy by promoting the expansion and effector function of CD8+ progenitor exhausted T cells"

## [ATAC-seq_data_processing](https://github.com/liyarubio/DP-treatment-Multi-omics-analysis/blob/main/ATAC-seq_data_processing)
## [WGBS_data_processing](https://github.com/liyarubio/DP-treatment-Multi-omics-analysis/blob/main/WGBS_data_processing)
## [In_vivo_assay_seurat.R](https://github.com/liyarubio/DP-treatment-Multi-omics-analysis/blob/main/In_vivo_assay_seurat.R)
## [ACT_assay_seurat.R](https://github.com/liyarubio/DP-treatment-Multi-omics-analysis/blob/main/ACT_assay_seurat.R)
## [Velocyto](https://github.com/liyarubio/DP-treatment-Multi-omics-analysis/blob/main/Velocyto)
## [pySCENIC](https://github.com/liyarubio/DP-treatment-Multi-omics-analysis/blob/main/pySCENIC)
## [Analysis and visualization.R](https://github.com/liyarubio/DP-treatment-Multi-omics-analysis/blob/main/Analysis%20and%20visualization.R)

# Analysis and visualization.R includes:
The code is ordered by Figure & Supplementary figure

### 0 PATH
### 0 Source data
### 0 color
### 1 Fig1G Fig.S3F Fig.S3H
#### 1.1 DEG
#### 1.2 volcano plot
### 2 Fig1G FigS3G
#### 2.1 GO
### 3 Fig1H
#### 3.1 expression heatmap
### 4 Fig S3A-l
#### 4.1 featureplot
#### 4.2 phase
#### 4.3 celltype
#### 4.4 celltype number & ratio
#### 4.5 Venn
### 5 FigS3D
#### 5.1 add TCR information
### 6 Fig S1E-H
#### 6.1 DMR: Differential Methylated Region
#### 6.2 DMR annotation
#### 6.3 KEGG
#### 6.4 DMR signal heatmap
##### 6.4.1 type transfer
##### 6.4.2 DMR ranger
#### 6.5 Boxplot+vlnplot Methy level
##### 6.5.1 TSS position
##### 6.5.2 DAC CpG in promoter
##### 6.5.3 Ctrl CpG in promoter
##### 6.5.4 Promoter CpG boxplot
### 7 Fig4 L-M
#### 7.1 tsne
#### 7.2 dotplot
#### 7.3 ratio
### 8 Fig S7 CD8+ T cells
#### 8.1 featureplot
#### 8.2 QC vInplot
#### 8.3 phase
### 9 Fig 4O Fig S8B
#### 9.1 monocle trjectory
#### 9.2 monocle gene exp
#### 9.3 add annotation to monocle gene exp
### 10 FigS8C-D
#### 10.1 destiny diffusion map
#### 10.2 get embeddings from destiny
#### 10.3 velocyto
### 11 Fig5A FigS9A
#### 11.1 immunarch
#### 11.2 diversity
#### 11.3 top proportion
### 12 Fig5B-D
#### 12.1 CLone size
#### 12.2 clone size TSNE
#### 12.3 size ratio by cluster
#### 12.4 highly expanded (size > 10)
### 13 Fig 5F
#### 13.1 top 10 clone by sample
#### 13 Fig 5E
#### 13.1 heatmap
### 14 FigS9H S9J
#### 14.1 GO
#### 14.2 CO Top50
### 15 Fig6A-B
#### 15.1 DEG
#### 15.2 volcano
#### 15.3 GO of all upregulated DEGs
### 16 Fig6C
#### 16.1 heatmap
#### 16.2 DEG Module
### 17 Fig6D FigS10A
#### 17.1 UP genes module GO
#### 17.2 UP genes module KEGG
### 18 FigS10B
### 19 Fig6E FigSB-E
#### 19.1 DP vs PD1 volcano
#### 19.2 PD1 vs Ctrl
#### 19.3 ATAC-seg heatmap correlation
#### 19.4 merge all sample peaks
#### 19.5 Venn for peaks
### 20 Fig6F FigS10-G
#### 20.1 SCENIC
#### 20.2 TF heatmap
#### 20.3 regulatory network
### 21 Fig7A FigS11A-B
#### 21.1 JunD exp in TSNE
#### 21.2 sclmpute
#### 21.3 JunD exp in public data
#### 22 Fig7B
#### 22.1 activation score
### 23 Fig7C-D
#### 23.1 GSEA
#### 23.2 GSEA Barplot
### 24 FigS12A-B
#### 24.1 GSEA
#### 24.2 Barplot
### 25 Fig7E & FigS12E
#### 25.1 Jund exp level class
#### 25.2 Jund & activation
#### 25.3 Jund level class in ACT assay
#### 25.4 Jund & activation score in ACT assay
