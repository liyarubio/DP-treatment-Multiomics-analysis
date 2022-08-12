# Analysis &  Visualization  
# Created by Yaru, Li 20220812 V1
# The code is ordered by Figure & Supplementary figure
# The processing data (csv, xlsx, rds ... ) of the code will be uploaded after the article is published

library(Seurat)
library(patchwork)
library(dplyr)
library(plyr)
library(ggplot2)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(monocle3)
library(data.table)
library(readxl)
library(data.table)
library(ggpubr)
library(rstatix)
library(grid)
library(pheatmap)
library(scales)
library(org.Mm.eg.db)
library(clusterProfiler)
library(Hmisc)
library(VennDiagram)
library(destiny)

####-----0 PATH------------------------------
# in vivo assay
setwd("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP")  

# ACT assay
setwd("/media/liyaru/LYR1/301project/3_clone")

####-----0 Source data------------------------
CD8 <- readRDS("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/RDS/new/CD8_C0-6.rds")
DimPlot(CD8)

merge_inter <- readRDS("/media/liyaru/LYR1/301project/3_clone/result/RDS/merge_inter.rds")
DimPlot(merge_inter)

####----0 color--------------------------------
color1 <- brewer.pal(8, "Set1")
color1 <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF")

color2 <- brewer.pal(8, "Set2")

rdwhbu <- colorRampPalette(c("navy", "white", "brown3"))
rdwhbu_re <- colorRampPalette(c("brown3", "white", "navy"))

color_4_sample <- c("black","#285291","#579125","#B9181E")
color_4_sample_1 <- c(rgb(23,23,23,160, maxColorValue = 255),
                      "#285291","#579125","#B9181E")

color3 <- c("#FEE8C8","#FDBB84","#FC8D59","#D7301F","#7F0000")

color4 <- c("#FDD49E","#EF6548","#7F0000")

color7 <- brewer.pal(7, "Set2")

####-----1 Fig1G Fig.S3F Fig.S3H----------------------------
####-----1.1 DEG-----------------------------
merge_inter@meta.data$cluster_sample <- paste0(merge_inter@meta.data$celltype,"_",merge_inter@meta.data$sample)
table(merge_inter$cluster_sample)
DefaultAssay(merge_inter) <- "RNA"

# DP vs p
Pro.markers <- FindMarkers(merge_inter,
                           group.by = "cluster_sample",
                           ident.1 = "Proliferating_DP",
                           ident.2 ="Proliferating_PD1",
                           #test.use = "MAST"
                           logfc.threshold=0
)
fwrite(Pro.markers,
       "./result/table/2_Pro_DEG.csv",
       row.names =T )

# P vs C
Pro.markers <- FindMarkers(merge_inter,
                           group.by = "cluster_sample",
                           ident.1 = "Proliferating_PD1",
                           ident.2 ="Proliferating_Con",
                           #test.use = "MAST"
                           logfc.threshold=0
)
fwrite(Pro.markers,
       "./result/table/2_Pro_DEG_P_C.csv",
       row.names =T )

# D vs C
Pro.markers <- FindMarkers(merge_inter,
                           group.by = "cluster_sample",
                           ident.1 = "Proliferating_DAC",
                           ident.2 ="Proliferating_Con",
                           #test.use = "MAST"
                           logfc.threshold=0
)
fwrite(Pro.markers,
       "./result/table/2_Pro_DEG_D_C.csv",
       row.names =T )

####-----1.2 volcano plot----------------
library(ggplot2)
library(ggrepel)

# from Supplimentary table S3
a <- fread("./result/table/2_Pro_DEG.csv")

# a <- fread("./result/table/2_Pro_DEG_P_C.csv")  # Fig.S3F & Fig.S3H
# a <- fread("./result/table/2_Pro_DEG_D_C.csv")
# a <- fread("./result/table/2_Non-Pro_DEG.csv")

Dat<- as.data.frame(a)

Dat$gene <- Dat$V1

fc = 0.2
fc1 = 0.5

Dat$threshold <- 'Other'
Dat[Dat$p_val_adj < 0.05 & Dat$avg_log2FC > fc,'threshold'] <- 'Up'
Dat[Dat$p_val_adj < 0.05 & Dat$avg_log2FC > fc1,'threshold'] <- 'Up-Top'
Dat[Dat$p_val_adj < 0.05 & Dat$avg_log2FC < -fc,'threshold'] <- 'Down'
Dat[Dat$p_val_adj < 0.05 & Dat$avg_log2FC < -fc1,'threshold'] <- 'Down-Top'
table(Dat$threshold)
Dat$threshold <- factor(Dat$threshold,levels = c('Up','Up-Top','Down','Down-Top','Other'))

ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c(  "#EEBBBB","#CD3333","#AAAAD4", "#000080","#808080"))+
  geom_text_repel(
    #data = Dat[(Dat$p_val_adj<0.05 & abs(Dat$avg_log2FC) > fc1)| Dat$V1 %in% c("Jund","Jun") ,],
    data = Dat[(Dat$p_val_adj<0.05 & abs(Dat$avg_log2FC) > fc1),],
    aes(label = gene),
    size = 5,
    segment.color = "black", show.legend = FALSE,
    max.overlaps=40)+#添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank(),#不显示图例标题
  )+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05



####------2 Fig1G  FigS3G-----------
####------2.1 GO -------------------
# from Supplimentary table S3
a <- fread("./result/table/2_Pro_DEG.csv")
#a <- fread("./result/table/2_Pro_DEG_P_C.csv")
#a <- fread("./result/table/2_Pro_DEG_D_C.csv")

a <- a[a$avg_log2FC > 0.2 & a$p_val_adj < 0.05,]
ego <- enrichGO(
  gene          = a$V1,
  keyType       = "SYMBOL",
  OrgDb         =  org.Mm.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  #readable      = TRUE
)

t <- ego@result
t <- t[as.numeric(t$p.adjust )<=0.05,]

#select terms
terms <-c ("T cell differentiation",
           "ribonucleoprotein complex biogenesis",
           "protein folding",
           "cellular response to chemical stress",
           "regulation of T cell activation",
           "nuclear transport",
           "regulation of immune effector process",
           "regulation of adaptive immune response",
           "regulation of DNA metabolic process",
           "protein polyubiquitination",
           "regulation of mRNA metabolic process",
           "T cell proliferation",
           "cell killing")

t <- t[t$Description %in% terms,]

ttt <- t$Description
barplot(ego,showCategory = ttt)
clusterProfiler::dotplot(ego,showCategory = ttt)+
  scale_color_gradientn(colors = c(rdwhbu_re(100)[1:45],rdwhbu_re(100)[55:100]))

####-------3 Fig1H---------------------
####-------3.1 expression heatmap-------
# get average expression
t <- AverageExpression(merge_inter,assays = "RNA",group.by = "cluster_sample")
t <- t[["RNA"]]
t <- as.data.frame(t)
t <- t[,c("Non-Proliferating_Con","Non-Proliferating_DAC","Non-Proliferating_PD1","Non-Proliferating_DP",
          "Proliferating_Con","Proliferating_DAC","Proliferating_PD1","Proliferating_DP")]

# get from table in averageExpression
a=as.data.frame(read_excel("./result/table/genes_for_heatmap.xlsx"))
rownames(a) <- a$...1
a <- a[,2:9]

pheatmap::pheatmap(a,scale = "row",
                   cluster_rows=F,cluster_cols = F,
                   color=rdwhbu(100))

####---------4 Fig S3A-I -----------------
####---------4.1 featureplot--------------
FeaturePlot(merge_inter,
            features = c("Cd3d","Cd3e",
                         "Cd8a","Cd8b1",
                         "Pdcd1","Ctlas","Tigit",
                         "Havcr2","Gzmb","Nr4a2","Lef1",
                         "Gzmk","Lef1","Tcf7",
                         "Ifng","Cx3cr1","Mki67"),
            
            ncol=4,
            pt.size = 1)& 
  scale_color_gradientn(colors = c(rdwhbu(100)[1:45],rdwhbu(100)[55:100]))

####---------4.2 phase--------------
g2m.genes<-c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2",
             "Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2",
             "Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1",
             "Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8",
             "Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2",
             "G2e3","Gas2l3","Cbx5","Cenpa")
s.genes<-c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl",
           "Prim1","Uhrf1","Mlf1ip","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76",
           "Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1",
           "Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8")

merge_inter <- CellCycleScoring(merge_inter, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(merge_inter,group.by = "Phase",
        cols = color1[3:5],
        pt.size = 0.6)

####---------4.3 celltype--------------
DimPlot(merge_inter,group.by = 'celltype',
        label = F,
        #pt.size = 1,
        label.size = 7,
        cols=color1[c(2,1)],
        pt.size = 0.6)

#####--------4.4 celltype number & ratio-----------------------
m <- merge_inter@meta.data
m$number <- 1

ggplot(m,aes(sample,number,fill=celltype))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 16)


m <- ddply(m,'sample',transform,percent = 1/sum(number)*100)
colnames(m)
m$percent

ggplot(m,aes(sample,percent,fill=celltype))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 20)+
  scale_fill_manual(values = color1[c(2,1)])

t<- table(merge_inter$orig.ident,merge_inter$celltype)
x <- matrix(t[3:4,1:2],nrow = 2,ncol = 2)
fisher.test(x)

###-------4.5 Venn--------------------------
a <- fread("./result/table/2_Pro_DEG.csv")
a <- a[a$avg_log2FC > 0.2 & a$p_val_adj < 0.05,]
pro <- a$V1

a <- fread("./result/table/2_Non-Pro_DEG.csv")
a <- a[a$avg_log2FC > 0.2 & a$p_val_adj < 0.05,]
no <- a$V1

t <- intersect(pro,no)

up <- c(list(pro),list(no))
names(up) <- c('Pro','NO-Pro')

p = venn.diagram(
  x = up,
  filename=NULL,
  fill= color1[1:2],
  width = 1000,
  height = 1000, 
)

grid.draw(p)

#####--------5 FigS3D ----------------------
###------ 5.1 add TCR information-----------------------
sup <- c("_1","_2","_3","_4")
j=1
for (i in c("CON-0830-TCR","D-0830-TCR","DP-0830-TCR","P-0830-TCR")){
  tcr <- read.csv(paste0("./cellranger/2_TCR/",i,"/outs/filtered_contig_annotations.csv"))
  tcr <- tcr[!duplicated(tcr$barcode), ]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  
  clono <- read.csv(paste0("./cellranger/2_TCR/",i,"/outs/clonotypes.csv"))
  
  tcr <- merge(tcr, clono)
  tcr$sample <- i
  tcr$barcode <- paste0(tcr$barcode,sup[j])
  j=j+1
  
  if (i=="CON-0830-TCR"){
    tcr_all <-tcr
  }else{
    tcr_all <- rbind(tcr_all,tcr)
  }
}
rownames(tcr_all) <- tcr_all$barcode

merge_inter <- AddMetaData(merge_inter,metadata = tcr_all)

merge_inter$sample_clonotype_id <- paste0(merge_inter$orig.ident,"_",merge_inter$clonotype_id)

merge_inter$clonotype <- merge_inter$clonotype_id 
merge_inter@meta.data[merge_inter@meta.data$clonotype_id %in%  c('clonotype1'),'clonotype'] <- 'Main'
merge_inter@meta.data[grepl('clonotype',merge_inter@meta.data$clonotype),'clonotype'] <- 'Other'

DimPlot(merge_inter,group.by = "clonotype",
        cols=color1[1:2],
        pt.size = 0.6)

t <- merge_inter@meta.data
fwrite(t,"../Table/ACT_assay_metadata.csv",row.names = T)

####--------6 Fig S1E-H---------------
####--------6.1 DMR: Differential Methylated Region----------
library("DSS")
library("bsseq") 
read.faster <- function(file, header = FALSE, sep = "\t", showProgress = TRUE, skip = 0){
  
  suppressPackageStartupMessages(library("data.table"))
  Tx = data.frame(fread(file = file, header = header, sep = sep, showProgress = showProgress, skip = skip))
  return(Tx)
}

read.bedgraph <- function(file, showProgress = FALSE){
  
  Tx = read.faster(file = file, header = FALSE, sep = "\t", skip = 1, showProgress = showProgress)
  return(Tx)
}

SRX = c("Ctrl","DAC")

#WGBS CpG bedgraph
files=c("/media/liyaru/LYR1/301project/1_MOUSE_T_Chi_DAC/methy_result/4_MethylDackel/CON_CpG.bedGraph",
        "/media/liyaru/LYR1/301project/1_MOUSE_T_Chi_DAC/methy_result/4_MethylDackel/DAC_CpG.bedGraph")

allDat = lapply(files, function(file){
  cat("Remian", length(files) - match(file, files), "\n")
  T1 = read.bedgraph(file)
  T2 = data.frame(chr = T1[,1], pos = T1[,3], N = T1[,5] + T1[,6], X = T1[,5]) 
  return(T2)
})
head(allDat[[1]])
t <- allDat[[1]]
t[t$chr=="chr1" & t$pos==49455995,]

BSobj   <- makeBSseqData(allDat, SRX) 
dmlTest <- DMLtest(BSobj, group1 = "Ctrl",group2 = "DAC", smoothing = T)

dmls    <- callDML(dmlTest, delta = 0.1, p.threshold = 0.05)
dmrs    <- callDMR(dmlTest, delta = 0.1, p.threshold = 0.05)

fwrite(dmrs, file="/media/liyaru/LYR/301project/1_MOUSE_T_Chi_DAC/methy_result/7_DSS/1_DMR.bed",sep="\t")
fwrite(dmls, file="/media/liyaru/LYR/301project/1_MOUSE_T_Chi_DAC/methy_result/7_DSS/2_DML.bed",sep="\t")

####--------6.2 DMR annotation----------
library("ChIPseeker")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
setwd("/media/liyaru/LYR1/301project/1_MOUSE_T_Chi_DAC/methy_result/")

a <- fread("7_DSS/1_DMR.bed")
up <- a[diff.Methy < 0]
down <- a[diff.Methy > 0]

fwrite(a[,1:3],"7_DSS/bed/all_DMR.bed",col.names=F,sep="\t")
fwrite(up[,1:3],"7_DSS/bed/up_DMR.bed",col.names=F,sep="\t")
fwrite(down[,1:3],"7_DSS/bed/down_DMR.bed",col.names=F,sep="\t")

setwd("/media/liyaru/LYR1/301project/1_MOUSE_T_Chi_DAC/methy_result/7_DSS/group")
files=list.files()
mypeaks <- list() 
for (i in files){
  peak<- readPeakFile(i)
  mypeaks <- c(mypeaks,peak)
}
names(mypeaks) <- files

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=5000, downstream=5000)

options(ChIPseeker.ignore_1st_exon = T)
options(ChIPseeker.ignore_1st_intron = T)
options(ChIPseeker.ignore_promoter_subcategory = T)

tagMatrixList <- lapply(mypeaks, getTagMatrix, windows=promoter)

peakAnnoList <- lapply(mypeaks, annotatePeak, TxDb=txdb,tssRegion=c(-5000, 5000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb = "org.Mm.eg.db")
plotAnnoBar(peakAnnoList)


for (i in 1:3){
  a <-as.data.frame(peakAnnoList[i])
  colnames(a)[6:11]<-c("width","strand","length","nCG","meanMethy1","meanMethy2","diff.Methy","areaStat")
  write.table(
    a,
    paste0(names(peakAnnoList)[i],".tsv"),
    sep="\t",
    row.names = F,
    quote = F)
}


###-------6.3 KEGG---------------------------
library(enrichR)
#down DMR region with annotation See in Supplementary Table S2
peaks<- as.data.frame(fread("/media/liyaru/LYR1/301project/1_MOUSE_T_Chi_DAC/methy_result/7_DSS/DMR_anno.tsv"))
loss <- peaks [which(peaks$diff.Methy>0 & peaks$annotation=="Promoter (<=1kb)" ),"SYMBOL"]

g=loss
dbs <- c("KEGG_2019_Mouse")
kegg <- enrichr(g, dbs)
kegg <- as.data.frame(kegg)

kegg <- kegg[kegg$KEGG_2019_Mouse.P.value<=0.05,]
fwrite(kegg,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/3_DMR_430_genes_KEGG.csv")

tf <- kegg[,c(1,3)]

colnames(tf) <- c("GO terms","P-value")

tf <- tf[as.numeric(tf$`P-value`)<=0.05,]
tf$`P-value` <- -log10(as.numeric(tf$`P-value`))
rownames(tf) <- tf$`GO terms`
tf <- tf[1:9,]
tf <- tf[order(as.numeric(tf$`P-value`)),]

barplot(as.numeric(tf$`P-value`),horiz=T,xlim=c(0,max(tf$`P-value`)+0.5),axes=T,col=rgb(171,205,233, maxColorValue = 255),xlab ="-log10(p-value)",
        cex.axis=1.3,cex.lab=1.5,border = NA) 
for (i in 1:nrow(tf)){
  text(0,(1.2*i-0.6),tf$`GO terms`[i],cex=1.6,pos=4)
}

####-------6.4 DMR signal heatmap ---------------------------
library(data.table)
library(readxl)
library(EnrichedHeatmap)
library(tidyverse)
library(rtracklayer)
library(ggplot2)
library(RColorBrewer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(cowplot)
library(ggsci)
library(clusterProfiler)
blues <- colorRampPalette(colors = brewer.pal(9,"Blues"))


####--------6.4.1 type transfer------
a <- fread("/media/liyaru/LYR1/301project/1_MOUSE_T_Chi_DAC/methy_result/7_DSS/down_methy_gene.csv")

tt <- bitr(a$`gene name` ,fromType = 'SYMBOL',
           toType = c('ENTREZID'),
           OrgDb='org.Mm.eg.db')
colnames(tt)[1] <- "gene"


####-------6.4.2 DMR ranger----------
a <- fread("/media/liyaru/LYR1/301project/1_MOUSE_T_Chi_DAC/methy_result/7_DSS/1_DMR_down.bed")

t <- GRanges(seqnames = a$V1,
             ranges = IRanges(start = a$V2,end=a$V3))
target.tss <- t

####---------methylation----------------------
setwd("/media/liyaru/LYR1/301project/1_MOUSE_T_Chi_DAC/bigwig/methy_bedgraph")

a <- c("CON_CpG.bedGraph","DAC_CpG.bedGraph")

result <- list()
result_ht <- c()
names <- c("Ctrl","DAC")
j=1

for (i in a){
  #i="CON_CpG.bedGraph"
  ICM_ATAC <- import.bedGraph(i)
  mat.ATAC.tss <- normalizeToMatrix(signal = ICM_ATAC,
                                    target=target.tss,
                                    extend = 500,
                                    w = 50,
                                    value_column = "score",
                                    keep = c(0,0.95),
                                    background = NA
  )
  print(i)
  
  dim(mat.ATAC.tss)
  t <-  as.data.frame(rownames(mat.ATAC.tss))
  
  ht1 <- EnrichedHeatmap(mat.ATAC.tss,
                         col = blues(10),
                         name = paste0(names[j]),
                         axis_name = c("-500bp","start","end","500bp"),
                         column_title = names[j], 
                         top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = colors(9)),
                         ylim=c(0,80)
                         )),
                         row_title_rot = 0,
                         use_raster=F,
  )
  ht1
  result_ht <- result_ht+ht1
  j=j+1
}

result_ht


####-------6.5 Boxplot+vlnplot Methy level -------
####-------6.5.1 TSS position---------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)

mm10.g <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
mm10.g.data <- data.frame(mm10.g)
head(mm10.g.data)

mm10.tss <- promoters(mm10.g,upstream = 1000,downstream = 1000)
mm10.tss <- as.data.frame(mm10.tss)
colnames(mm10.tss)[6] <- "ENTREZID"

t <- mm10.tss$ENTREZID
tt <- bitr(t,fromType = 'ENTREZID',
           toType = c('SYMBOL'),
           OrgDb='org.Mm.eg.db')

m <- merge(mm10.tss,tt,by="ENTREZID")

fwrite(m,"/home/liyaru/public_Data/mm10_TSS_1kb_Txdb.csv")

m <- m[,c(2:4,7,1)]
fwrite(m,"/home/liyaru/public_Data/mm10_TSS_1kb_Txdb.bed",col.names = F,sep="\t")

###-----6.5.2 DAC CpG in promoter---------------------
peak.dt <- fread("/home/liyaru/public_Data/mm10_TSS_1kb_Txdb.csv")
summary(peak.dt) 
table(peak.dt$seqnames)
colnames(peak.dt)[2] <- "chrom"

methyl.dt <- fread("/media/liyaru/LYR1/301project/1_MOUSE_T_Chi_DAC/methy_result/4_MethylDackel/DAC_CpG.bedGraph", 
                   skip=1, header=FALSE, 
                   col.names=c("chrom", "start", "end", "methyl.pct", "methyl.read.count", "unmethyl.read.count"))
methyl.dt[1:5,1:5]

chroms.vector <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                   "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")

all.peaks.dt.list <- lapply(chroms.vector, FUN=function(temp.chrom){
  temp.peak.dt <- peak.dt[chrom==temp.chrom]
  temp.methyl.dt <- methyl.dt[chrom==temp.chrom]
  aa=temp.methyl.dt[, {
    index <- .I
    if ((index+1) %% 10000 == 0){
      cat(date(), " : processing index", index, "\n")
    }
    temp.peak.dt[CpG.start>= start & CpG.end <=end][,'ENTREZID']
  },  list(chrom, CpG.start=start, CpG.end=end,methyl.pct,methyl.read.count,unmethyl.read.count)]
})

result <- rbindlist(all.peaks.dt.list)
result[1:5]
fwrite(result,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/result/10_mehty/1_DAC_promoter.bed",sep = "\t")

###-----6.5.3 Ctrl CpG in promoter---------------------
peak.dt <- fread("/home/liyaru/public_Data/mm10_TSS_1kb_Txdb.csv")
summary(peak.dt) 
table(peak.dt$seqnames)
colnames(peak.dt)[2] <- "chrom"


methyl.dt <- fread("/media/liyaru/LYR1/301project/1_MOUSE_T_Chi_DAC/methy_result/4_MethylDackel/CON_CpG.bedGraph", 
                   skip=1, header=FALSE, 
                   col.names=c("chrom", "start", "end", "methyl.pct", "methyl.read.count", "unmethyl.read.count"))
methyl.dt[1:5,1:5]
summary(methyl.dt)
chroms.vector <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                   "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")
all.peaks.dt.list <- lapply(chroms.vector, FUN=function(temp.chrom){
  temp.peak.dt <- peak.dt[chrom==temp.chrom]
  temp.methyl.dt <- methyl.dt[chrom==temp.chrom]
  aa=temp.methyl.dt[, {
    index <- .I
    if ((index+1) %% 10000 == 0){
      cat(date(), " : processing index", index, "\n")
    }
    temp.peak.dt[CpG.start>= start & CpG.end <=end][,'ENTREZID']
  },  list(chrom, CpG.start=start, CpG.end=end,methyl.pct,methyl.read.count,unmethyl.read.count)]
})

result <- rbindlist(all.peaks.dt.list)
result[1:5]
fwrite(result,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/result/10_mehty/1_CON_promoter.bed",sep = "\t")


###------6.5.4 Promoter CpG boxplot-----------------------
# merge
con <- fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/result/10_mehty/1_CON_promoter.bed")
con$class <- 'Ctrl'
dac <- fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/result/10_mehty/1_DAC_promoter.bed")
dac$class <- 'DAC'

result <- rbind(con,dac)
table(result$class)

P <- ggplot(result,aes(class,`methyl.pct`))+
  geom_violin(aes(fill=class),cex=1.2)+  
  scale_fill_manual(values = color_4_sample[1:2])+
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none')+
  stat_compare_means(aes(group=class),
                     ref="Ctrl",
                     method = "wilcox.test",
                     paired=F,
                     label = "p.signif")

P

####-------7 Fig4 L-M--------------------
####-------7.1 tsne-------------------
pdf("./PDF/1.tsne.pdf",height = 6)
DimPlot(CD8,reduction = "tsne",group.by = "RNA_snn_res.0.5",pt.size = 2,
        cols = color1,
        #label = T,
        label.size = 8)
dev.off()

####----7.2 dotplot----------------
n=length(unique(CD8@meta.data$RNA_snn_res.0.5))
celltype=data.frame(ClusterID=0:(n-1),
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0),2]='Exp - c0'
celltype[celltype$ClusterID %in% c(1),2]='Ex - c1'
celltype[celltype$ClusterID %in% c(2),2]='Ex prolif. - c2'
celltype[celltype$ClusterID %in% c(3),2]='Em prolif. - c3'
celltype[celltype$ClusterID %in% c(4),2]='Ex - c4' 
celltype[celltype$ClusterID %in% c(5),2]='Active - c5' 
celltype[celltype$ClusterID %in% c(6),2]='Naive - c6'


CD8@meta.data$celltype = "unknown"
for(i in 1:nrow(celltype)){
  CD8@meta.data[which(CD8@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(CD8@meta.data$celltype)

mylevel <- c(
  'Naive - c6','Active - c5','Em prolif. - c3','Exp - c0','Ex prolif. - c2','Ex - c1','Ex - c4' 
)
CD8@meta.data$celltype <- factor(CD8@meta.data$celltype,levels=mylevel)


pdf("./PDF/2.CD8_dotplot2.pdf",height = 5)
DotPlot(CD8,
        features=c(
          "Pdcd1","Havcr2","Tigit",
          "Nkg7",
          "Gzmb","Gzmk","Mki67","Tcf7","Lef1","Il7r",
          "Tox","Nr4a2"),
        group.by = "celltype",
        scale = T)+
  RotatedAxis()& 
  scale_color_gradientn(colors = c(rdwhbu(100)[1:45],rdwhbu(100)[55:100]))
dev.off()


####------7.3 ratio------------
m <- CD8@meta.data
m$orig.ident <- factor(m$orig.ident,levels = c("Con","DAC","PD1","DP"))
m$number=1

pdf("./PDF/3.ratio.pdf")
ggplot(m,aes(orig.ident,number,fill=RNA_snn_res.0.5))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 16)+
  scale_fill_manual(
    values = color1
  )


m <- ddply(m,'orig.ident',transform,percent = 1/sum(number)*100)

ggplot(m,aes(orig.ident,percent,fill=RNA_snn_res.0.5))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 20)+
  scale_fill_manual(values = color1)

dev.off()


####------8 Fig S7 CD8+ T cells---------------
###------8.1 featureplot----------------------
pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/14.featureplot.pdf",width = 14,height = 10)
FeaturePlot(CD8,features=c("Cd3d","Cd3e",
                           "Cd8a","Cd8b1",
                           "Cd4",
                           "Pdcd1","Havcr2",
                           "Tigit", "Nr4a2",
                           "Nkg7","Gzmb","Gzmk","Mki67",
                           "Tcf7","Lef1",
                           "Il7r"),
            pt.size = 0.5,
            ncol = 4)& 
  scale_color_gradientn(colors = c(rdwhbu(100)[1:45],rdwhbu(100)[55:100]))
dev.off()

####------8.2 QC vlnplot-----------
pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/14.QC.pdf",width = 4,height = 4)
VlnPlot(CD8,features = c("nFeature_RNA"),pt.size=0,
        group.by = "sample")+
  scale_fill_manual(values = color_4_sample)
VlnPlot(CD8,features = c( "nCount_RNA"),pt.size=0,
        group.by = "sample")+
  scale_fill_manual(values = color_4_sample)
VlnPlot(CD8,features = c( "percent.mt"),pt.size=0,
        group.by = "sample")+
  scale_fill_manual(values = color_4_sample)
dev.off()

####-------8.3 phase -----------------
g2m.genes<-c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2",
             "Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2",
             "Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1",
             "Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8",
             "Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2",
             "G2e3","Gas2l3","Cbx5","Cenpa")
s.genes<-c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl",
           "Prim1","Uhrf1","Mlf1ip","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76",
           "Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1",
           "Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8")

CD8 <- CellCycleScoring(CD8, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/14.Phase.pdf",width = 5,height = 4)
DimPlot(CD8,group.by = "Phase",
        cols = color1[3:5],
        pt.size = 1)
DimPlot(CD8,group.by = "orig.ident",
        cols = color_4_sample ,
        pt.size = 1)
dev.off()


####--------9 Fig 4O Fig S8B-----------------
####-------9.1 monocle trjectory-------------
library(monocle)
Idents(CD8) <- CD8$RNA_snn_res.0.5

CD8_2 <- subset(CD8,idents=c("0","1","2","4"))

expr <- GetAssayData(CD8_2,assay = "RNA",slot = "count")
gene_anno<-data.frame(gene_id=rownames(expr),gene_short_name=rownames(expr))
rownames(gene_anno) <- rownames(expr)
sample_sheet<- CD8_2@meta.data

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_anno)

HSMM <- newCellDataSet(expr,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily=negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

disp_table <- dispersionTable(HSMM)

disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id

HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)


HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')


HSMM <- orderCells(HSMM)

saveRDS(HSMM,"/media/liyaru/LYR1/301project/data/monocle.rds")

plot_cell_trajectory(HSMM, color_by = "RNA_snn_res.0.5",cell_size = 1,
                     show_tree = T)+
  scale_colour_manual(values = color1[c(1,2,3,5)])

plot_cell_trajectory(HSMM, color_by = "RNA_snn_res.0.5",
                     show_branch_points = F) + 
  facet_wrap("~RNA_snn_res.0.5", nrow = 1)+
  scale_colour_manual(values = color1[c(1,2,3,5)])

plot_cell_trajectory(HSMM, color_by = "RNA_snn_res.0.5",
                     show_branch_points = F) + 
  facet_wrap("~orig.ident", nrow = 1)+
  scale_colour_manual(values = color1[c(1,2,3,5)])

plot_cell_trajectory(HSMM,color_by="Pseudotime", size=1,show_backbone=TRUE)& 
  scale_color_gradientn(colors = c(rdwhbu(100)[1:45],rdwhbu(100)[55:100]))

####--------9.2 monocle gene exp-------------
markers <-c("Il7r","Tcf7","Lef1","Gzmk","Prf1","Cx3cr1","Gzmb","Havcr2","Pdcd1","Tox","Tigit")

my_genes <- markers

t <- HSMM@phenoData@data
t <- t[t$RNA_snn_res.0.5 != "2",] 
t <- rownames(t)

#drop C2(proliferating)
cds_subset <- HSMM[my_genes,t]

#each gene in trajectory
#pdf("./PDF/5.monocle_gene.pdf",height = 8,width = 9)
plot_genes_in_pseudotime(cds_subset, 
                         color_by="RNA_snn_res.0.5",
                         ncol=2)+
  scale_colour_manual(values = color1[c(1,2,3,5)])
#dev.off()

#pdf("./PDF/5.monocle_gene_pseudotime.pdf",height = 8,width = 9)
plot_genes_in_pseudotime(cds_subset, 
                         color_by="Pseudotime",
                         ncol=2)& 
  scale_color_gradientn(colors = c(rdwhbu(100)[1:45],rdwhbu(100)[55:100]))
#dev.off()

#gene in pseudotime heatmap
#pdf("./PDF/5.monocle_gene_heatmap.2.pdf",height = 3,width = 3)
plot_pseudotime_heatmap(cds_subset,
                        show_rownames = T,
                        cluster_rows = F)
#dev.off()

####--------9.3 add annotation to monocle gene exp -------------
#add annotation col in pseudotime heatmap
#from plot_pseudotime_heatmap function
num_clusters <- min(num_clusters, nrow(cds_subset))
pseudocount <- 1

#change length.out = length(cds_subset$Pseudotime)
#each cell as a column rather than only have 100 coloum
newdata <- data.frame(Pseudotime = seq(min(cds_subset$Pseudotime),
                                       max(cds_subset$Pseudotime), 
                                       #length.out = 100,
                                       length.out = length(cds_subset$Pseudotime)))
# get gene exp 
trend_formula = "~sm.ns(Pseudotime, df=3)"
m <- genSmoothCurves(cds_subset, cores = 1, 
                     trend_formula = trend_formula, 
                     relative_expr = T, new_data = newdata)

m = m[!apply(m, 1, sum) == 0, ]

# norm & scale
norm_method = c("log", "vstExprs")
norm_method <- match.arg(NULL,norm_method)

if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == 
    FALSE) {
  m = vstExprs(cds_subset, expr_matrix = m)
}else if (norm_method == "log") {
  m = log10(m + pseudocount)
}

m = m[!apply(m, 1, sd) == 0, ]
m = Matrix::t(scale(Matrix::t(m), center = TRUE))
m = m[is.na(row.names(m)) == FALSE, ]
m[is.nan(m)] = 0
m[m > scale_max] = scale_max
m[m < scale_min] = scale_min

heatmap_matrix <- m

row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1
hmcols = NULL

if (is.null(hmcols)) {
  bks <- seq(-3.1, 3.1, by = 0.1)
  hmcols <- monocle:::blue2green2red(length(bks) - 1)
}else {
  bks <- seq(-3.1, 3.1, length.out = length(hmcols))
}

df <- cds_subset@phenoData@data
df <- df[,c("Pseudotime", "RNA_snn_res.0.5","orig.ident")]
df <- df[order(df$Pseudotime, decreasing = F),]

tt <- df
rownames(tt) <- colnames(heatmap_matrix)
colnames(tt) <- c("Pseudotime","Cluster","Sample")
summary(tt)

Cluster_color <- color1[c(1,2,5)]
names(Cluster_color) <- c("0","1","4")

Sample_color <- color_4_sample
names(Sample_color) <- levels(tt$Sample)

ann_colors = list(
  Pseudotime = rdwhbu(50),
  Cluster = Cluster_color,
  Sample = Sample_color
)

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0711/2.pheatmap.pdf",height = 4,width = 5)
ph_res <- pheatmap(heatmap_matrix,
                   #useRaster = T,
                   cluster_cols = F, 
                   cluster_rows = F, 
                   show_rownames = T, 
                   show_colnames = F, 
                   annotation = tt,
                   # treeheight_row = 20, breaks = bks, fontsize = 6, 
                   color = hmcols, 
                   #border_color = NA, silent = TRUE, filename = NA
                   annotation_colors=ann_colors
)
ph_res
dev.off()

####---------10 FigS8C-D------------------------
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(destiny)
library(rgl)

####-------10.1 destiny diffusion map -----------------
ct <-GetAssayData(CD8,assay = "RNA",slot = "data")
ct<-ct[VariableFeatures(CD8),]
ct <- as.ExpressionSet(as.data.frame(t(ct)))

#add annotation
ct$celltype <- CD8@meta.data[,c("RNA_snn_res.0.5")]
dm <- DiffusionMap(ct,k = 10)

dpt <- DPT(dm,tips=1)
saveRDS(dm,"/media/liyaru/LYR1/301project/data/destiny.rds")
saveRDS(dpt,"/media/liyaru/LYR1/301project/data/destiny_diffusion.map.rds")


gg <- TSNEPlot(CD8,group.by="RNA_snn_res.0.5",cols=color1)
gg
c <- ggplot_build(gg)$data[[1]]$colour

plot.DPT(dpt,col_by = "RNA_snn_res.0.5",col=c)

plot3d(
  eigenvectors(dm)[,1:3],
  col = c,
  size = 10,
  type = 'p'
)

####--------10.2 get embeddings from destiny -------------------
t <- CD8@reductions$tsne@cell.embeddings
tt <- as.data.frame(t)

t[,1] <- dm$DC1[rownames(tt)]
t[,2] <- dm$DC2[rownames(tt)]

CD8@reductions$tsne@cell.embeddings <- t

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0707/diffusion.map.pdf",width = 5,height = 5)
DimPlot(CD8,group.by = "RNA_snn_res.0.5",
        cols = color1,
        pt.size = 1)
dev.off()

####-------10.3 velocyto ------------------
ldat <- ReadVelocity(file = "/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/1_cellranger_result/RNA/Con-5-lib/velocyto/Con-5-lib.loom")
bm <- as.Seurat(x = ldat)

Idents(CD8) <-  CD8@meta.data$orig.ident
Ctrl <- subset(CD8,idents="Con")

#barcode name 
a <- gsub("x","-1_1",colnames(bm$spliced))
a <- gsub("Con-5-lib:","",a)

ttt <- bm@assays$spliced
colnames(ttt)
colnames(ttt@counts) <- a
colnames(ttt@data) <- a
bm@assays$spliced <- ttt

ttt <- bm@assays$unspliced
colnames(ttt)
colnames(ttt@counts) <- a
colnames(ttt@data) <- a
bm@assays$unspliced <- ttt

ttt <- bm@assays$ambiguous
colnames(ttt)
colnames(ttt@counts) <- a
colnames(ttt@data) <- a
bm@assays$ambiguous<- ttt

# filter cell (use same cell in seurat)
rownames(bm@meta.data) <- gsub("x","-1_1",rownames(bm@meta.data))
rownames(bm@meta.data) <- gsub("Con-5-lib:","",rownames(bm@meta.data))

sp <- bm$spliced[,rownames(Ctrl@meta.data)]
unsp <- bm$unspliced[,rownames(Ctrl@meta.data)]
WTumap <- Ctrl@reductions$tsne@cell.embeddings

cell.dist <- as.dist(1-armaCor(t(WTumap)))
fit.quantile <- 0.02

rvel.cd <- gene.relative.velocity.estimates(sp,unsp,deltaT=2,kCells=10,
                                            cell.dist=cell.dist,fit.quantile=fit.quantile,
                                            n.cores=10)

gg <- TSNEPlot(Ctrl,group.by="RNA_snn_res.0.5",label=T,label.size=10,cols=color1)
gg
ggplot_build(gg)$data
colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(WTumap)

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0707/Velocyto.pdf",width = 5,height = 5)
show.velocity.on.embedding.cor(WTumap,
                               rvel.cd,
                               n=100,
                               scale='sqrt',
                               cell.colors=ac(colors,alpha=0.5),
                               cex=1.2,
                               arrow.scale=0.08,
                               show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5,
                               grid.n=20,
                               arrow.lwd=1,
                               do.par=F,cell.border.alpha = 0.1)
dev.off()
saveRDS(bm,"/media/liyaru/LYR1/301project/data/velocyto.rds")


####-----11 Fig5A FigS9A --------------
####-----11.1 immunarch--------
library(immunarch)  # Load the package into R

# # get CD8 TCR 
# {
#   file_path = "/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/1_cellranger_result/TCR/merge/"
#   cells <- colnames(CD8)
#   #整合各个tcr文件 在barcode上加上相应的序号
#   sup <- c("-1_1","-1_2","-1_3","-1_4")
#   j=1
#   for (i in c("Con","DAC","DP","PD1")){
#     #i="Con"
#     tcr <- read.csv(paste0(file_path,i,".csv"))
#     tcr$barcode <- gsub("-1", "", tcr$barcode)
#     tcr$barcode <- paste0(tcr$barcode,sup[j])
#     tcr <- tcr[tcr$barcode %in% cells,]
#     fwrite(tcr,paste0("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/1_cellranger_result/TCR/CD8/",i,".csv"))
#     j=j+1
#   }
# }

file_path = "/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/1_cellranger_result/TCR/CD8/"
immdata_10x <- repLoad(file_path)

names(immdata_10x$data)

####-----11.2 diversity--------
# Hill numbers
div_hill <- repDiversity(immdata_10x$data, "hill")
p1 <- vis(div_hill, .by = c("Sample"), .meta = immdata_10x$meta)
p1 + scale_color_manual(values=c(color_4_sample))


# D50
div_d50 <- repDiversity(immdata_10x$data, "d50")
p2 <- vis(div_d50)
p2+scale_fill_manual(values=c(color_4_sample))


####----11.3 top proportion----------------
#FigS9I
imm_top <- repClonality(immdata_10x$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_top
imm_top %>% vis()



####-----12 Fig5B-D --------------------------
####-----12.1 CLone size----------------------
CD8@meta.data$sample_clonetype <- paste0(CD8@meta.data$orig.ident,"_",CD8@meta.data$clonotype_id)
CD8@meta.data$cluster_sample_clonetype <- paste0(CD8@meta.data$RNA_snn_res.0.5,"_",CD8@meta.data$orig.ident,"_",CD8@meta.data$clonotype_id)

CD8@meta.data[CD8@meta.data$frequency >= 10, 'freq']="10+"
CD8@meta.data[CD8@meta.data$frequency >= 2 & CD8@meta.data$frequency < 10, 'freq']="2~9"
CD8@meta.data[CD8@meta.data$frequency == 1, 'freq']="1"
CD8@meta.data$freq <- factor(CD8@meta.data$freq,levels = c("1","2~9","10+"))

###-----12.2 clone size TSNE ----------------
Idents(CD8) <- CD8$RNA_snn_res.0.5
CD8_0 <- subset(CD8,ident="0")

DimPlot(CD8,group.by = "freq",pt.size = 1,split.by = "orig.ident",
        cols = color4,ncol = 4)
DimPlot(CD8_0,group.by = "freq",pt.size = 1,split.by = "orig.ident",
        cols = color4,ncol = 4)

###-----12.3 size ratio by cluster -----------
# by cell
m <- CD8@meta.data
m$number=1

pdf("./PDF/TCR/1.clonetype_by_cluster.pdf")
ggplot(m,aes(RNA_snn_res.0.5,number,fill=freq))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 16)+
  scale_fill_manual(
    values = color3[c(1,3,5)]
  )

m <- ddply(m,'RNA_snn_res.0.5',transform,percent = 1/sum(number)*100)
m$percent
ggplot(m,aes(RNA_snn_res.0.5,percent,fill=freq))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 20)+
  scale_fill_manual(values = color3[c(1,3,5)])

dev.off()

# clonotype by Cluster
# FigS9 C-D
m <- CD8@meta.data
m <- m[!duplicated(m$cluster_sample_clonetype),]
m$number=1
ggplot(m,aes(RNA_snn_res.0.5,number,fill=freq))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 16)+
  scale_fill_manual(
    values = color3[c(1,3,5)]
  )

m <- ddply(m,'RNA_snn_res.0.5',transform,percent = 1/sum(number)*100)
m$percent
ggplot(m,aes(RNA_snn_res.0.5,percent,fill=freq))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 20)+
  scale_fill_manual(values = color3[c(1,3,5)])

dev.off()


# cell by sample
# FigS9 E-F
m <- CD8@meta.data
m$number=1
ggplot(m,aes(orig.ident,number,fill=freq))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 16)+
  scale_fill_manual(
    values = color3[c(1,3,5)]
  )

m <- ddply(m,'orig.ident',transform,percent = 1/sum(number)*100)
m$percent
ggplot(m,aes(orig.ident,percent,fill=freq))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 20)+
  scale_fill_manual(values = color3[c(1,3,5)])

dev.off()


###-----12.4 highly expanded (size > 10)------
m$cluster_clone <- paste0(m$RNA_snn_res.0.5,"_",m$sample_clonetype)
table(m$cluster_clone)

m <- m[m$frequency >=10,]
table(m$orig.ident)

m$number=1

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF/TCR/5.more_than_10_by_sample_cluster.pdf")
ggplot(m,aes(orig.ident,number,fill=RNA_snn_res.0.5))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 16)+
  scale_fill_manual(
    values = color1
    #values = my_color_palette
  )+
  RotatedAxis()

m <- ddply(m,'orig.ident',transform,percent = 1/sum(number)*100)
colnames(m)
m$percent
ggplot(m,aes(orig.ident,percent,fill=RNA_snn_res.0.5))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 20)+
  scale_fill_manual(values = color1)+
  RotatedAxis()

dev.off()

####------13 Fig 5F ----------------------------
###-------13.1 top 10 clone by sample-----------
pdf("./PDF/TCR/6.top10_clone.pdf",width = 10)
for (i in c("Con","DAC","PD1","DP")){

  m <- CD8@meta.data
  m <- m[m$orig.ident==i,]
  t <- table(m$sample_clonetype)
  t <- t[order(-t)]
  tt <- t[1:10]
  
  cc <- names(tt)
  
  m <- m[m$sample_clonetype %in% cc,]

  m$sample_clonetype <- factor(m$sample_clonetype,levels = cc)
  m$number=1

  ttt <- ggplot(m,aes(sample_clonetype,number,fill=RNA_snn_res.0.5))+
    geom_bar(stat="identity",position="stack")+
    theme_classic(base_size = 16)+
    scale_fill_manual(
      values = color1)+
    RotatedAxis()
  print(ttt)
}
dev.off()


####------13 Fig 5E ----------------------------
####------13.1 heatmap--------------------------
a <- fread("/media/liyaru/LYR1/301project/Table/target_gene_2.csv",header = F)
a <- a$V1
a <- tolower(a)
a <- capitalize(a)
a[27] <- 'Tbx21'

CD8@meta.data$freq <- as.character(CD8@meta.data$freq)
CD8@meta.data[CD8@meta.data$frequency >= 10, 'freq']="10+"
CD8@meta.data[CD8@meta.data$frequency >= 2 & CD8@meta.data$frequency < 10, 'freq']="2~10"
CD8@meta.data[CD8@meta.data$frequency == 1, 'freq']="1"

CD8@meta.data$sample_freq <- paste0(CD8@meta.data$orig.ident,"_",CD8@meta.data$freq)
table(CD8@meta.data$sample_freq)

t <- AverageExpression(CD8,assays = "RNA",features = a,group.by = "sample_freq")
t <- t[["RNA"]]
t <- as.data.frame(t)
t <- t[,c("Con_2~10","DAC_2~10","PD1_2~10","DP_2~10","Con_10+","DAC_10+","PD1_10+","DP_10+")]
 
pheatmap(t,scale = "row",
         cluster_rows=F,cluster_cols = F,
         color=rdwhbu(100))


####------14 FigS9H S9J ---------------------
library(clusterProfiler)
library(data.table)
library(org.Mm.eg.db)
library(openxlsx)

####------14.1 GO---------------------
CD8@meta.data$sample_freq <- paste0(CD8@meta.data$orig.ident,"_",CD8@meta.data$freq)
DefaultAssay(CD8) <- "RNA"
CD8.markers <- FindMarkers(CD8,
                           group.by = "sample_freq",
                           ident.1 = "DP_10+",
                           ident.2 ="PD1_10+",
                           #test.use = "MAST"
                           logfc.threshold=0
)
#fwrite(CD8.markers,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0620/12.expand10.DEG.csv",row.names = T)


a <- CD8.markers
a <- a[a$p_val< 0.05 & a$avg_log2FC > 0,]

ego <- enrichGO(
  gene          = rownames(a),
  keyType       = "SYMBOL",
  OrgDb         =  org.Mm.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  #readable      = TRUE
)

clusterProfiler::dotplot(ego, showCategory = 10)+
  scale_color_gradientn(colors = c(rdwhbu_re(100)[1:45],rdwhbu_re(100)[55:100]))

####------14.2 C0 Top50 ---------------------
m <- CD8@meta.data
#unique(m$sample_clonetype)
m <- m[m$RNA_snn_res.0.5 == "0",]
t <- table(m$sample_clonetype)
t <- t[order(-t)]
tt <- t[1:50]

cc <- names(tt)
m <- m[m$sample_clonetype %in% cc,]

m$sample_clonetype <- factor(m$sample_clonetype,levels = cc)

m$number=1

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0620/14.C0.clonotype.pdf",width = 10,height = 6)
ttt <- ggplot(m,aes(sample_clonetype,number,fill=orig.ident))+
  geom_bar(stat="identity",position="stack")+
  theme_classic(base_size = 16)+
  scale_fill_manual(
    values = color_4_sample
    #values = my_color_palette
  )+
  RotatedAxis()
ttt
dev.off()
print(ttt)

####------15 Fig6A-B ---------------------
####------15.1 DEG-----------------------
CD8@meta.data$cluster_sample <- paste0(CD8@meta.data$RNA_snn_res.0.5,"_",CD8@meta.data$orig.ident)
table(CD8$cluster_sample)
DefaultAssay(CD8) <- "RNA"

#DP VS P DEG
CD8.markers <- FindMarkers(CD8,
                           group.by = "cluster_sample",
                           ident.1 = "0_DP",
                           ident.2 ="0_PD1",
                           #test.use = "MAST"
                           logfc.threshold=0)
fwrite(CD8.markers,
       "/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/C0.DP_PD1.csv",
       row.names =T )


####------15.2 volcano-----------------------
library(ggplot2)
library(ggrepel)
library(org.Mm.eg.db)
a <- fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/C0.DP_PD1.csv")

target <- c("Jund","Eif3f","H2afj","Tomm20","Ube2s","Pgls","Ets1","Eif5","Xbp1","Grina","Dusp2")

Dat<- as.data.frame(a)
Dat$gene <- Dat$V1

Dat$threshold <- 'Other'
Dat[Dat$p_val < 0.05 & Dat$avg_log2FC > 0,'threshold'] <- 'Up'
Dat[Dat$gene %in% target,'threshold'] <- 'Label'
Dat[Dat$p_val < 0.05 & Dat$avg_log2FC < 0,'threshold'] <- 'Down'
table(Dat$threshold)

Dat$threshold <- factor(Dat$threshold,levels = c('Up','Down','Label','Other'))
show_col(c(  "#EEBBBB","#AAAAD4", "#CD3333","#808080"))
pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/4_volcano_DP_vs_P_C0.pdf",width = 5,height = 4)
ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val),color=threshold))+
  geom_point()+
  #scale_color_manual(values=c("#DC143C","#00008B",'black',"#808080"))+
  #scale_color_manual(values=c( "#FB6A4A","#67000D","#9ECAE1","#08519C","#808080"))+
  #scale_color_manual(values=c(  "#EEBBBB","#CD3333","#AAAAD4", "#000080","#808080"))+
  scale_color_manual(values=c(  "#EEBBBB","#AAAAD4", "#CD3333","#808080"))+
  geom_text_repel(
    #data = Dat[(Dat$p_val_adj<0.05 & abs(Dat$avg_log2FC) > fc1)| Dat$V1 %in% c("Jund","Jun") ,],
    #data = Dat[(Dat$p_val_adj<0.05 & abs(Dat$avg_log2FC) > fc),],
    data = Dat[Dat$gene %in% target,],
    aes(label = gene),
    size = 5,
    segment.color = "black", show.legend = FALSE,
    max.overlaps=40)+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    #text = element_text(family="Arial",size = 17),
    text = element_text(size = 17),
    # axis.title.x= element_text(size=14 , family="Arial"),
    # axis.title.y = element_text(size = 14, family="Arial"),
    # axis.text = element_text(size = 14, family="Arial")
  )+
  ylab('-log10 (p-value)')+
  xlab('log2 (FoldChange)')+
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)
dev.off()

####------15.3 GO of all upregulated DEGs-----------------------
library(org.Mm.eg.db)
library(clusterProfiler)

a <- fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/C0.DP_PD1_0.05.csv")
a <- a[a$avg_log2FC > 0]

ego <- enrichGO(
  gene          = a$V1,
  keyType       = "SYMBOL",
  OrgDb         =  org.Mm.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  #readable      = TRUE
)

terms <-c ("T cell differentiation",
           "ribonucleoprotein complex biogenesis",
           "protein folding",
           "cellular response to chemical stress",
           "regulation of T cell activation",
           "nuclear transport",
           "regulation of immune effector process",
           "regulation of adaptive immune response",
           "regulation of DNA metabolic process",
           "protein polyubiquitination",
           "regulation of mRNA metabolic process",
           "T cell proliferation",
           "cell killing",
           "mitochondrion organization"
)

t <- ego@result
#fwrite(t,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/11.C0.DP.GO.csv",row.names = T)
t <- t[t$Description %in% terms,]
t <- t[as.numeric(t$pvalue)<=0.05,]
t$`P-value` <- -log10(as.numeric(t$p.adjust))
rownames(t) <- t$Description
t <- t[order(as.numeric(t$`P-value` )),]

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/11.C0.DP.GO.pdf",width = 7,height = 5)
barplot(as.numeric(t$`P-value`),
        horiz=T,
        xlim=c(0,max(t$`P-value`)+0.5),
        axes=T,
        #col="lightblue",
        col = "#EEBBBB",
        xlab ="-log10(p.adjust)",
        cex.axis=1.3,cex.lab=1,border = NA) 
for (i in 1:nrow(t)){
  text(0,(1.2*i-0.6),t$Description[i],cex=1.2,pos=4)
}
dev.off()


####--------16 Fig6C---------------
####--------16.1 heatmap----------
library(pheatmap)
library(RColorBrewer)
setwd("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP")

DefaultAssay(CD8) <- "RNA"
CD8@meta.data$cluster_sample <- paste0(CD8@meta.data$RNA_snn_res.0.5,"_",CD8@meta.data$orig.ident)
a<-AverageExpression(CD8,group.by = "cluster_sample",slot = "data")
a <- a[["RNA"]]

t <-as.data.frame(fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/C0.DP_PD1.csv"))  
up <- t[t$avg_log2FC > 0 & t$p_val < 0.05,'V1']
aa <- a[up,]
colnames(aa)
aa <- aa[,c("0_Con","0_PD1","0_DP","0_DAC")]

# show genes
genes <- c("NFATC3","MAPK1","SOCS1","IFNG",
           "JAK1","AKT2","RUNX2","RUNX3","ZBTB1",
           "MAP2K","UBE2B","IL7R","MEF2D","IKZF1",
           "NFKB1","NFKB2","PRF1","GZMK",
           "GZMA","GZMB","JUND","ETS1","CD47",
           "Fyn",
           "Eif5","Eif3b","Mapkapk3",
           "Xbp1","Atp1b3","Atp5md","Cox7c","Tomm7",
           "Ndufa3",
           "Arf5",
           "Tomm20","Cox5a","Pak2","Arf6")
genes <- genes %>%
  tolower() %>%
  Hmisc::capitalize()


t <- t(scale(t(aa)))

p1 <- pheatmap(t,
               border_color=NA,
               color = rdwhbu(100),
               #scale = "row",
               #cutree_rows = N,
               cluster_row = T,
               cluster_col = F,
               show_rownames=T,
               show_colnames=T,
               clustering_distance_rows = "euclidean",
               clustering_method='complete',
               #breaks = breaksList,
               #annotation_row=annotation_row,
               #annotation_colors=ann_colors
)

###------16.2 DEG Module------------
d = dist(t, method = 'euclidean')
tree = hclust(d, method = 'complete')

N=7
v = cutree(tree, N)[tree$order]
gaps = which((v[-1] - v[-length(v)]) != 0)

gene.cluster <- as.data.frame(v)
gene.cluster$gene <- rownames(gene.cluster)
table(gene.cluster$v)

{
  annotation_row <- data.frame(Cut = gene.cluster$v,
                               gene = rownames(gene.cluster)
  )
  
  rownames(annotation_row) <- annotation_row$gene
  
  annotation_row[annotation_row$Cut == "2",'Module'] = 'G1'
  annotation_row[annotation_row$Cut == "7",'Module'] = 'G2'
  annotation_row[annotation_row$Cut == "6",'Module'] = 'G3'
  annotation_row[annotation_row$Cut == "4",'Module'] = 'G4'
  annotation_row[annotation_row$Cut == "3",'Module'] = 'G5'
  annotation_row[annotation_row$Cut == "5",'Module'] = 'G6'
  annotation_row[annotation_row$Cut == "1",'Module'] = 'G7'
}

annotation_row <- annotation_row[,c(2:3)]
annotation_row[!(annotation_row$gene %in% genes),'gene']=""
annotation_row$gene <- factor(annotation_row$gene,levels = c(unique(annotation_row$gene),NA))

c1 <- brewer.pal(N, "Set2")
names(c1) <- unique(annotation_row$Module)

colnames(annotation_row)

ann_colors = list(
  Module = c1)

p2 <- pheatmap(t[,c(1,4,2,3)],
               border_color=NA,
               color = rdwhbu(100),
               #scale = "row",
               cutree_rows = N,
               cluster_row = T,
               cluster_col = F,
               show_rownames=F,
               show_colnames=T,
               clustering_distance_rows = "euclidean",clustering_method='complete',
               #breaks = breaksList,
               annotation_row=annotation_row,
               annotation_colors=ann_colors
)

print(p2)

colnames(annotation_row)[1] <- "Labeled" 
annotation_row$gene <- rownames(annotation_row)
fwrite(annotation_row,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/13.up.DEG.Module.csv")

DEG <- readxl::read_excel("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/4_C0.DP_PD1_up.xlsx")
m <- merge(DEG,annotation_row,by="gene",all.x=T,all.y=T)
m <- m[order(m$Module),]
fwrite(m,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/13.up.DEG.Module.FC.P.all.csv")


####-------17 Fig6D FigS10A ---------------------
####-------17.1 UP genes module GO---------------
#split genes in each module
a <- readxl::read_excel("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/13.up.DEG.Module.FC.P.all.xlsx")
table(a$Module)
target <- NULL
for(i in unique(a$Module)){
  #i=="G1"
  ii <- a[a$Module==i,]
  ii <- ii$gene
  target <- c(target,list(ii))
}
names(target) <- paste0("Module",1:7)


#GO for each module
for (i in c(1:7)){
  ii=target[[i]]  
  ego <- enrichGO(
    gene          = ii,
    keyType       = "SYMBOL",
    OrgDb         =  org.Mm.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    #readable      = TRUE
  )
  
  kegg <- ego@result
  
  fwrite(kegg,paste0("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF/DEG/KEGG/Cut",i,".GO.csv"))
  
  pdf(paste0("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF/DEG/KEGG/Cut",i,".GO.pdf"))
  
  { # barplot
    tf <- kegg[,c(3,6)]
    colnames(tf) <- c("GO terms","P-value")
    tf <- tf[as.numeric(tf$`P-value`)<=0.05,]
    tf$`P-value` <- -log10(as.numeric(tf$`P-value`))
    rownames(tf) <- tf$`GO terms`
    tf <- tf[1:10,]
    
    tf <- tf[order(as.numeric(tf$`P-value`)),]
    
    p <- barplot(as.numeric(tf$`P-value`),horiz=T,xlim=c(0,max(tf$`P-value`)+0.5),axes=T,
                 col=color7[i],
                 xlab ="-log10(p-value)",
                 cex.axis=1.3,cex.lab=1.5,border = NA) 
    p
    for (iii in 1:nrow(tf)){
      p+text(0,(1.2*iii-0.6),tf$`GO terms`[iii],cex=1.6,pos=4)
    }
  }
  print(p)
  dev.off()
  
}

#dotplot each module of selected GO terms 
setwd("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF/DEG/KEGG")
a <- list.files("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF/DEG/KEGG",pattern="GO.csv")
a
result <- as.data.frame(matrix(nrow = 0,ncol = 11))
for (i in 1:7){
  #i=1
  file <- paste0("Cut",i,".GO.csv")
  t <- as.data.frame(fread(file))
  t$Module <- i
  result <- rbind(result,t)
}
table(result$Module)

term <- fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/5_selected_GO.csv")

term[term$Cut == "2",'Module'] = 'G1'
term[term$Cut == "7",'Module'] = 'G2'
term[term$Cut == "6",'Module'] = 'G3'
term[term$Cut == "4",'Module'] = 'G4'
term[term$Cut == "3",'Module'] = 'G5'
term[term$Cut == "5",'Module'] = 'G6'
term[term$Cut == "1",'Module'] = 'G7'

term <- term[order(term$Module,term$p.adjust),]
term <- term$Description

result <- result[result$Description %in% term,]
result$Module <- paste0("G",result$Module)

result <- separate(result,col = GeneRatio,into = c("R1","R2"),sep = "[/]",remove = F)
result$Ratio <- as.numeric(result$R1) / as.numeric(result$R2)

result$Description <- factor(result$Description,levels = term[length(term):1])

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/5_Module_7_selected.GO.pdf")
ggplot(result,aes(x = Module,y = Description))+
  geom_point(aes(
    #color = pvalue,
    color=p.adjust,
    size=Ratio))+
  scale_size()+
  scale_color_gradientn(colors = c(rdwhbu_re(100)[1:45],rdwhbu_re(100)[55:100]))+
  theme_bw()
dev.off()

####-------17.2 UP genes module KEGG---------------
#split genes in each module
a <- readxl::read_excel("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/13.up.DEG.Module.FC.P.all.xlsx")
table(a$Module)
target <- NULL
for(i in unique(a$Module)){
  #i=="G1"
  ii <- a[a$Module==i,]
  ii <- ii$gene
  target <- c(target,list(ii))
}
names(target) <- paste0("Module",1:7)

# each module KEGG
dbs <- c("KEGG_2019_Mouse")
color7 <- brewer.pal(7, "Set2")
color7 <- c("#E5C494","#66C2A5","#A6D854","#E78AC3","#FFD92F","#8DA0CB","#FC8D62")
show_col(color7)

for (i in c(1:7)){
  #i=1
  ii=target[[i]]
  
  #KEGG
  kegg <- enrichr(ii, dbs)
  kegg <- as.data.frame(kegg)
  kegg <- kegg[kegg$KEGG_2019_Mouse.P.value<=0.05,]
  
  fwrite(kegg,paste0("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF/DEG/KEGG/Cut",i,".KEGG.csv"))
  
  pdf(paste0("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF/DEG/KEGG/Cut",i,".KEGG.pdf"))
  
  { # barplot
    tf <- kegg[,c(1,3)]
    colnames(tf) <- c("GO terms","P-value")
    tf <- tf[as.numeric(tf$`P-value`)<=0.05,]
    tf$`P-value` <- -log10(as.numeric(tf$`P-value`))
    rownames(tf) <- tf$`GO terms`
    tf <- tf[1:10,]
    
    tf <- tf[order(as.numeric(tf$`P-value`)),]
    
    p <- barplot(as.numeric(tf$`P-value`),horiz=T,xlim=c(0,max(tf$`P-value`)+0.5),axes=T,
                 col=color7[i],
                 xlab ="-log10(p-value)",
                 cex.axis=1.3,cex.lab=1.5,border = NA) 
    p
    for (iii in 1:nrow(tf)){
      p+text(0,(1.2*iii-0.6),tf$`GO terms`[iii],cex=1.6,pos=4)
    }
  }
  print(p)
  dev.off()
} 

#dotplot KEGG selected term in each module
setwd("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF/DEG/KEGG")
a <- list.files("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF/DEG/KEGG",pattern="KEGG.csv")
a
result <- as.data.frame(matrix(nrow = 0,ncol = 10))
for (i in 1:7){
  #i=1
  file <- paste0("Cut",i,".KEGG.csv")
  t <- as.data.frame(fread(file))
  t$Module <- i
  result <- rbind(result,t)
}
table(result$Module)
colnames(result) <- gsub("KEGG_2019_Mouse[.]","",colnames(result))

term <- readxl::read_excel("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/5_selected_KEGG.xlsx") %>% as.data.frame(term)

term[term$Cut == "2",'Module'] = 'G1'
term[term$Cut == "7",'Module'] = 'G2'
term[term$Cut == "6",'Module'] = 'G3'
term[term$Cut == "4",'Module'] = 'G4'
term[term$Cut == "3",'Module'] = 'G5'
term[term$Cut == "5",'Module'] = 'G6'
term[term$Cut == "1",'Module'] = 'G7'

term <- term[order(term$Module,term$KEGG_2019_Mouse.Adjusted.P.value),]
term <- term$KEGG_2019_Mouse.Term

result <- result[result$Term %in% term,]
result <- result[result$Adjusted.P.value < 0.05,]
result$Module <- paste0("G",result$Module)

result$Term <- factor(result$Term,levels = unique(term[length(term):1]))

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/5_Module_7_selected.KEGG.pdf")
ggplot(result,aes(x = Module,y = Term))+
  geom_point(aes(
    color=Adjusted.P.value,
    size=Odds.Ratio))+
  scale_size()+
  scale_color_gradientn(colors = c(rdwhbu_re(100)[1:45],rdwhbu_re(100)[55:100]))+
  theme_bw()
dev.off()

####------18 FigS10B-----------------
library(tidyr)
library(reshape2)

t <- readxl::read_excel("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/13.up.DEG.Module.FC.P.all.xlsx") %>% as.data.frame()
t <- t[,c("gene","Module")]
tt <- dcast(t,gene~Module,value.var = "gene")

#output to metascape
fwrite(tt[,2:8],"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/Metascape/Module7.genelist.csv")

#https://metascape.org/gp/index.html

####------19 Fig6E FigSB-E-----------
# caculated by HOMER findMotifsGenome.pl 
# see in ATAC-seq_data_processing 8 homer

####------19.1 DP vs PD1 volcano-------------
library(ggplot2)
library(ggrepel)
library(data.table)
library(tidyverse)

setwd("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/ATAC_result/7_peak_new/4_homer/2_motif")

DP_gain <- fread("gain/knownResults.txt")

DP_gain$FC <- (as.numeric(gsub("%","",DP_gain$`% of Target Sequences with Motif`)))  / (as.numeric(gsub("%","",DP_gain$`% of Background Sequences with Motif` ))) 

colnames(DP_gain) <- c("Motif Name","Consensus","P-value","Log P-value","q-value (Benjamini)",
                       "# of Target Sequences with Motif",
                       "% of Target Sequences with Motif",
                       "# of Background Sequences with Motif",
                       "% of Background Sequences with Motif",
                       "FC")
fwrite(DP_gain,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/7_DP_vs_P_gain.csv")

DP_gain <- DP_gain[,c(1,4,5,10)]

DP_loss <- fread("loss/knownResults.txt")
DP_loss$FC <- (as.numeric(gsub("%","",DP_loss$`% of Target Sequences with Motif`)))  / (as.numeric(gsub("%","",DP_loss$`% of Background Sequences with Motif` ))) 
summary(DP_loss$FC)
DP_loss$FC <- -DP_loss$FC
colnames(DP_loss) <- c("Motif Name","Consensus","P-value","Log P-value","q-value (Benjamini)",
                       "# of Target Sequences with Motif",
                       "% of Target Sequences with Motif",
                       "# of Background Sequences with Motif",
                       "% of Background Sequences with Motif",
                       "FC")
fwrite(DP_loss,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/7_DP_vs_P_loss.csv")
DP_loss <- DP_loss[,c(1,4,5,10)]

#merge
r <- rbind(DP_gain,DP_loss)
Dat<- r

colnames(Dat)[1] <- "Motif"
colnames(Dat)[2] <- "logP"
Dat <- as.data.frame(Dat)
Dat <- separate(Dat,Motif,into=c("Motif","Motif_info"),sep="/")

tf <- c(
  #DP LOSS
  "SpiB(ETS)","PU.1:IRF8(ETS:IRF)","PU.1(ETS)",
  "ELF5(ETS)","IRF8(IRF)","CTCF(Zf)",
  "BORIS(Zf)","ELF3(ETS)","BIM1(bHLH)",
  "Elf4(ETS)","PU.1-IRF(ETS:IRF)","IRF3(IRF)",
  #DP GAIN
  "Cux2(Homeobox)","TEAD2(TEA)","TEAD(TEA)","TEAD1(TEAD)",
  "NF1(CTF)","RUNX1(Runt)","RUNX2(Runt)","Fos(bZIP)","JunB(bZIP)",
  "Atf3(bZIP)","BATF(bZIP)","Fosl2(bZIP)")

Dat$threshold <- 'Other'
Dat[Dat$`q-value (Benjamini)` < 0.05 & Dat$FC > 0,'threshold'] = 'Gain Peak Enriched Motif'
Dat[Dat$`q-value (Benjamini)` < 0.05 & Dat$FC > 0 & Dat$Motif %in% tf,'threshold'] = 'Gain Labeled'
Dat[Dat$`q-value (Benjamini)` < 0.05 & Dat$FC < 0,'threshold'] = 'Loss Peak Enriched Motif'
Dat[Dat$`q-value (Benjamini)` < 0.05 & Dat$FC < 0 & Dat$Motif %in% tf,'threshold'] = 'Loss Labeled'
table(Dat$threshold)

Dat$threshold <- factor(Dat$threshold,levels = c('Gain Peak Enriched Motif','Gain Labeled','Loss Peak Enriched Motif',
                                                 'Loss Labeled','Other'))

data_text <- Dat[Dat$Motif %in% tf & Dat$`q-value (Benjamini)` < 0.05 ,]
data_text <- separate(data_text,col = Motif,into = c("Motif","Sup"),sep = "[(]")

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/7_DP_vs_P_Motif.pdf",
    width = 7,height = 5)
ggplot(Dat,aes(x=FC,y=-logP,color=threshold))+
  geom_point()+
  scale_color_manual(values=c(  "#EEBBBB","#CD3333","#AAAAD4", "#000080","#808080"))+
  geom_text_repel(
    data = data_text ,
    aes(label = Motif),
    size = 4,
    segment.color = "black", show.legend = FALSE,
    max.overlaps=40)+
  theme_bw()+
  theme(
    legend.title = element_blank()
  )+
  ylab('-log P value')+
  xlab('Fold change')
dev.off()

fwrite(Dat,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/7_DP_vs_P_Motif.csv")

####------19.2 PD1 vs Ctrl  -------------
PD1_gain <- fread("PD1_Ctrl_gain/knownResults.txt")
PD1_gain$FC <- (as.numeric(gsub("%","",PD1_gain$`% of Target Sequences with Motif`)))  / (as.numeric(gsub("%","",PD1_gain$`% of Background Sequences with Motif` ))) 
summary(PD1_gain$FC)

colnames(PD1_gain) <- c("Motif Name","Consensus","P-value","Log P-value","q-value (Benjamini)",
                        "# of Target Sequences with Motif",
                        "% of Target Sequences with Motif",
                        "# of Background Sequences with Motif",
                        "% of Background Sequences with Motif",
                        "FC")
fwrite(PD1_gain,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/7_P_vs_C_gain.csv")
PD1_gain <- PD1_gain[,c(1,4,5,10)]

PD1_loss <- fread("PD1_Ctrl_loss/knownResults.txt")
PD1_loss$FC <- (as.numeric(gsub("%","",PD1_loss$`% of Target Sequences with Motif`)))  / (as.numeric(gsub("%","",PD1_loss$`% of Background Sequences with Motif` ))) 
PD1_loss$FC <- -PD1_loss$FC

colnames(PD1_loss) <- c("Motif Name","Consensus","P-value","Log P-value","q-value (Benjamini)",
                        "# of Target Sequences with Motif",
                        "% of Target Sequences with Motif",
                        "# of Background Sequences with Motif",
                        "% of Background Sequences with Motif",
                        "FC")
fwrite(PD1_loss,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/7_P_vs_C_loss.csv")
PD1_loss <- PD1_loss[,c(1,4,5,10)]

#merge
r <- rbind(PD1_gain,PD1_loss)
Dat<- r

colnames(Dat)[1] <- "Motif"
colnames(Dat)[2] <- "logP"

Dat <- as.data.frame(Dat)

Dat <- tidyr::separate(Dat,col=Motif,into=c("Motif","Motif_info"),sep="\\/")

fc = 0
Dat$threshold = factor(ifelse(Dat$`q-value (Benjamini)` < 0.05 & abs(Dat$FC) >= fc, ifelse(Dat$FC>=fc ,'Up','Down'),'Other'),levels=c('Up','Down','Other'))

tf <- c(
  #PD1 LOSS
  "NF1(CTF)","TEAD(TEA)","CTCF(Zf)","TEAD2(TEA)",
  "BORIS(Zf)","TEAD1(TEAD)","Jun-AP1(bZIP)","Fosl2(bZIP)",
  "Fra2(bZIP)","JunB(bZIP)","Fos(bZIP)","Atf3(bZIP)","BATF(bZIP)",
  "NFAT:AP1(RHD,bZIP)",
  #PD1 GAIN
  "SpiB(ETS)","PU.1:IRF8(ETS:IRF)","PU.1(ETS)",
  "IRF8(IRF)","IRF1(IRF)","Elf4(ETS)","ELF3(ETS)",
  "ETS(ETS)","IRF3(IRF)","IRF2(IRF)","Nur77(NR)")


Dat$threshold <- 'Other'
Dat[Dat$`q-value (Benjamini)` < 0.05 & Dat$FC > 0,'threshold'] = 'Gain Peak Enriched Motif'
Dat[Dat$`q-value (Benjamini)` < 0.05 & Dat$FC > 0 & Dat$Motif %in% tf,'threshold'] = 'Gain Labeled'
Dat[Dat$`q-value (Benjamini)` < 0.05 & Dat$FC < 0,'threshold'] = 'Loss Peak Enriched Motif'
Dat[Dat$`q-value (Benjamini)` < 0.05 & Dat$FC < 0 & Dat$Motif %in% tf,'threshold'] = 'Loss Labeled'
table(Dat$threshold)

Dat$threshold <- factor(Dat$threshold,levels = c('Gain Peak Enriched Motif','Gain Labeled','Loss Peak Enriched Motif',
                                                 'Loss Labeled','Other'))

data_text <- Dat[Dat$Motif %in% tf & Dat$`q-value (Benjamini)` < 0.05 ,]
data_text <- separate(data_text,col = Motif,into = c("Motif","Sup"),sep = "[(]")

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/7_P_vs_C_Motif.pdf",
    width = 7,height = 5)
ggplot(Dat,aes(x=FC,y=-logP,color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#EEBBBB","#CD3333","#AAAAD4", "#000080","#808080"))+
  geom_text_repel(
    data = data_text,
    aes(label = Motif),
    size = 4,
    segment.color = "black", show.legend = FALSE,
    max.overlaps=40)+
  theme_bw()+
  theme(
    legend.title = element_blank()
  )+
  ylab('-log P value')+
  xlab('Fold change')
dev.off()

fwrite(Dat,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/7_P_vs_C_Motif.csv")

####----------19.3 ATAC-seq heatmap correlation----------
# calculated from deeptools
# see in ATAC-seq_data_processing 9 RPKM nomolized 
d <- as.data.frame(fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/ATAC_result/6_deeptools_TSS/outRawCounts_RPKM_bw.tab"))
colnames(d) <- c("chr","start","end","Ctrl1","Ctrl2","P1","P2","DP1","DP2")

d <- d[,4:9]
t <- cor(d,method = "pearson")
pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0620/10.ATAC.correlation.pdf",width = 3.5,height = 3)
pheatmap(t,color = rdwhbu(100))
dev.off()

####----------19.4 merge all sample peaks----
library(dplyr)
library(tidyverse)
setwd("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/ATAC_result/7_peak_new/1_peak")
peak.dt <- fread("all_1_overlap.narrowPeak")
colnames(peak.dt)=c("chr","start","end")

peak.dt2 <- fread("Con_overlap_1_peak")
peak.dt2 <- peak.dt2[,1:3]
colnames(peak.dt2)=c("chr2","start2","end2")

chroms.vector <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                   "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")

all.peaks.dt.list <- lapply(chroms.vector, FUN=function(temp.chrom){
  #temp.chrom="chrY"
  temp.peak.dt <- peak.dt[chr==temp.chrom]
  temp.peak2.dt <- peak.dt2[chr2==temp.chrom]
  aa=temp.peak.dt[, {
    index <- .I
    if ((index+1) %% 1000 == 0){
      cat(date(), " : processing index", index, "\n")
    }
    nrow(temp.peak2.dt[start <= end2 & start2 <=end])
  },  list(chr, start,end)]
  
})
result <- rbindlist(all.peaks.dt.list)
table(result$V1)
result[which(result$V1 >1),'V1'] <- 1
colnames(result)[4] <- "Ctrl"

## PD1
peak.dt <- result

peak.dt2 <- fread("PD1_overlap_1_peak")
peak.dt2 <- peak.dt2[,1:3]
colnames(peak.dt2)=c("chr2","start2","end2")

chroms.vector <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                   "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")

summary(peak.dt)
class(peak.dt$start)
all.peaks.dt.list <- lapply(chroms.vector, FUN=function(temp.chrom){
  #temp.chrom="chrY"
  temp.peak.dt <- peak.dt[chr==temp.chrom]
  temp.peak2.dt <- peak.dt2[chr2==temp.chrom]
  aa=temp.peak.dt[, {
    index <- .I
    if ((index+1) %% 1000 == 0){
      cat(date(), " : processing index", index, "\n")
    }
    nrow(temp.peak2.dt[start <= end2 & start2 <=end])
  },  list(chr, start,end,Ctrl)]
  
})
result <- rbindlist(all.peaks.dt.list)
table(result$V1)
result[which(result$V1 >1),'V1'] <- 1
colnames(result)[5] <- "PD1"

##DP
peak.dt <- result

peak.dt2 <- fread("DP_overlap_1_peak")
peak.dt2 <- peak.dt2[,1:3]
colnames(peak.dt2)=c("chr2","start2","end2")

chroms.vector <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                   "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")

all.peaks.dt.list <- lapply(chroms.vector, FUN=function(temp.chrom){
  #temp.chrom="chrY"
  temp.peak.dt <- peak.dt[chr==temp.chrom]
  temp.peak2.dt <- peak.dt2[chr2==temp.chrom]
  aa=temp.peak.dt[, {
    index <- .I
    if ((index+1) %% 1000 == 0){
      cat(date(), " : processing index", index, "\n")
    }
    nrow(temp.peak2.dt[start <= end2 & start2 <=end])
  },  list(chr, start,end,Ctrl,PD1)]
  
})
result <- rbindlist(all.peaks.dt.list)
table(result$V1)
result[which(result$V1 >1),'V1'] <- 1
colnames(result)[6] <- "DP"

fwrite(result,"/media/liyaru/LYR/301project/2_MOUSE_PD1_DP/ATAC_result/7_peak_new/2_peak_matrix/all_peak.csv")

###-----------19.5 Venn for peaks------------
library(VennDiagram)
peak.dt <- fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/ATAC_result/7_peak_new/2_peak_matrix/all_peak.csv")
head(peak.dt)
peak.dt <- unite(peak.dt,peak,chr,start,end)

Ctrl <- peak.dt [which(peak.dt$Ctrl==1),]
Ctrl <- Ctrl$peak

PD1 <- peak.dt [which(peak.dt$PD1==1),]
PD1 <- PD1$peak

DP <- peak.dt [which(peak.dt$DP==1),]
DP<- DP$peak

venn_list <- list(Ctrl,PD1,DP)
names(venn_list) <- c("Ctrl","PD1","DP")

p=venn.diagram(
  x = venn_list,
  filename = NULL,
  col = "black",
  fill=color_4_sample_1[c(1,3,4)],
  alpha = 0.5
)
pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0620/24.ATAC.venn.pdf",width = 4,height = 4)
grid.draw(p)
dev.off()

####---------20 Fig6F FigS10-G------------
####---------20.1 SCENIC------------------
# calculated by pySCENIC
library(SCopeLoomR)
library(SCENIC)
library(ComplexHeatmap)

setwd("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/result/6_SCENIC/")
pyScenicDir <- "/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/result/6_SCENIC"
pyScenicLoomFile <- file.path(pyScenicDir, "CD8.SCENIC.loom")


loom <- open_loom(pyScenicLoomFile, mode="r")
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
dim(regulons_incidMat)

regulons <- regulonsToGeneLists(regulons_incidMat)
names(regulons)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)

CD8@meta.data$cluster_sample <- paste0(CD8@meta.data$RNA_snn_res.0.5,"_",CD8@meta.data$orig.ident)

# by cluster_sample
cellClusters <- as.data.frame(CD8@meta.data$cluster_sample)
rownames(cellClusters) <- rownames(CD8@meta.data)
colnames(cellClusters) <- "celltype_sample"

selectedResolution <- "celltype_sample" # select resolution
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,selectedResolution])

###-------- 20.2 TF heatmap-----------
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))


options(repr.plot.width=8, repr.plot.height=10)
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size

rownames(regulonActivity_byCellType)

TF <- c("Jun(+)","Junb(+)","Jund(+)",
        "Fos(+)","Fosb(+)","Fosl1(+)","Fosl2(+)",
        "Cux1(+)",
        "Runx1(+)","Runx2(+)",
        "Nfkb1(+)","Nfkb2(+)",
        "Ets1(+)","Hivep1(+)",
        "Xbp1(+)" ,"Irf4(+)",
        "Batf(+)",
        "Eomes(+)", "Tbx21(+)","Bach2(+)",
        "Lef1(+)")

samples <- c( "0_Con","0_DAC","0_PD1","0_DP",
              "1_Con","1_DAC","1_PD1","1_DP",
              "4_Con","4_DAC","4_PD1","4_DP")

pheatmap(regulonActivity_byCellType_Scaled[TF,samples],
         cluster_cols = F,
         color = rdwhbu(100))

t <- regulonActivity_byCellType_Scaled
fwrite(as.data.frame(t),"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/19.SCENIC.activity.scaled.csv",
       row.names = T)

t <- regulonActivity_byCellType
fwrite(as.data.frame(t),"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/19.SCENIC.activity.csv",
       row.names = T)

###-------- 20.3 regulatory network-----------
loom <- open_loom(pyScenicLoomFile, mode="r")
# Read information from loom file:
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
SCENIC_TF <- rownames(regulons_incidMat)

rownames(regulons_incidMat)
tf <- c("Jund(+)")

t <- regulons_incidMat[tf,]
t <- as.data.frame(t)
t$gene <- rownames(t)

tt <- regulons$`Jund(+)`

t <- as.data.frame(tt)
colnames(t) <- "gene"

a <- fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/result/6_SCENIC/adj.CD8.tsv")
head(a)

aa <- a[a$TF %in% c("Jund") & a$importance > 30,] #final use

tt <- intersect(aa$target,t$gene)

target <- aa[aa$target %in% tt,]

# DP vs P FC
a <- as.data.frame(fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/C0.DP_PD1_0.05.csv"))  
a <- a[order(-a$avg_log2FC),]
colnames(a)[1] <- "gene" 
colnames(target)[2] <- "gene"

m <- merge(target,a,by="gene",all.x=T)

m <- m[,c(1:5)]
colnames(m) <- c("gene","TF","importance","p_val (DP vs P in progenitor Tex)","avg_log2FC (DP vs P in progenitor Tex)")

# P vs Ctrl FC
a <- as.data.frame(fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table/DEG/C0.PD1_Con.csv"))  
a <- a[order(a$avg_log2FC),]
colnames(a)[1] <- "gene" 

m <- merge(m,a[,c(1:3)],by="gene",all.x=T)
colnames(m)[6:7] <- c("p_val (P vs C in progenitor Tex)","avg_log2FC (P vs C in progenitor Tex)")

#fwrite(m,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table/network/JUND_target_importance30.csv")

###---------21 Fig7A FigS11A-B----------------
#####-------21.1 JunD exp in TSNE-------------
pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF/16.Jund_featureplot.pdf",
    height = 3,width = 12)
FeaturePlot(CD8,features = "Jund",split.by = "orig.ident",pt.size = 1.2, ncol=2)& 
  scale_color_gradientn(colors = c(rdwhbu(100)[1:45],rdwhbu(100)[55:100]))
dev.off()


#####-------21.2 scImpute-------------
library(scImpute)
library(data.table)
setwd("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/result/12_impute")
count <- fread("./count_matrix.csv")
cells <- colnames(count)[-1]
head(cells)
meta <- CD8@meta.data
labels <- as.character(meta[cells,"RNA_snn_res.0.5"])
head(labels)

scimpute(# full path to raw count matrix
  count_path = "./count_matrix.csv", 
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = "./",           # full path to output directory
  labeled = T,          # cell type labels not available
  labels = labels,
  drop_thre = 0.5,          # threshold set on dropout probability
  #Kcluster = 2,             # 2 cell subpopulations
  ncores = 10)              # number of cores used in parallel computation

#impute result
a <- fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/result/12_impute/scimpute_count.csv")
a <- as.data.frame(a)
row.names(a) <- a$V1
a <- a[,-1]
colnames(a) <- gsub("[.]","-",colnames(a))
a[1:5,1:5]

object <- CreateSeuratObject(counts = a,min.cells = 0, min.features = 0)
object@meta.data
colnames(object)

m <- CD8@meta.data
m$cell <- rownames(m)

object <- AddMetaData(object,m)

mm <- object@meta.data


object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(object)
object <- ScaleData(object,features = all.genes)


Idents(object) <- object$RNA_snn_res.0.5
object2 <- subset(object,idents=c("0","1","2","4"))

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/6_Impute_expression_Jund.pdf",width = 5,height = 3)
VlnPlot(object2,
        #features = c("Junb","Jund","Jun"),
        features = "Jund",
        split.by = "orig.ident",
        group.by = "RNA_snn_res.0.5",
        pt.size = 0,
        cols = color_4_sample,
        ncol = 1,
        #slot = "count"
        slot="data"
        #slot = "scale.data"
)
dev.off()

CD8_2 <- subset(CD8,ident=c("0","1","2","4"))

VlnPlot(CD8_2,features = "Jund",
        split.by = "orig.ident",
        group.by = "RNA_snn_res.0.5",
        pt.size = 0,
        cols = color_4_sample,
        ncol = 1,
        slot = "data")

####-------21.3 JunD exp in public data--------------
# data from GSE179994 
# processed expression data from http://nsclcpd1.cancer-pku.cn/

a <- fread("/media/liyaru/LYR1/301project/public_data/nature cancer 21/JUN.exp.csv")
a <- as.data.frame(a)
summary(a$JUND)

a$Treatment <- factor(a$Treatment,levels=c("pre","post"))
a$Sub_Cluster <- factor(a$Sub_Cluster,levels = c("Non-exhausted","Terminal Tex","Proliferative"))


stat_wilcox <- wilcox_test(group_by(a, Sub_Cluster), JUND ~ Treatment)  
stat_wilcox <- add_significance(stat_wilcox, 'p')  
stat_wilcox.test <-  add_xy_position(stat_wilcox, x = 'Sub_Cluster')

p <- ggboxplot(a,x="Sub_Cluster",y="JUND" ,
               fill='Treatment',
               #color = "Treatment",
               outlier.shape=NA,
               #add="jitter"
)+scale_fill_manual(values = color_4_sample_1[c(1,3)])


p+ stat_pvalue_manual(stat_wilcox.test, 
                      #label = 'p', 
                      label = 'p.signif',
                      tip.length = 0.005,
                      hide.ns=T)  

stat_wilcox <- wilcox_test(group_by(a, Sub_Cluster), JUN ~ Treatment)  
stat_wilcox <- add_significance(stat_wilcox, 'p')  
stat_wilcox.test <-  add_xy_position(stat_wilcox, x = 'Sub_Cluster')

p <- ggboxplot(a,x="Sub_Cluster",y="JUN" ,
               fill='Treatment',
               #color = "Treatment",
               outlier.shape=NA,
               #add="jitter"
)+scale_fill_manual(values = color_4_sample_1[c(1,3)])

p+ stat_pvalue_manual(stat_wilcox.test, 
                      #label = 'p', 
                      label = 'p.signif',
                      tip.length = 0.005,
                      hide.ns=T)  


###--------22 Fig7B-------------------------- 
###--------22.1 activation score--------------- 
library(ggpubr)
library(rstatix)

# T cell activation
# GO:0042110
t <- read_excel("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GO/GO_term_summary_20220103_120331_T_cell_activation.xlsx")
tt <- unique(t$Symbol)

CD8 <- AddModuleScore(CD8,features = list(tt),name = "T cell activation")
colnames(CD8[[]])


t <- CD8@meta.data
colnames(t)


p <- ggboxplot(t, x="RNA_snn_res.0.5",y="T.cell.activation1",
               fill = 'orig.ident',
               outlier.shape=NA,
               palette =  color_4_sample_1,)
p

stat_wilcox <- wilcox_test(group_by(t, RNA_snn_res.0.5), T.cell.activation1 ~ orig.ident)  
stat_wilcox <- add_significance(stat_wilcox, 'p')  
stat_wilcox.test <-  add_xy_position(stat_wilcox, x = 'RNA_snn_res.0.5')

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0620/4.activation.score.pdf",width = 5,height = 5)
p+ stat_pvalue_manual(stat_wilcox.test, 
                      #label = 'p', 
                      label = 'p.signif',
                      tip.length = 0.005,
                      hide.ns=T)  
dev.off()

#check
tt <- compare_means(T.cell.activation1 ~ orig.ident, 
                    data = t, 
                    group.by = "RNA_snn_res.0.5" ,
                    paired = F)


###--------23 Fig7C-D----------------------------
# GSEA of ACT assay
###--------23.1 GSEA----------------------------- 
library(clusterProfiler)
colnames(merge_inter@meta.data)
DimPlot(merge_inter,group.by = "cluster_sample")
unique(merge_inter$cluster_sample)

DefaultAssay(merge_inter) <- "RNA"
markers <- FindMarkers(merge_inter,
                       group.by = "cluster_sample",
                       ident.1 = "Proliferating_DP",
                       ident.2 = "Proliferating_PD1",
                       #test.use = "MAST"
                       logfc.threshold=0
)
fwrite(markers,
       "/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/Pro.DP_PD1.csv",
       row.names =T )

a <- fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/Pro.DP_PD1.csv")
glist <- a$avg_log2FC
names(glist) <- a$V1
glist <- sort(glist,decreasing = T)
names(glist) <- toupper(names(glist))

geneset <- read.gmt("/home/liyaru/public_Data/GSEA.gmt/c7.all.v7.5.1.symbols.gmt")

#p.adjust.methods
gsea <- GSEA(glist, TERM2GENE=geneset, verbose=FALSE,pvalueCutoff=1,seed=618)
t <- gsea@result
fwrite(t,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/15.Proliferating.GSEA.c7.immunologic signature gene sets.csv")

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/17.selected_GSEA_Pro.2.pdf",width = 5,height = 4)
enrichplot::gseaplot2(gsea,geneSetID = "GSE41867_DAY8_EFFECTOR_VS_DAY30_EXHAUSTED_CD8_TCELL_LCMV_CLONE13_UP")
enrichplot::gseaplot2(gsea,geneSetID = "GSE41867_DAY6_EFFECTOR_VS_DAY30_MEMORY_CD8_TCELL_LCMV_ARMSTRONG_UP")
enrichplot::gseaplot2(gsea,geneSetID = "GSE20754_WT_VS_TCF1_KO_MEMORY_CD8_TCELL_DN")
dev.off()

####-----23.2 GSEA Barplot----------------
t <- readxl::read_excel("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/17.selected_GSEA_Pro.xlsx")
t$LogP <- log10(t$pvalue)
summary(t)
t <- t[order(-t$NES),]

t$Description <- factor(t$Description,levels=t$Description)

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/17.selected_GSEA_Pro.pdf",
    width = 12,height = 4)

ggplot(t, aes(x=NES, y=Description, fill=-LogP)) + 
  geom_bar(stat='identity') + 
  scale_fill_gradientn(colors = c(rdwhbu(100)[1:45],rdwhbu(100)[55:100]))+
  theme_classic()+ ylab(NULL)
dev.off()


###--------24 FigS12A-B --------------
# GSEA of in  vivo Tex 
###--------24.1 GSEA------------------
library(clusterProfiler)
library(data.table)
library(org.Mm.eg.db)
library(openxlsx)
library(GSEABase) 
set.seed(618)


#C0
a <- fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/C0.DP_PD1.csv")
keytypes(org.Mm.eg.db)

glist <- a$avg_log2FC
names(glist) <- a$V1
glist <- sort(glist,decreasing = T)
names(glist) <- toupper(names(glist))

# c7: immunologic signature gene sets
geneset <- read.gmt("/home/liyaru/public_Data/GSEA.gmt/c7.all.v7.5.1.symbols.gmt")

#p.adjust.methods
gsea <- GSEA(glist, TERM2GENE=geneset, verbose=FALSE,pvalueCutoff=1,seed=618)
t <- gsea@result

fwrite(t,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/15.C0.GSEA.c7.immunologic signature gene sets.csv")


#C1
DefaultAssay(CD8) <- "RNA"
CD8.markers <- FindMarkers(CD8,
                           group.by = "cluster_sample",
                           ident.1 = "1_DP",
                           ident.2 ="1_PD1",
                           #test.use = "MAST"
                           logfc.threshold=0
)
fwrite(CD8.markers,
       "/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/C1.DP_PD1.csv",
       row.names =T )

a <- fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/C1.DP_PD1.csv")
keytypes(org.Mm.eg.db)

glist <- a$avg_log2FC
names(glist) <- a$V1
glist <- sort(glist,decreasing = T)
names(glist) <- toupper(names(glist))


geneset <- read.gmt("/home/liyaru/public_Data/GSEA.gmt/c7.all.v7.5.1.symbols.gmt")

#p.adjust.methods
gsea <- GSEA(glist, TERM2GENE=geneset, verbose=FALSE,pvalueCutoff=1,seed=618)
t <- gsea@result
fwrite(t,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/15.C1.GSEA.c7.immunologic signature gene sets.csv")


#C4
DefaultAssay(CD8) <- "RNA"
CD8.markers <- FindMarkers(CD8,
                           group.by = "cluster_sample",
                           ident.1 = "4_DP",
                           ident.2 ="4_PD1",
                           #test.use = "MAST"
                           logfc.threshold=0
)
fwrite(CD8.markers,
       "/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/C4.DP_PD1.csv",
       row.names =T )



a <- fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/C4.DP_PD1.csv")
keytypes(org.Mm.eg.db)

glist <- a$avg_log2FC
names(glist) <- a$V1
glist <- sort(glist,decreasing = T)
names(glist) <- toupper(names(glist))

geneset <- read.gmt("/home/liyaru/public_Data/GSEA.gmt/c7.all.v7.5.1.symbols.gmt")

#p.adjust.methods
gsea <- GSEA(glist, TERM2GENE=geneset, verbose=FALSE,pvalueCutoff=1,seed=618)
t <- gsea@result
fwrite(t,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/15.C4.GSEA.c7.immunologic signature gene sets.csv")

pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/17.selected_GSEA_C4_2.pdf",width = 5,height = 4)
enrichplot::gseaplot2(gsea,geneSetID = "GSE41867_DAY15_EFFECTOR_VS_DAY30_EXHAUSTED_CD8_TCELL_LCMV_CLONE13_UP")
dev.off()


#C2
DefaultAssay(CD8) <- "RNA"
CD8.markers <- FindMarkers(CD8,
                           group.by = "cluster_sample",
                           ident.1 = "2_DP",
                           ident.2 ="2_PD1",
                           #test.use = "MAST"
                           logfc.threshold=0
)
fwrite(CD8.markers,
       "/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/C2.DP_PD1.csv",
       row.names =T )



a <- fread("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GSEA/C2.DP_PD1.csv")
keytypes(org.Mm.eg.db)

glist <- a$avg_log2FC
names(glist) <- a$V1
glist <- sort(glist,decreasing = T)
names(glist) <- toupper(names(glist))

geneset <- read.gmt("/home/liyaru/public_Data/GSEA.gmt/c7.all.v7.5.1.symbols.gmt")

#p.adjust.methods
gsea <- GSEA(glist, TERM2GENE=geneset, verbose=FALSE,pvalueCutoff=1,seed=618)
t <- gsea@result
fwrite(t,"/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/15.C2.GSEA.c7.immunologic signature gene sets.csv")
enrichplot::gseaplot2(gsea,geneSetID = 1,pvalue_table=T)


###---------24.2 Barplot------------------
t <- readxl::read_excel("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/16.selected_GSEA.xlsx")
t$LogP <- log10(t$pvalue)
summary(t)

t <- t[t$Cluster=="0",]
t <- t[order(-t$NES),]
t$Description <- factor(t$Description,levels=t$Description)
pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/17.selected_GSEA_C0.pdf",
    width = 10,height = 2)
ggplot(t, aes(x=NES, y=Description, fill=-LogP)) + 
  geom_bar(stat='identity') + 
  scale_fill_gradientn(colors = c(rdwhbu(100)[1:45],rdwhbu(100)[55:100]))+
  theme_classic()+ylab(NULL)
dev.off()

t <- readxl::read_excel("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/16.selected_GSEA.xlsx")
t$LogP <- log10(t$pvalue)
summary(t)
t <- t[t$Cluster=="1",]
t <- t[order(-t$NES),]
t$Description <- factor(t$Description,levels=t$Description)
pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/17.selected_GSEA_C1.pdf",
    width = 10,height = 3)
ggplot(t, aes(x=NES, y=Description, fill=-LogP)) + 
  geom_bar(stat='identity') + 
  scale_fill_gradientn(colors = c(rdwhbu(100)[1:45],rdwhbu(100)[55:100]))+
  theme_classic()+ylab(NULL)
dev.off()

t <- readxl::read_excel("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/16.selected_GSEA.xlsx")
t$LogP <- log10(t$pvalue)
summary(t)
t <- t[t$Cluster=="2",]
t <- t[order(-t$NES),]
t$Description <- factor(t$Description,levels=t$Description)
pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/17.selected_GSEA_C2.pdf",
    width = 10,height = 3)
ggplot(t, aes(x=NES, y=Description, fill=-LogP)) + 
  geom_bar(stat='identity') + 
  scale_fill_gradientn(colors = c(rdwhbu(100)[1:45],rdwhbu(100)[55:100]))+
  theme_classic()+ylab(NULL)
dev.off()

t <- readxl::read_excel("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/16.selected_GSEA.xlsx")
t$LogP <- log10(t$pvalue)
summary(t)
t <- t[t$Cluster=="4",]
t <- t[order(-t$NES),]
t$Description <- factor(t$Description,levels=t$Description)
pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0531/17.selected_GSEA_C4.pdf",
    width = 10,height = 1)
ggplot(t, aes(x=NES, y=Description, fill=-LogP)) + 
  geom_bar(stat='identity') + 
  scale_fill_gradientn(colors = c(rdwhbu(100)[1:45],rdwhbu(100)[55:100]))+
  theme_classic()+ylab(NULL)
dev.off()

####------25 Fig7E & FigS12E -----------------------
# Fig7E in vivo 
####------25.1 Jund exp level class------
t <- as.data.frame(GetAssayData(CD8,assay = "RNA"))

tt <-as.data.frame(t(t["Jund",]))
summary(tt)
tt$cell <- rownames(tt)

quantile(tt$Jund)
#cut by quantile
tt <- tt[order(tt$Jund),]
tt$Jund_class <- as.character(as.numeric(cut(1:nrow(tt),breaks = 4)))
table(tt$Jund_class)

#cut by median
# ttt <- median(tt$Jund)
# tt[tt$Jund > ttt,'Jund_class'] <- "JunD High"
# tt[tt$Jund < ttt,'Jund_class'] <- "JunD Low"
# table(tt$Jund_class)

CD8 <- AddMetaData(CD8,tt)

#CD8$Jund_class <- factor(CD8$Jund_class,levels = c("JunD Low","JunD High"))

####-----25.2 Jund & activation---------
t <- CD8@meta.data
pdf("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/PDF_0711/1.Jund.Activation.Proliferation.pdf")
p <- ggboxplot(t, x="RNA_snn_res.0.5",y="T.cell.activation1",
               fill = 'Jund_class',
               outlier.shape=NA,
               palette =  c("#38389C","#AAAAD4","#EEBBBB","#D86060")
)

stat_wilcox <- wilcox_test(group_by(t, RNA_snn_res.0.5), T.cell.activation1 ~ Jund_class)  
stat_wilcox <- add_significance(stat_wilcox, 'p')  
stat_wilcox.test <-  add_xy_position(stat_wilcox, x = 'RNA_snn_res.0.5')

p+ stat_pvalue_manual(stat_wilcox.test, 
                      #label = 'p', 
                      label = 'p.signif',
                      tip.length = 0.005,
                      hide.ns=T)

#ACT assay #FigS12E
####--------25.3 Jund level class in ACT assay------------------
t <- as.data.frame(GetAssayData(merge_inter,assay = "RNA"))

tt <-as.data.frame(t(t["Jund",]))
summary(tt)
tt$cell <- rownames(tt)

quantile(tt$Jund)
#cut by quantile
tt <- tt[order(tt$Jund),]
tt$Jund_class <- as.character(as.numeric(cut(1:nrow(tt),breaks = 4)))
table(tt$Jund_class)

#cut by median
# ttt <- median(tt$Jund)
# tt[tt$Jund > ttt,'Jund_class'] <- "JunD High"
# tt[tt$Jund < ttt,'Jund_class'] <- "JunD Low"
# table(tt$Jund_class)

merge_inter <- AddMetaData(merge_inter,tt)

###------25.4 Jund & activation score in ACT assay ----------
DefaultAssay(merge_inter) <- "RNA"
t <- read_excel("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/table_old/GO/GO_term_summary_20220103_120331_T_cell_activation.xlsx")
tt <- unique(t$Symbol)
merge_inter <- AddModuleScore(merge_inter,features = list(tt),name = "T cell activation")
t <- merge_inter@meta.data
p <- ggboxplot(t, x="celltype",
               y="T.cell.activation1",
               fill = 'Jund_class',
               outlier.shape=NA,
               palette =  c("#38389C","#AAAAD4","#EEBBBB","#D86060")
)
stat_wilcox <- wilcox_test(group_by(t, celltype), T.cell.activation1 ~ Jund_class)  
stat_wilcox <- add_significance(stat_wilcox, 'p')  
stat_wilcox.test <-  add_xy_position(stat_wilcox, x = 'celltype')

p+ stat_pvalue_manual(stat_wilcox.test, 
                      #label = 'p', 
                      label = 'p.signif',
                      tip.length = 0.005,
                      hide.ns=T)  


