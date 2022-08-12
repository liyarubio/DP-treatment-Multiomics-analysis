library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggsci)
library(monocle3)
library(data.table)
library(readxl)
library(data.table)

setwd("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP")
getwd()

####----------read from 10X---------------------
# The raw data and matrix from 10X are available from  GEO:GSE209965 

####------------Build Seurat-------------------------------
#read data
path_all="./cellranger1"
files=list.files(path_all)
class(files)
paths=NULL
mylist=list()
class(mylist)

for (i in c(1:length(files))){
  file=paste0("./cellranger1/",files[i],"/outs/filtered_feature_bc_matrix")
  paths[i]=file
  
  data<-Read10X(data.dir = file)
  # default min.cell = 3
  object <- CreateSeuratObject(counts = data,project = files[i],min.cells = 5, min.features = 200)
  
  # object[["sample"]] <- sample[[1]][1]
  # object[["sample15"]] <- paste(sample[[1]][1],sample[[1]][2],sep=" ")
  # object[["type"]] <- substring(sample[[1]][1],1,2)
  
  mylist<- c(mylist,object)
}
names(mylist)<-files
paths

mylist[[1]]@meta.data

####-------------merge data-------------------------------------
merge<- merge(mylist[[1]],
              y=mylist[2:length(files)],
              add.cell.ids=files,
              project="DP")
head(merge@meta.data)
#merge@meta.data$sample <- sub("-5-lib","",merge@meta.data$orig.ident)
merge@meta.data$sample <- merge@meta.data$orig.ident
unique(merge$sample)
merge@meta.data$sample <- factor(merge@meta.data$sample,
                                 levels = c("Con","DAC","PD1","DP"))

#vlnplot
merge[["percent.mt"]] <- PercentageFeatureSet(merge,pattern = "^mt")
VlnPlot(merge,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3,pt.size=0,
        group.by = "sample")
table(merge@meta.data$sample)

orig_merge=c(dim(merge),table(merge$sample))
orig_merge

names(orig_merge)[1] <- "gene"
names(orig_merge)[2] <- "cell"

df = as.data.frame(matrix(nrow=0,ncol=(length(files)+2)))

df[1,]<- orig_merge
colnames(df)<-names(orig_merge)

####------filter-------------------------------------
merge <- subset(merge, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(merge,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3,pt.size=0,
        group.by = "sample")

filter_merge=c(dim(merge),table(merge$sample))
filter_merge

df[2,]<- filter_merge

rownames(df)<-c("before_filter","after_filter")

write.csv(df,file="result/2_seurat/1_merge_cell_number.csv")

#df_t <- t(df)
df_plot<-df[,3:ncol(df)]
df_plot_t <- t(df_plot)

df_long<-melt(df_plot_t,measure.vars = c('before_filter','after_filter'),
              variable.name = 'filter',
              value.name = 'n')
colnames(df_long) <- c ("sample","filter","cell_number")

scales::show_col(pal_npg()(4))
ttt <- colorRampPalette(colors = pal_npg()(10))
scales::show_col(ttt(20))

# check filter cell number
ggplot(df_long,aes(x=sample,y=cell_number))+
  geom_bar(stat = 'identity',aes(fill = filter),position = position_dodge(),width = 0.5) +
  theme_classic()+
  theme(text=element_text(family="Songti SC",size=26,face = "bold"), 
        axis.text.x = element_text(size=20,angle=45,hjust = 1))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = ttt(4))

#ggsave("./FIG/1/1_barplot_cell_number.png",width = 16,height = 9,dpi=350)

# pipeline
merge <- NormalizeData(merge, normalization.method = "LogNormalize", scale.factor = 10000)
merge <- FindVariableFeatures(merge, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merge)
merge <- ScaleData(merge,features = all.genes)
merge <- RunPCA(merge, npcs = 30, verbose = FALSE)
ElbowPlot(merge, ndims = 30)
merge<- FindNeighbors(merge, dims = 1:20)
merge<- FindClusters(merge, resolution = 0.5)
merge <- RunTSNE(merge,dims = 1:20)
DimPlot(merge, reduction = "tsne",label=TRUE, pt.size = 1.2)
DimPlot(merge,group.by = "sample",pt.size = 0.6)

####---------------------integrate data (batch effect)--------------------
mylist2 <- list()

for (i in mylist){
  i[["percent.mt"]] <- PercentageFeatureSet(i,pattern = "^mt")
  i <- subset(i, subset = nFeature_RNA > 1000& nFeature_RNA < 6000 & percent.mt < 5)
  i <- NormalizeData(i, normalization.method = "LogNormalize", scale.factor = 10000)
  i <- FindVariableFeatures(i, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(i)
  i <- ScaleData(i,features = all.genes)
  mylist2 <- c(mylist2,i)
}

anchors <- FindIntegrationAnchors(object.list = mylist2,dims=1:20)
merge_inter <- IntegrateData(anchorset = anchors, dims = 1:20)

merge_inter@assays
merge_inter@meta.data

table(merge_inter$orig.ident)

DefaultAssay(merge_inter) <- "integrated"
merge_inter <- ScaleData(merge_inter, verbose = FALSE)

####-------basic pipeline----------------------------
DefaultAssay(merge_inter) <- "integrated"
merge_inter <- RunPCA(merge_inter, npcs = 30, verbose = FALSE)
ElbowPlot(merge_inter, ndims = 30)
merge_inter <- FindNeighbors(merge_inter, dims = 1:20)
merge_inter <- FindClusters(merge_inter, resolution = c(0.3,0.5,0.7,0.9,1))
merge_inter <- RunTSNE(merge_inter,dims = 1:20)

DimPlot(merge_inter, reduction = "tsne",label=TRUE, pt.size = 0.6)
DimPlot(merge_inter,group.by = "orig.ident", pt.size = 0.6)

#merge_inter@meta.data$sample <- sub("-5-lib","",merge_inter@meta.data$orig.ident)
merge_inter@meta.data$sample <- factor(merge_inter@meta.data$orig.ident,
                                       levels = c("Con","DAC","PD1","DP"))
DimPlot(merge_inter,group.by = "sample", pt.size = 0.6)

colnames(merge_inter@meta.data)
DimPlot(merge_inter,group.by = "integrated_snn_res.0.5",pt.size = 0.6,
        label = T,label.size = 6)
DimPlot(merge_inter,group.by = "integrated_snn_res.1",pt.size = 0.6,
        label = T,label.size = 6)

FeaturePlot(merge_inter,
            features = c("Cd3d","Cd3e","Cd3g",
                         "Cd8a","Cd8b1","Cd4",
                         "Il7r","Nkg7","Gzmb",
                         "Fcgr3","Ptprc","Ncam1","Klrb1c"),
            ncol = 4,pt.size = 1)

###------add TCR information-------------------
# merge TCR file , adjust barcode  
sup <- c("-1_1","-1_2","-1_3","-1_4")
j=1
for (i in c("Con","DAC","DP","PD1")){
  #i="Con"
  tcr <- read.csv(paste0("./1_cellranger_result/TCR/",i,"-TCR/outs/filtered_contig_annotations.csv"))
  tcr$barcode <- gsub("-1", "", tcr$barcode)
  tcr <- tcr[!duplicated(tcr$barcode), ]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  
  clono <- read.csv(paste0("./1_cellranger_result/TCR/",i,"-TCR/outs/clonotypes.csv"))
  tcr <- merge(tcr, clono)
  tcr$sample <- i
  tcr$barcode <- paste0(tcr$barcode,sup[j])
  j=j+1
  
  if (i=="Con"){
    tcr_all <-tcr
  }else{
    tcr_all <- rbind(tcr_all,tcr)
  }
  #TCR <- c(TCR,list(tcr))
}

unique(tcr_all$barcode)
rownames(tcr_all) <- tcr_all$barcode
merge_inter <- AddMetaData(merge_inter,metadata = tcr_all)

DimPlot(merge_inter)

#####------RDS------------------
#saveRDS(merge_inter,"/media/liyaru/LYR/301project/2_MOUSE_PD1_DP/RDS/new/merge.rds")
merge_inter <- readRDS("/media/liyaru/LYR1/301project/2_MOUSE_PD1_DP/RDS/new/merge.rds")

#####------CD8 T cell------------------
Idents(merge_inter) <- merge_inter@meta.data$integrated_snn_res.0.5
CD8  <- subset(merge_inter,idents=c("2","5","7"))
DimPlot(CD8)
CD8@meta.data$orig.ident <- factor(CD8@meta.data$orig.ident,
                                   levels = c("Con","DAC","PD1","DP"))

#remove cell without information
t <-CD8@meta.data
tt <- rownames(t[which(!is.na(t$chain)),])
CD8<- subset(CD8,cells=tt)

DefaultAssay(CD8) <- "RNA"
all.genes <- rownames(CD8)
CD8 <- CD8 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(features = all.genes) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution=c(0.3,0.5,0.8,1.0,1.2)) %>%
  RunTSNE(dims = 1:20) 

DimPlot(CD8,group.by = "sample")
DimPlot(CD8,group.by = "integrated_snn_res.0.5")


DefaultAssay(CD8) <- "integrated"
CD8 <- CD8 %>%
  #NormalizeData() %>%
  #FindVariableFeatures() %>%
  #ScaleData(features = all.genes) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution=c(0.3,0.5,0.8,1.0,1.2)) %>%
  RunTSNE(dims = 1:20) 
DimPlot(CD8,group.by = "sample")
DimPlot(CD8,group.by = "integrated_snn_res.0.5")

####-----filter CD8 cluster 7 8 ------------------
FeaturePlot(CD8,features = c("Cd8a","Cd8b1","Cd4","Cd3d"),
            pt.size = 1)

Idents(CD8) <- CD8@meta.data$RNA_snn_res.0.5
CD8 <- subset(CD8,idents=c("0","1","2","3","4","5","6"))
DimPlot(CD8,
        pt.size = 1,
        label = T,
        label.size = 10)

saveRDS(CD8,"/media/liyaru/LYR/301project/2_MOUSE_PD1_DP/RDS/new/CD8_C0-6.rds")
















