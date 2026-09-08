
rm(list=ls())

####one.创建seurat对象####
library(Seurat)
library(harmony)
library(tidyverse)
library(dplyr)
library(patchwork)
library(tidydr)
library(ggplot2)
library(cowplot)
library(data.table)
setwd("D:\\R language\\5.ADAM17\\NO2.scRNA\\GSE117891")
getwd()

# 读取 TSV 文件
count_matrix <- read.delim("D:\\R language\\5.ADAM17\\NO2.scRNA\\all_6148.umi.count.matrix.tsv", row.names = 1)
# 转换为稀疏矩阵（节省内存）
library(Matrix)
count_matrix <- as(as.matrix(count_matrix), "dgCMatrix")
# 创建 Seurat 对象
pbmc <- CreateSeuratObject(
  counts = count_matrix,  # UMI 计数矩阵
  project = "scRNAseq",   # 项目名称（自定义）
  min.cells = 3,          # 保留在至少 3 个细胞中表达的基因
  min.features = 200)     # 保留至少检测到 200 个基因的细胞

saveRDS(pbmc, file = "scRNAseq_6148_seurat_object.rds")


####two.数据质控####
pbmc <- readRDS("scRNAseq_6148_seurat_object.rds")
#向pbmc新增一列percent.mt数据
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#使用小提琴图可视化QC指标
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0, ncol = 3)
dev.off()
##FeatureScatter通常用于可视化 feature-feature 相关性，
#nCount_RNA 与percent.mt的相关性
plot1<-FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#nCount_RNA与nFeature_RNA的相关性
plot2<-FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#合并两图
plot1 + plot2 
pbmc<- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 10,pt.size=0)
pbmc


####three.标准化####
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
##鉴定高变基因##
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# 查看最高变的10个基因
top10 <- head(VariableFeatures(pbmc), 10)
top10
# 画出不带标签或带标签基因点图
plot1<-VariableFeaturePlot(pbmc)
plot2<-LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
pbmc


####four.细胞分类####
##1.归一化数据：数据缩放
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

##2.PCA降维
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#查看PCA结果
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
#可视化#
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
a<-DimPlot(pbmc, reduction = "pca")
#DimHeatmap()可以方便地探索数据集中异质性的主要来源
#并且可以确定哪些PC维度可以用于下一步的下游分析。细胞和基因根据PCA分数来排序。
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) #1个PC 500个细胞
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) #15个PC

#harmony去批次#
library(harmony)
pbmc_harmony <- RunHarmony(pbmc, group.by.vars = "orig.ident")
pbmc_harmony@reductions[["harmony"]][[1:5,1:5]]
b=DimPlot(pbmc_harmony,reduction = "harmony",group.by = "orig.ident")
#PCA图看到还有一点的批次效应（融合得比较好批次就弱）
b
pca_harmony_integrated <- CombinePlots(list(a,b),ncol=1)
pca_harmony_integrated

##3.确立数据的维度
#JackStraw
pbmc <- JackStraw(pbmc, num.replicate = 100)  #最开始是50
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)    #最开始是dims = 1:15
JackStrawPlot(pbmc, dims = 1:20)
#Elbow plot
ElbowPlot(pbmc)

##4.细胞聚类(Seurat使用KNN算法进行聚类。)
#dims = 1:10 即选取前30个主成分来分类细胞。
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.5)  #最开始是1
#查看前5个细胞的分类ID
head(Idents(pbmc), 5)

##5.非线性降维##
pbmc <- RunUMAP(pbmc, dims = 1:15)
DimPlot(pbmc, reduction = "umap")
# 显示在聚类标签
DimPlot(pbmc, reduction = "umap", label = TRUE)
# 使用TSNE聚类
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")
# 显示在聚类标签
DimPlot(pbmc, reduction = "tsne", label = TRUE)
#保存rds，用于后续分析
saveRDS(pbmc, file = "pbmc_UMAP+TSNE.rds")



####five.寻找mark####
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 1)
pbmc<-JoinLayers(pbmc)
markers <- FindAllMarkers(object = pbmc, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.5,
                          min.pct = 0.25)  #要限制一下在pct里面的一个最低表达
#要注意加这个min.pct=0.5的限制，不然可能会出现假阳性高变基因的情况

saveRDS(markers,file = "FindMarker.rds")
#直接双击load吧

#寻找每个簇的高变基因
library(dplyr)
all.markers = markers %>% 
  dplyr::select(gene, everything()) %>% 
  subset(p_val_adj<0.01)
#将每个cluster lgFC排在前10的marker基因挑选出来
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
View(top10)
#将top gene放在一起
gene_column <- top10[ , "gene", drop = FALSE]
gene_vector <- top10[ , "gene"]
view(gene_vector)


cell_type_markers = c("MOG",#Non-tumor1
                      "KRT5",#LC-meta
                      "MKI67","HYDIN","IGFBPL1","CLU","APQ4","PDGFRA","OLIG2","PTPRZ1", #Glioma
                      "FCN1","FCGR3B","CXCL1","CD163","P2RY12","CD3E","PTPRC") #Non-tumor2"



###ds给的marker
cell_type_markers = c("EGFR", "PDGFRA", "SOX2", "VIM", "TNC", "S100A4", "PTPRZ1", #Malignant
                      "CD68", "AIF1", "CSF1R", "CD163", "MRC1", "ITGAM", "CX3CR1", #Macrophage
                      "CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD8A", "CD4", "GZMK", "GZMB", #T_cell
                      "NCAM1", "KLRB1", "NCR1", "NKG7", #NK_cell
                      "MBP", "MOG", "MOBP", "PLP1", "OLIG1", "OLIG2", # Oligodendrocyte
                      "GFAP", "AQP4", "ALDH1L1", "SLC1A2", # Astrocyte
                      "CLDN5", "PECAM1", "VWF", "CD34", "ENG", # Endothelial
                      "PDGFRB", "COL1A1", "ACTA2", "RGS5", "DES") # Fibroblast



#简单修饰的气泡热图,修改ggplot参数即可
DotPlot(pbmc, features =cell_type_markers) + coord_flip() +  # 使用DotPlot函数生成气泡图，使用coord_flip函数使坐标轴翻转
  theme_bw() +  # 应用黑白主题
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +  # 移除面板网格，调整x轴文本的水平和垂直对齐方式
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +  # 移除x和y轴标签，调整图例的顺序
  theme(axis.text.x = element_text(angle = 90)) +  # 将x轴文本旋转90度
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))  # 应用自定义的颜色渐变

####six.第一层注释####
#细胞大群注释
celltype=data.frame(ClusterID=0:19,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0,15),2]='Oligodendrocyte'
celltype[celltype$ClusterID %in% c(7),2]='Endothelial'
celltype[celltype$ClusterID %in% c(10,14),2]='T Cell'
celltype[celltype$ClusterID %in% c(1,2,18,19),2]='Macrophage'
celltype[celltype$ClusterID %in% c(3,4,5,6,8,9,11,12,13,16,17),2]='Glioma'
#celltype[celltype$ClusterID %in% c(15,18),2]='X'

celltype
table(celltype$celltype)
sce.in=pbmc
#先加一列celltype所有值为空，用于存放注释信息
sce.in@meta.data$celltype = "NA"
#> View(sce.in@meta.data) 可见多了列celltype为NA
###注释
for(i in 1:nrow(celltype)){
  sce.in@meta.data[which(sce.in@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.in@meta.data$celltype)



#整理各种类型细胞在各个cluster中的数量
table(sce.in@meta.data$celltype,sce.in@meta.data$seurat_clusters)


library(ggplot2)
sce=sce.in
#sce=sce_pca
p.dim.cell=DimPlot(sce, reduction = "tsne", group.by = "celltype",label = T,pt.size = 1) 
p.dim.cell
ggsave(plot=p.dim.cell,filename="tsne_celltype(新).pdf",width=8, height=8)

p.dim.cell=DimPlot(sce, reduction = "umap", group.by = "celltype",label = T,pt.size = 1,repel = T) 
p.dim.cell
ggsave(plot=p.dim.cell,filename="umap_celltype（新）.pdf",width=9, height=7)
saveRDS(sce,file = "GSE117891_注释后（新）.rds")



####查看基因的表达情况####
FeaturePlot(sce,features=c('ADAM17'),
            cols = c("lightgrey", 'red'),
            reduction = "tsne"
)

# 小提琴图（按亚群展示分布）
library(Seurat)
p1<-VlnPlot(sce, features =c('ADAM17'), group.by = "celltype")
p1
write.csv(p1, "ADAM17.csv", row.names = FALSE)


####一、提取成胶质瘤细胞亚聚类####
#别忘了默认分组为celltype
scRNA_harmony=readRDS("\\GSE117891_注释后.rds")
sce_pca=scRNA_harmony
Idents(sce_pca) = "celltype"
Idents(sce_pca)

#提取胶质瘤细胞
Glioma=subset(sce_pca,idents=c("Glioma"))
table(Glioma@meta.data$celltype)
#Glioma 
# 3034

#亚群聚类常规方法是提取counts来走作标准流程，这样做肯定没问题；
#其他方法也有很多，但具有争议这里不讲。

# seuratV5要这样提取counts
GetAssayData(Glioma, assay="RNA", layer='counts')
# V5这样提取count是没有行名的(另外一种方法，不建议)
Glioma@assays$RNA@layers$counts  #这样没有行名和列名，V5的锅

# 创建SeuratObject对象
Glioma_sce = CreateSeuratObject(counts = GetAssayData(Glioma, assay="RNA", layer='counts'),  # 使用提取的细胞构建新的Seurat对象
                             meta.data = Glioma@meta.data)  # 保留meta.data所有信息


#标注化、归一化、高变基因、pca
Glioma_sce <- NormalizeData(Glioma_sce) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)#对数据进行归一化、找高变基因、均一化、进行PCA降维
http://127.0.0.1:42699/graphics/plot_zoom_png?width=1200&height=900

#重点！！！！！！！！！
#这里要注意是否选择去批次，可选可不选都能解释通（建议都跑一遍）
# 比如我这里是上皮细胞里边有癌细胞，本身存在很强的样本异质性，去了批次反而去掉了样本间差异.
# 免疫细胞也可以不去批次，也存在较强异质性
# 我这里就不运行
# FB_sce <- RunHarmony(FB_sce, group.by.var = "orig.ident")

#降维聚类
Glioma_sce
ElbowPlot(Glioma_sce, ndims=30, reduction="pca")
Glioma_sce <- FindNeighbors(Glioma_sce, reduction = "pca", dims = 1:10)#reduction="harmony"
#分辨率树状图
library(clustree)
obj <- FindClusters(Glioma_sce, resolution = seq(0.1,1,by=0.1))
clustree(obj)

Glioma_sce = FindClusters(Glioma_sce,resolution = 0.1)  #0.477有13个亚组
table(Glioma_sce@meta.data$seurat_clusters)
#    0   1   2   3   4   5 
#   900 654 562 482 226 210 
Glioma_sce <- RunUMAP(Glioma_sce, reduction = "pca", dims = 1:10)##reduction="harmony"
Glioma_sce <- RunTSNE(Glioma_sce, reduction = "pca", dims = 1:10)

#Glioma_umap图_未注释
plot1 =DimPlot(Glioma_sce, reduction = "tsne",label = T,raster=FALSE) 
plot1
plot2 = DimPlot(Glioma_sce, reduction = "tsne", group.by='orig.ident',raster=FALSE) 
plot2
#split.by =""  通过什么信息来划分组别
plot3 = DimPlot(Glioma_sce, reduction = "tsne",split.by = "orig.ident",label = T,raster=FALSE) 
plot3
plot4 = DimPlot(Glioma_sce, reduction = "tsne",group.by = "orig.ident",shuffle = T,raster=FALSE) 
plot4
#自己手动保存图片吧
dev.off()


##（二）亚聚类注释###
table(Glioma_sce@meta.data$seurat_clusters)
#查看表达情况
cell_type_markers = c("SEC61G","TNFRSF12A",
                      "CCNB1","CENPF",
                      "FAM183A","C1orf194",
                      "CHCHD2","FKBP5",
                      "DLL3","GRIA2",
                      "SPP1",'HAMP',
                      'BCAS1','PLP1',
                      'TOP2A','HIST1H4C',
                      'PDGFRA','SOX11',
                      'MYBPHL','NNAT',
                      'CHI3L1','SOD2',
                      'AQP1','AQP4',
                      'SNRPE','CCT6A',"ADAM17"
)
duplicated_genes <- cell_type_markers[duplicated(cell_type_markers) | duplicated(cell_type_markers, fromLast = TRUE)]
#1.画气泡图
DotPlot(Glioma_sce, features = cell_type_markers)
#画热图
# 假设你已经有一个 Seurat 对象 Fib_sub
# 并且有一个数据框 top3，包含每个细胞亚群的前 3 个高表达基因

#2.绘制热图
DoHeatmap(Glioma_sce, features = cell_type_markers, size = 3) + NoLegend()

DoHeatmap(Glioma_sce,
          features = cell_type_markers,
          assay = 'RNA',
          group.colors = c("#00BFC4","#AB82FF","#00CD00","#C77CFF"))+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

DotPlot(Glioma_sce, features = cell_type_markers) + coord_flip() +  # 使用DotPlot函数生成气泡图，使用coord_flip函数使坐标轴翻转
  theme_bw() +  # 应用黑白主题
  theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5)) +  # 移除面板网格，调整x轴文本的水平和垂直对齐方式
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +  # 移除x和y轴标签，调整图例的顺序
  theme(axis.text.x = element_text(angle = 90)) +  # 将x轴文本旋转90度
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))  # 应用自定义的颜色渐变

ggsave(filename="17_0.477_markerHeatMap_Glioma.pdf",width=7, height=6)


#先直接用top gene注释吧
#M-c1-THBS1
celltype=data.frame(ClusterID=0:5,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0),2]='G1'
celltype[celltype$ClusterID %in% c(1),2]='G2'
celltype[celltype$ClusterID %in% c(2),2]='G3'
celltype[celltype$ClusterID %in% c(3),2]='G4'
celltype[celltype$ClusterID %in% c(4),2]='G5'
celltype[celltype$ClusterID %in% c(5),2]='G6'
celltype[celltype$ClusterID %in% c(6),2]='G7'
celltype[celltype$ClusterID %in% c(7),2]='G8'
celltype[celltype$ClusterID %in% c(8),2]='G9'
celltype[celltype$ClusterID %in% c(9),2]='G10'
celltype[celltype$ClusterID %in% c(10),2]='G11'
celltype[celltype$ClusterID %in% c(11),2]='G12'
celltype[celltype$ClusterID %in% c(12),2]='G13'

celltype
table(celltype$celltype)
sce.in=Glioma_sce
#先加一列celltype所有值为空，用于存放注释信息
sce.in@meta.data$celltype = "NA"
View(sce.in@meta.data) #可见多了列celltype为NA
###注释
for(i in 1:nrow(celltype)){"http://127.0.0.1:45663/graphics/plot_zoom_png?width=1200&height=900&index="
  sce.in@meta.data[which(sce.in@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.in@meta.data$celltype)
View(sce.in@meta.data)

library(ggplot2)
sce=sce.in

p.dim.cell=DimPlot(sce, reduction = "umap", group.by = "celltype",label = T,pt.size = 1,repel = T) 
p.dim.cell
ggsave(plot=p.dim.cell,filename="Gliama_umap6clu_注释.pdf",width=6, height=4)


sce <- RunTSNE(sce, reduction = "pca", dims = 1:25)
p.dim.cell=DimPlot(sce, reduction = "tsne", group.by = "celltype",label = T,pt.size = 1) 
p.dim.cell
ggsave(plot=p.dim.cell,filename="Gliama_tsne_6clu_注释.pdf.pdf",width=8, height=8)

saveRDS(sce,file = "ADAM17_Gliama_tsne+umap_6clu_注释.rds")#得要saveRDS!!


##featurePlot函数查看单个基因的feature图
FeaturePlot(sce,features = "ADAM17",pt.size = 1,
            cols = c("lightgrey", 'red'),
            min.cutoff="q10",max.cutoff="q90")

FeaturePlot(sce,features = "ADAM17",pt.size = 1,
            cols = c("lightgrey", 'red'))

FeaturePlot(sce,features = "ADAM17",pt.size = 1)

# 1. 提取表达数据和分组信息（例如按 "seurat_clusters" 分组）
plot_df <- FetchData(sce, vars = c("ADAM17", "seurat_clusters"))

# 2. 使用 ggplot2 绘制箱线图
library(ggplot2)
ggplot(plot_df, aes(x = seurat_clusters, y = ADAM17, fill = seurat_clusters)) +
  geom_boxplot() +
  labs(title = "ADAM17 expression by cluster", x = "Cluster", y = "Log-normalized expression") +
  theme_minimal() +
  theme(legend.position = "none")


##使用Seurat包查看某个基因在细胞亚群中的表达量
# 检查基因是否在数据中
gene_name <- "ADAM17"  # 例如 "CD4", "SOX2"
if (!gene_name %in% rownames(sce)) {
  stop("基因不存在于数据中！")
}
# 提取表达矩阵（使用标准化后的数据：data；原始计数：counts）
expression_matrix <- GetAssayData(sce, slot = "counts")  # 或 slot = "counts"

# 按亚群分组统计（均值、中位数、检出率等）
cell_groups <- sce$celltype  # 替换为你的亚群列名
table(cell_groups)
gene_expression <- data.frame(
  CellType = cell_groups,
  Expression = expression_matrix[gene_name, ]
)

# 汇总统计
library(dplyr)
summary_stats <- gene_expression %>%
  group_by(CellType) %>%
  summarise(
    Mean = mean(Expression), #均值
    Median = median(Expression),  #中位数
    DetectionRate = sum(Expression > 0) / n()  #检出率
  )
print(summary_stats)
# 小提琴图（按亚群展示分布）
library(Seurat)
p1<-VlnPlot(sce, features =gene_name, group.by = "celltype")
p1

ggsave(plot=p1,filename="ADAM17_VlnPlot_Ast_3clu.pdf.pdf",width=6, height=4)
# 点图（展示平均表达和检出率）
DotPlot(sce, features = gene_name, group.by = "celltype") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



####seven.火山图+GO+KEGG+GSEV####
rm(list=ls())
library(msigdbr) #提供MSigdb数据库基因集
library(gplots)
library(ggplot2)
library(clusterProfiler)#Y叔富集分析
library(org.Hs.eg.db)
library(Seurat) # 单细胞数据处理核心R包
scRNA_harmony = scRNA_6148_注释后
# 先看一下注释好的umap图
DimPlot(scRNA_harmony, reduction = "umap", group.by = "celltype",label = T,pt.size = 1,repel = T)
#查看各群细胞数量
table(scRNA_harmony$seurat_clusters)
# 0   1   2   3   4   5   6   7   8   9  10  11  12  13 
#889 861 607 494 458 429 352 308 304 268 230 210 209 186 


##(一).差异分析##
# 1、所有cluster整体的差异分析
Idents(scRNA_harmony)="celltype"#先把Idents改回来
# 找出每个cluster的标记与所有剩余的细胞相比较，only.pos = T只报告阳性细胞
cell.markers <- FindAllMarkers(object = scRNA_harmony, 
                               only.pos = FALSE, # 是否只保留表达相对上调的基因，设置FALSE则会保留下调的
                               test.use = "wilcox", # 默认使用 wilcox 非参数检验，其它选项可以查看说明
                               slot = "data", # 需要注意的是，默认使用 data，而不是 counts
                               min.pct = 0.5, # 设置表达比例的阈值，没有统一标准，差异基因很多的情况下可以把阈值调高，差异基因有1000个够用
                               logfc.threshold = 0.5 # 设置 log2FC 即差异倍率的阈值，没有统一标准，差异基因很多的情况下可以把阈值调高
)
#看戏细胞和cluter分布
table(scRNA_harmony@meta.data$celltype,scRNA_harmony@meta.data$seurat_clusters)
# 2、两个cluster差异分析
#pct.1表示在聚类簇1中表达该基因的单细胞样本占总单细胞样本数的比例，而pct.2则表示在聚类簇2中表达该基因的单细胞样本占总单细胞样本数的比例

G_degs = FindMarkers( scRNA_harmony, 
                                logfc.threshold = 0.25,
                                min.pct = 0.1, # 表达比例的阈值设置，小一点找出更多差异基因
                                only.pos = FALSE, #false既保留上调的基因，也保留下调的基因
                                ident.1 = "G4", ident.2 = "G3") %>% # ident 1 vs 2
  mutate( gene = rownames(.) ) # 将行名新增为列，方便检索基因



# 自定义筛选
G_degs_fil = G_degs %>% 
  filter( pct.1 > 0.2 & p_val_adj < 0.05 ) %>% # &：和；|：或。 
  filter( abs( avg_log2FC ) > 0.5 ) #%>% # 按 logFC 绝对值筛选
# filter( avg_log2FC > 0 ) # 按 logFC 正负筛选


# 火山图
library(ggrepel) 
# 查看 T 细胞差异基因数据的列名
colnames(G_degs_fil)
# 统计绝对平均对数折叠变化大于2的基因数量，变化的数值可以自己调
table(abs(G_degs_fil$avg_log2FC) > 2) 
# 创建用于绘制火山图的数据集
plotdt = G_degs_fil %>% 
  mutate(gene = ifelse(abs(avg_log2FC) >= 2, gene, NA))
# 选择平均对数折叠变化大于等于2的基因，创建新的数据集
# 绘制火山图
ggplot(plotdt, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                   size = pct.1, # 按照表达比例 pct.1 设置散点大小
                   color = avg_log2FC)) +
  geom_point() + # 绘制散点图
  ggtitle(label = "G4 vs G2", subtitle = "Antibody: elevated vs control") + # 添加标题和副标题
  geom_text_repel(aes(label = gene), size = 3, color = "black") + # 加上 top 基因的标签
  theme_bw() + # 设置主题为白色背景
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), # 设置标题字体样式
    plot.background = element_rect(fill = "transparent",colour = NA) # 设置图表背景透明
  ) +
  scale_color_gradient2(low = "olivedrab", high = "salmon2", 
                        mid = "grey", midpoint = 0) + # 设置颜色渐变
  scale_size(range = c(1,3)) # 设置点的大小范围




####GO/KEGG富集####
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
#为每个基因添加对应的ENTREZID
G_degs_fil$gene <- rownames(G_degs_fil)
ids=bitr(G_degs_fil$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
#合并数据，cluser3.markers中没有ENTREZID的基因将被过虑掉
G_degs_fil=merge(G_degs_fil,ids,by.x='gene',by.y='SYMBOL')
#查看数据结构
head(G_degs_fil)
# > head(Tcells_degs)
# gene         p_val avg_log2FC pct.1 pct.2     p_val_adj ENTREZID
# 1  AAK1  3.622514e-11  1.2998420 0.572 0.378  7.540987e-07    22848
# 2 AAMDC  4.991934e-43 -0.3475556 0.043 0.433  1.039171e-38    28971
# 3  AASS 2.328826e-101 -0.7607128 0.021 0.567  4.847917e-97    10157
# 4 ABCA3 1.055439e-185 -1.7731186 0.011 0.756 2.197108e-181       21
# 5 ABCC6  8.980491e-86 -0.2696661 0.001 0.300  1.869469e-81      368
# 6 ABCD3  1.160593e-80 -0.6544183 0.039 0.600  2.416007e-76     5825
#将基因按照avg_log2FC的大小进行降序排列
G_degs_fil <- G_degs_fil[order(G_degs_fil$avg_log2FC,decreasing = T),]
#生成仅含有ENTREZID名字和avg_log2FC值的gene list
G_degs_list <- as.numeric(G_degs_fil$avg_log2FC)
names(G_degs_list) <- G_degs_fil$ENTREZID
head(G_degs_list)
#  5320    27306    10261    57817    10631   120939 
# 5.972006 5.776493 5.761866 5.472381 5.422898 5.342027  

# enrichGO
#筛选差异较大的基因集
cluster3_de <- names(G_degs_list)[abs(G_degs_list) > 2]
head(cluster3_de)
# [1]"5320"   "27306"  "10261"  "57817"  "10631"  "120939"
#go富集
cluster3_ego <- enrichGO(cluster3_de, OrgDb = "org.Hs.eg.db", ont="ALL", readable=TRUE)
head(cluster3_ego)
#气泡图
dotplot(cluster3_ego, showCategory=10,title="G4 vs G3 GO")
barplot(cluster3_ego, split ="ONTOLOGY", font.size =10,  showCategory=10) +facet_grid(ONTOLOGY ~ ., space ="free_y",scales ="free_y")

#KEGG富集
cluster3_ekg <- enrichKEGG(gene= cluster3_de, organism = "hsa",pvalueCutoff = 0.05)
head(cluster3_ekg)
#气泡图
barplot(cluster3_ekg)
dotplot(cluster3_ekg, showCategory=10,title="G4 vs G3 KEGG")


##GSEA分析##
G_degs = FindMarkers( scRNA_harmony, 
                                logfc.threshold = 0.25,
                                min.pct = 0.1, # 表达比例的阈值设置，小一点找出更多差异基因
                                only.pos = FALSE, #false既保留上调的基因，也保留下调的基因
                                ident.1 = "G4", ident.2 = "G2") %>% # ident 1 vs 2
  mutate( gene = rownames(.) ) # 将行名新增为列，方便检索基因

# GSEA需要使用整个gene list，这里进行kegg分析
#为每个基因添加对应的ENTREZID
G_degs$gene <- rownames(G_degs)
ids=bitr(G_degs$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
#合并数据，cluser3.markers中没有ENTREZID的基因将被过虑掉
G_degs=merge(G_degs,ids,by.x='gene',by.y='SYMBOL')
#查看数据结构
head(G_degs)

#将基因按照avg_log2FC的大小进行降序排列
G_degs <- G_degs[order(G_degs$avg_log2FC,decreasing = T),]
#生成仅含有ENTREZID名字和avg_log2FC值的gene list
cluster3.markers_list <- as.numeric(G_degs$avg_log2FC)
names(cluster3.markers_list) <- G_degs$ENTREZID
head(cluster3.markers_list)



cluster3_gsekg <- gseKEGG(cluster3.markers_list,organism = "hsa",pvalueCutoff = 0.05)
head(cluster3_gsekg)
#将富集结果按照NES绝对值降序排列
cluster3_gsekg_arrange <- arrange(cluster3_gsekg,desc(abs(NES)))
head(cluster3_gsekg_arrange)
library(clusterProfiler)


#作图 install.packages("ggupset")
library(ggupset)

color <- c("#f7ca64", "#43a5bf", "#86c697", "#a670d6", "#ef998a")
gsekp1 <- gseaplot2(cluster3_gsekg_arrange, 1:5, color = color, pvalue_table=F, base_size=14)
gsekp1
gsekp2 <- upsetplot(cluster3_gsekg_arrange, n=5)
gsekp2
cowplot::plot_grid(gsekp1, gsekp2, rel_widths=c(1, .6), labels=c("A", "B"))
gsekp1+gsekp2

# 使用 cowplot::plot_grid 组合
cowplot::plot_grid(gsekp1_subplot, gsekp2_subplot, rel_widths = c(1, .6), labels = c("A", "B"))

