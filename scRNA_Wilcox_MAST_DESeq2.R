
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
setwd("D:\\R language\\5.ADAM17\\NO2.scRNA")
getwd()

# 读取 TSV 文件
count_matrix <- read.delim("D:\\R language\\5.ADAM17\\NO2.scRNA\\CGGA.scRNAseq_6148.count.matrix.tsv", row.names = 1)
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
pbmc <- JackStraw(pbmc, num.replicate = 50)
pbmc <- ScoreJackStraw(pbmc, dims = 1:15)
JackStrawPlot(pbmc, dims = 1:15)
#Elbow plot
ElbowPlot(pbmc)

##4.细胞聚类(Seurat使用KNN算法进行聚类。)
#dims = 1:10 即选取前10个主成分来分类细胞。
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
#查看前5个细胞的分类ID
head(Idents(pbmc), 5)

##5.非线性降维##
pbmc <- RunUMAP(pbmc, dims = 1:10)
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
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc<-JoinLayers(pbmc)
markers <- FindAllMarkers(object = pbmc, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25,
                          min.pct = 0.5)  #要限制一下在pct里面的一个最低表达
#要注意加这个min.pct=0.5的限制，不然可能会出现假阳性高变基因的情况

saveRDS(markers,file = "FindMarker.rds")
#直接双击load吧

#寻找每个簇的高变基因
library(dplyr)
all.markers = markers %>% 
  dplyr::select(gene, everything()) %>% 
  subset(p_val_adj<0.05)
#将每个cluster lgFC排在前10的marker基因挑选出来
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
View(top10)
#将top gene放在一起
gene_column <- top10[ , "gene", drop = FALSE]
gene_vector <- top10[ , "gene"]
view(gene_vector)


####20.singleR自动注释####
# SingleR安装包下载地址：https://www.bioconductor.org/packages/release/bioc/html/SingleR.html
# 安装并加载 SingleR 和 celldex 包
BiocManager::install("SingleR")        
BiocManager::install("celldex")        

# 加载 celldex 和 SingleR 库
library(celldex)
library(SingleR)

# 列出 celldex 包中的所有数据集
ls("package:celldex")
# [1] "BlueprintEncodeData" 人类的免疫细胞             "DatabaseImmuneCellExpressionData"人类免疫细胞
# [3] "HumanPrimaryCellAtlasData"人广泛的细胞类型注释        "ImmGenData" 小鼠免疫细胞                     
# [5] "MonacoImmuneData" 人类免疫细胞                "MouseRNAseqData"  小鼠单细胞数据的广泛细胞类型注释               
# [7] "NovershternHematopoieticData" 人类血液 

# 加载 HumanPrimaryCellAtlasData 作为参考数据集
HumanPrimaryCellAtlas <- HumanPrimaryCellAtlasData()
HumanPrimaryCellAtlas <- celldex::HumanPrimaryCellAtlasData()


# 加载之前降维聚类后未注释的数据
load("scRNA_harmony2.Rdata")
p1=DimPlot(pbmc,label = T)#umap图
p1
scRNA <- pbmc   # 将 scRNA_harmony 数据赋值给 scRNA

# 提取data数据
test <- GetAssayData(scRNA, layer="data")  # 从 scRNA 对象中提取data数据

# 显示当前 scRNA 对象中各 cluster 的细胞数
table(scRNA@active.ident)
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13 
#889 861 607 494 458 429 352 308 304 268 230 210 209 186 


# 使用 SingleR 进行细胞类型注释
Rename_scRNA <- SingleR(
  test = test,                                # 目标数据
  ref = HumanPrimaryCellAtlas,                # 参考数据集
  labels = HumanPrimaryCellAtlas$label.main,  # 使用粗粒度的标注信息（label.fine 可用于更细粒度的标注，但是很容易出错）
  clusters = scRNA@active.ident               # 根据 scRNA 中的 cluster 注释
)

# 查看注释后的细胞类型标签
Rename_scRNA$labels          # 输出每个 cluster 的标注结果
Rename_scRNA$delta.next      # 输出 delta.next 值，表示每个 cluster 的第二高匹配度
# 较大的 delta.next 值表示分类的置信度较高，因为第一和第二标签之间的分数差异大；较小的 delta.next 值则表示置信度较低，分类结果可能不太确定。
Rename_scRNA$pruned.labels   # 输出去除低置信度的标注标签
# 当 delta.next 的值过小时，意味着该细胞的分类不够可靠，SingleR 会将这些分类结果设置为 NA，即未分类


# 更新 scRNA 对象中的 cluster 名称
new.cluster.ids <- Rename_scRNA$pruned.labels  # 获取最终的注释结果
names(new.cluster.ids) <- levels(scRNA)        # 将 cluster 名称赋予 new.cluster.ids
scRNA <- RenameIdents(scRNA, new.cluster.ids)  # 重命名 scRNA 对象中的 cluster
Idents(scRNA)
#新增一列celltype保存注释结果
scRNA$celltype_singleR=Idents(scRNA)

# 绘制 UMAP 图，显示注释结果
p2 <- DimPlot(scRNA, reduction="umap", label=T, pt.size=0.5,group.by = "celltype_singleR") + NoLegend()
p2  # 显示 UMAP 图
p1+p2
# 查看注释结果的热图
plotScoreHeatmap(Rename_scRNA)   # 绘制细胞类型匹配得分的热图，以便可视化注释的置信度



####six.第一层注释####
#细胞大群注释
celltype=data.frame(ClusterID=0:13,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(10,13),2]='T cell'
celltype[celltype$ClusterID %in% c(1),2]='Neuron cell'
celltype[celltype$ClusterID %in% c(0),2]='Macrophage'
celltype[celltype$ClusterID %in% c(5,7),2]='Monocyte'
celltype[celltype$ClusterID %in% c(2,3,4,6,8,9,11,12),2]='Astrocyte'
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
View(sce.in@meta.data)


#整理各种类型细胞在各个cluster中的数量
table(sce.in@meta.data$celltype,sce.in@meta.data$seurat_clusters)


library(ggplot2)
sce=sce.in
p.dim.cell=DimPlot(sce, reduction = "tsne", group.by = "celltype",label = T,pt.size = 1) 
p.dim.cell
ggsave(plot=p.dim.cell,filename="tsne_celltype.pdf",width=9, height=7)

p.dim.cell=DimPlot(sce, reduction = "umap", group.by = "celltype",label = T,pt.size = 1,repel = T) 
p.dim.cell
ggsave(plot=p.dim.cell,filename="umap_celltype.pdf",width=9, height=7)
saveRDS(sce,file = "scRNA 6148_注释后.rds")



####查看基因的表达情况####
FeaturePlot(sce,features=c('CHI3L1','OLIG2'),
            cols = c("lightgrey", 'red'),
           )

# 小提琴图（按亚群展示分布）
library(Seurat)
p1<-VlnPlot(sce, features =c('CHI3L1','OLIG2'), group.by = "celltype")
p1




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

Macrophages_degs = FindMarkers( scRNA_harmony, 
                                logfc.threshold = 0.25,
                                min.pct = 0.1, # 表达比例的阈值设置，小一点找出更多差异基因
                                only.pos = FALSE, #false既保留上调的基因，也保留下调的基因
                                ident.1 = "Macrophage", ident.2 = "Astrocyte") %>% # ident 1 vs 2
  mutate( gene = rownames(.) ) # 将行名新增为列，方便检索基因



# 自定义筛选
Macrophages_degs_fil = Macrophages_degs %>% 
  filter( pct.1 > 0.2 & p_val_adj < 0.05 ) %>% # &：和；|：或。 
  filter( abs( avg_log2FC ) > 0.5 ) #%>% # 按 logFC 绝对值筛选
# filter( avg_log2FC > 0 ) # 按 logFC 正负筛选


# 火山图
library(ggrepel) 
# 查看 T 细胞差异基因数据的列名
colnames(Macrophages_degs_fil)
# 统计绝对平均对数折叠变化大于2的基因数量，变化的数值可以自己调
table(abs(Macrophages_degs_fil$avg_log2FC) > 2) 
# 创建用于绘制火山图的数据集
plotdt = Macrophages_degs_fil %>% 
  mutate(gene = ifelse(abs(avg_log2FC) >= 2, gene, NA))
# 选择平均对数折叠变化大于等于2的基因，创建新的数据集
# 绘制火山图
ggplot(plotdt, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                   size = pct.1, # 按照表达比例 pct.1 设置散点大小
                   color = avg_log2FC)) +
  geom_point() + # 绘制散点图
  ggtitle(label = "Macrophage vs Astrocyte", subtitle = "Antibody: elevated vs control") + # 添加标题和副标题
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
Macrophages_degs_fil$gene <- rownames(Macrophages_degs_fil)
ids=bitr(Macrophages_degs_fil$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
#合并数据，cluser3.markers中没有ENTREZID的基因将被过虑掉
Macrophages_degs_fil=merge(Macrophages_degs_fil,ids,by.x='gene',by.y='SYMBOL')
#查看数据结构
head(Macrophages_degs_fil)
# > head(Tcells_degs)
# gene         p_val avg_log2FC pct.1 pct.2     p_val_adj ENTREZID
# 1  AAK1  3.622514e-11  1.2998420 0.572 0.378  7.540987e-07    22848
# 2 AAMDC  4.991934e-43 -0.3475556 0.043 0.433  1.039171e-38    28971
# 3  AASS 2.328826e-101 -0.7607128 0.021 0.567  4.847917e-97    10157
# 4 ABCA3 1.055439e-185 -1.7731186 0.011 0.756 2.197108e-181       21
# 5 ABCC6  8.980491e-86 -0.2696661 0.001 0.300  1.869469e-81      368
# 6 ABCD3  1.160593e-80 -0.6544183 0.039 0.600  2.416007e-76     5825
#将基因按照avg_log2FC的大小进行降序排列
Macrophages_degs_fil <- Macrophages_degs_fil[order(Macrophages_degs_fil$avg_log2FC,decreasing = T),]
#生成仅含有ENTREZID名字和avg_log2FC值的gene list
Macrophages_degs_list <- as.numeric(Macrophages_degs_fil$avg_log2FC)
names(Macrophages_degs_list) <- Macrophages_degs_fil$ENTREZID
head(Macrophages_degs_list)
# 7852    28755     3575     5552     9235     4050 
# 3.574470 3.070149 3.062951 2.980478 2.917090 2.816810 

# enrichGO
#筛选差异较大的基因集
cluster3_de <- names(Macrophages_degs_list)[abs(Macrophages_degs_list) > 1]
head(cluster3_de)
# [1] "7852"  "28755" "3575"  "5552"  "9235"  "4050" 
#go富集
cluster3_ego <- enrichGO(cluster3_de, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
head(cluster3_ego)
#气泡图
dotplot(cluster3_ego, showCategory=10,title="Macrophage vs Astrocyte GO")


#KEGG富集
cluster3_ekg <- enrichKEGG(gene= cluster3_de, organism = "hsa",pvalueCutoff = 0.05)
head(cluster3_ekg)
#气泡图
dotplot(cluster3_ekg, showCategory=10,title="Macrophages vs Fibroblasts KEGG")


##GSEA分析##
Macrophages_degs = FindMarkers( scRNA_harmony, 
                                logfc.threshold = 0.25,
                                min.pct = 0.1, # 表达比例的阈值设置，小一点找出更多差异基因
                                only.pos = FALSE, #false既保留上调的基因，也保留下调的基因
                                ident.1 = "Macrophages", ident.2 = "Fibroblasts") %>% # ident 1 vs 2
  mutate( gene = rownames(.) ) # 将行名新增为列，方便检索基因

# GSEA需要使用整个gene list，这里进行kegg分析
#为每个基因添加对应的ENTREZID
Macrophages_degs$gene <- rownames(Macrophages_degs)
ids=bitr(Macrophages_degs$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
#合并数据，cluser3.markers中没有ENTREZID的基因将被过虑掉
Macrophages_degs=merge(Macrophages_degs,ids,by.x='gene',by.y='SYMBOL')
#查看数据结构
head(Macrophages_degs)

#将基因按照avg_log2FC的大小进行降序排列
Macrophages_degs <- Macrophages_degs[order(Macrophages_degs$avg_log2FC,decreasing = T),]
#生成仅含有ENTREZID名字和avg_log2FC值的gene list
cluster3.markers_list <- as.numeric(Macrophages_degs$avg_log2FC)
names(cluster3.markers_list) <- Macrophages_degs$ENTREZID
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
gsekp2 <- upsetplot(cluster3_gsekg_arrange, n=5)
cowplot::plot_grid(gsekp1, gsekp2, rel_widths=c(1, .6), labels=c("A", "B"))
gsekp1+gsekp2

# 使用 cowplot::plot_grid 组合
cowplot::plot_grid(gsekp1_subplot, gsekp2_subplot, rel_widths = c(1, .6), labels = c("A", "B"))


