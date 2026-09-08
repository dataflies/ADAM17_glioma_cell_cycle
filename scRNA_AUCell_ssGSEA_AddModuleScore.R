rm(list=ls())
setwd("D:\\R language\\5.ADAM17\\NO2.scRNA\\单细胞评分")


#####单细胞各种评分GSVA/AUCcell/AddModuleScore/ssGSEA等等####
####17、单细胞评分一：GSVA/ssGSEA####
library(Seurat)
library(gplots)
library(ggplot2) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA) 
library(GSEABase)
sce <- readRDS("D:/R language/5.ADAM17/NO2.scRNA/(1)-analyze data/scRNA_6148_注释后.rds")
setwd("单细胞打分/")
DimPlot(sce, reduction = "umap", group.by = "celltype",label = T,pt.size = 1,repel = T) 



#1、准备单细胞平均表达值矩阵
colnames(sce@meta.data)#查看meta.data信息
#获取celltype_epi的平均表达值
expr <-AverageExpression(sce ,
                         group.by = "celltype",#这里还可以是group/分组等等
                         assays = "RNA",
                         slot = "data")
expr=expr[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)#转换成matrix，GSVA分析需要用
head(expr)
write.csv(expr,file = 'AverageExpression.csv')
# 单细胞的GSVA还有一种方法就是直接取单细胞的data数据作GSVA，不过样本量大不推荐使用，其他打分可以代替
# expr <- as.data.frame(GetAssayData(sce, assay="RNA", layer ='data'))



cg = names(tail(sort(apply(expr, 1, sd)), 1000))  # 找到基因表达矩阵 `expr` 中标准差最高的 1000 个基因的名称
# `apply(expr, 1, sd)` 计算 `expr` 矩阵中每个基因（行）的标准差
pheatmap::pheatmap(cor(expr[cg, ]))  # 绘制热图显示这些基因表达的相关性

# `cor(...)` 计算这些基因表达数据之间的相关性矩阵
pheatmap::pheatmap(cor(expr[cg,]),
                   file = 'AverageExpression.pdf')

dev.off()



# 2、准备通路基因集；利用GSVA进行分析
# 可以直接下载gmt文件：
# https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H
library(msigdbr)
# install.packages("msigdbr")
all_gene_sets = msigdbr(species = "Homo sapiens",#Mus musculus
                        category='H') #从Msigdb中获取信号通路(C1-C7;H)，'H' 类别 (Hallmark gene sets)

# 将 all_gene_sets 数据框中的基因符号按基因集名称分组
gs = split(all_gene_sets$gene_symbol, all_gene_sets$gs_name)  
# 对每个基因集中的基因符号去重
gs = lapply(gs, unique)  # 对基因符号去重
head(gs)  
# 将基因集转换为 GeneSetCollection 对象
genesets  <- GeneSetCollection(mapply(function(geneIds, keggId) {  # 将基因集转换为 GeneSetCollection 对象
  GeneSet(geneIds, geneIdType = EntrezIdentifier(),         # 创建 GeneSet 对象，使用 EntrezIdentifier 作为基因 ID 类型
          collectionType = KEGGCollection(keggId),          # 使用 KEGGCollection 作为基因集类型
          setName = keggId)                                  # 设置基因集名称为 keggId
}, gs, names(gs)))                                          # gs 为基因 ID 列表，names(gs) 为基因集名称
genesets  # geneset


# 3、GSVA分析
gsva_result <- gsva(as.matrix(expr), genesets,  # expr 为表达矩阵，geneset 为基因集集合
                    mx.diff = FALSE,  # 使用标准化富集分数（NES），而不是最大差异
                    parallel.sz = 5)  # 设置并行计算的线程数为 8
#method参数可以选取其他方法，默认是gsva：
# method = c("gsva", "ssgsea", "zscore", "plage")

# 标准化后的数据用Gaussian(高斯分布)，原始的count用“Poisson”

# 目前下载GSVA都是更新之后版本，最新是V1.52，运行方式不同会报错，
# 推荐还是使用老版本：（下载包到本地安装）：
# https://mghp.osn.xsede.org/bir190004-bucket01/archive.bioconductor.org%2Fpackages%2F3.17%2Fbioc%2Fsrc%2Fcontrib%2FArchive%2FGSVA%2FGSVA_1.48.0.tar.gz



# 保存gsva分析所需要的结果
write.csv(gsva_result,file='gsva_result.csv')
head(gsva_result)
#gsva_result=read.csv("gsva_result.csv")

# 热图可视化1
library(pheatmap)
pheatmap(gsva_result)
pheatmap(gsva_result,show_rownames = F)#可视化
ggsave('pheatmap_by_GSVA_all.pdf',width = 7,height = 5)
# 可视化2
pheatmap(gsva_result, show_colnames = T, 
         scale = "row",angle_col = "45",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
ggsave('pheatmap_by_GSVA_all（2）.pdf',width = 7,height = 5)
dev.off()
# 气泡图
library( reshape2)
gsva_long <- melt(gsva_result, id.vars = "Genesets")
ggplot(gsva_long, aes(x = Var2, y = Var1, size = value, color = value)) +
  geom_point(alpha = 0.7) +  # 使用散点图层绘制气泡，alpha设置点的透明度
  scale_size_continuous(range = c(1, 6)) +  # 设置气泡大小的范围
  theme_minimal() + 
  scale_color_gradient(low='#008020',high='#08519C') +
  labs(x = "Gene Set", y = "Sample", size = "GSVA Score")+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
ggsave('pheatmap_by_GSVA_all（3）.pdf',width = 7,height = 5)


# 4.挑选最显著的top5条通路作图
cg = names(tail(sort(apply(gsva_result, 1, sd) ),5)) # 计算每个通路（行）的标准差，排序后选择标准差最大的10个通路名称
pheatmap(gsva_result[cg,]) # 绘制选中通路的热图
ggsave('pheatmap_by_GSVA_sd5.pdf', width = 7, height = 5) # 保存热图为PDF文件，设置宽度为7，高度为5



####后续可以用GSVA的结果做预后分析（score高低分组做生存分析），
# 或者差异分析(通路表达矩阵做差异分析)


####可以看到GSVA由于算法的局限性常规还是用每个celltype的平均表达值做评分，
# 并没有每一个细胞作评分，接着我们用AddModuleScore作每个细胞的评分


####如何自定义基因集gmt文件进行打分?####
#比较简单，我们用excel在电脑上随便打开一个gmt文件
#然后我们在修改里边的内容，按照gmt的内容修改，第二列设置为NA就行，保存
#最后读取我们自定义的基因集：
GeneSet_1 <- getGmt("./gmt文件/h.all.v6.2.symbols.gmt")

#最后参照上边代码正常GSVA分析即可





####18、单细胞评分二：AddModuleScore####
load("../Epi_已注释.Rdata")
# setwd("单细胞打分/")
DimPlot(sce)#umap图


#加载h.all.v6.2.symbols的gmt基因集
library(clusterProfiler)  
Hallmarker <- read.gmt("gmt文件/h.all.v6.2.symbols.gmt")  # 读取 GMT 文件，包含基因集和基因信息
GeneSets <- split(Hallmarker$gene, Hallmarker$term)  # 将 Hallmarker 数据框中的基因按照基因集名称进行分组
names(GeneSets) <- unique(Hallmarker$term)  # 将 GeneSets 列表的名称设置为每个基因集的名称

# AddModuleScore评分
# 计算每个基因集中的基因在每个细胞中的平均表达值
sce <- AddModuleScore(object = sce, features = GeneSets, name = names(GeneSets))

# 可视化
FeaturePlot(sce,features = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB1')
FeaturePlot(sce,features = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB1',split.by = 'group' )
DotPlot(sce, features = c('HALLMARK_TNFA_SIGNALING_VIA_NFKB1'))

names(sce@meta.data)#
#后两个通路为 顶端连接和 顶端表面，根据自己研究选择通路(也可以自定义，根据其他文献的基因集来打分等)
DotPlot(sce, features = c('HALLMARK_TNFA_SIGNALING_VIA_NFKB1','HALLMARK_P53_PATHWAY37',
                          'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION30','HALLMARK_APICAL_SURFACE21'))


#通过gsea数据库选择C2的GOBP_CELL_CYCLE打分试试
cell_cycling <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::filter(gs_name=="GOBP_CELL_CYCLE")

cycling_l <- list(cell_cycling=unique(cell_cycling$gene_symbol)) 
sce_T <- AddModuleScore(sce,features = cycling_l,name = "GOBP_CELL_CYCLE")
names(sce_T@meta.data)
VlnPlot(sce_T,features = 'GOBP_CELL_CYCLE1', pt.size = 0, adjust = 2)
DotPlot(sce_T, features = c('GOBP_CELL_CYCLE1'))



####自定义基因集打分####
proliferation = list(c('MKI67','IGF1','ITGB2','PDGFC','JAG1','PHGDH'))#增殖
migration = list(c('VIM','SNAI1','MMP9','AREG','ARID5B' ,'FAT1'))#转移
#增殖打分
proliferation_sce <- AddModuleScore(sce,
                                    features = proliferation,
                                    ctrl = 100,
                                    name = "proliferation1")
names(proliferation_sce@meta.data)
VlnPlot(proliferation_sce,features = 'proliferation11', pt.size = 0, adjust = 2)
DotPlot(proliferation_sce, features = c('proliferation11'))

#转移打分
migration_sce <- AddModuleScore(sce,
                                features = migration,
                                ctrl = 100,
                                name = "migration")
names(migration_sce@meta.data)
VlnPlot(migration_sce,features = 'migration1', pt.size = 0, adjust = 2)
DotPlot(migration_sce, features = c('migration1'))


#画UMAP图，将score映射到细胞
library(ggplot2)
mydata<- FetchData(migration_sce,vars = c("umap_1","umap_2","migration1","group"))

a <- ggplot(mydata,aes(x = umap_1,y =umap_2,colour = migration1))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) 
# facet_wrap(~group)

VlnPlot(sce,features = 'AUC', pt.size = 0, adjust = 2)





####19、单细胞评分三：AUcell####
# https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html
# BiocManager::install("AUCell")
library(AUCell)
library(clusterProfiler)
library(ggplot2)
setwd("单细胞打分/")
load("../Epi_已注释.Rdata")
DimPlot(sce)#umap图


#选择基因集(这里只列举了3种选取GeneSets的方法)
# 1、(下载gmt文件或自定义)
Hallmarker <- read.gmt("gmt文件/h.all.v6.2.symbols.gmt")  # 读取 GMT 文件，包含基因集和基因信息
GeneSets <- split(Hallmarker$gene, Hallmarker$term)  # 将 Hallmarker 数据框中的基因按照基因集名称进行分组
#选取任意基因集进行打分
geneSets = GeneSets[c(2,4,6,8,10,12,14,16,18)]
geneSets = GeneSets[c(36:50)]
geneSets = GeneSets[c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                      "HALLMARK_HYPOXIA",
                      "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
                      "HALLMARK_TGF_BETA_SIGNALING",
                      "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                      "HALLMARK_DNA_REPAIR",
                      "HALLMARK_G2M_CHECKPOINT",
                      "HALLMARK_APOPTOSIS",
                      "HALLMARK_NOTCH_SIGNALING",
                      "HALLMARK_P53_PATHWAY",
                      "HALLMARK_MYC_TARGETS_V1",
                      "HALLMARK_MYC_TARGETS_V2")]

#2、选取C7通路(包含了GO term)
library("msigdbr")
C7 = msigdbr(species = "Homo sapiens",#Mus musculus
             category='C7') #从Msigdb中获取信号通路(C1-C7;H)，'H' 类别 (Hallmark gene sets)
# 将 C7 数据框中的基因符号按基因集名称分组
geneSets = split(C7$gene_symbol, C7$gs_name) 
# geneSets = geneSets[c(1,3,5,7,9,11,13,15,17)]
geneSets = geneSets[c(10,15,20,25,80,120,500,600,701,855)]



# 3、选取上皮细胞signature进行打分
# 制作基因集文件
# devtools::install_github('immunogenomics/presto')
library(presto)#加快运行
celltype.markers <- FindAllMarkers(sce, assay = "RNA", only.pos = T)#差异分析
celltype.markers = celltype.markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
write.csv(celltype.markers,"Epi_all_allmarkers.csv",row.names = T)#保存
celltype.markers=read.csv("Epi_all_allmarkers.csv",header = T,row.names = 1)
top1000 = celltype.markers %>% group_by(cluster) %>% top_n(n = 1000, wt = avg_log2FC)
#拆分成基因集list
signatures <- split(top1000$gene, top1000$cluster)  




#获取data数据
exprMatrix =  GetAssayData(sce, assay = "RNA", layer = "data")
# 调用自定义函数plotGeneCount来绘制基因计数图，参数exprMatrix是表达矩阵
plotGeneCount(exprMatrix) # exprMatrix: 基因表达矩阵
dev.off()


# AUCell_buildRankings函数来构建细胞排名
cells_rankings <- AUCell_buildRankings(exprMatrix, splitByBlocks = TRUE)
#splitByBlocks=TRUE将数据分块以减少内存消耗
cells_rankings
# 显示细胞的排名结果
# 这样看更清楚
View(cells_rankings@assays@data@listData$ranking)



#AUCell分析
#为了应对大数据的计算以及去噪，只有每个细胞表达量最高的前5%或10%的基因参与计算
#这一比例可以通过'aucMaxRank'参数调整
cells_AUC <- AUCell_calcAUC(geneSets, #geneSets OR signatures
                            cells_rankings, 
                            aucMaxRank=nrow(cells_rankings)*0.1)
# geneSets: 基因集列表，包含要分析的基因集合
# cells_rankings: 细胞排名矩阵，之前使用AUCell_buildRankings函数生成
# aucMaxRank: 用于计算AUC的最大排名阈值，默认设置为细胞总数的5%，可以增加


par(mfrow=c(3,3)) #3X3绘图
# 给每个基因集分配活性阈值
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
# plotHist: 是否绘制直方图; assign: 是否自动分配细胞群
dev.off()

#热图
library(pheatmap)
aucMat <- getAUC(cells_AUC)
pheatmap(aucMat, show_colnames = F, scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))




#返回的数据结构包含软件选择的最佳阈值，其他阈值的返回结果(阈值参数与当前阈值下"on"的细胞数量)，软件返回的comment
cells_assignment$CAO_BLOOD_FLUZONE_AGE_05_14YO_7DY_DN$aucThr


warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)#把comment循环取出
warningMsg[which(warningMsg!="")]#查看非空的comment
#看看是否判定为正态分布？
### "The AUC might follow a normal distribution (random gene-set?).  
###  The global distribution overlaps the partial distributions. "

# 展示每个geneset的最佳阈值
sapply(cells_assignment, function(x) x$aucThr$selected)



#画出对应geneset直方图
geneSetName <- rownames(cells_AUC)[grep("CAO_BLOOD_FLUZONE_AGE_05_14YO_7DY_DN", rownames(cells_AUC))]
par(mfrow=c(1,1)) 
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=cells_assignment$BUCASAS_PBMC_FLUARIX_FLUVIRIN_CAUCASIAN_MALE_AGE_18_40YO_HIGH_RESPONDERS_1DY_3DY_POSITIVE_PREDICTIVE_OF_TITER$aucThr$selected)
abline(v=cells_assignment$`BUCASAS_PBMC_FLUARIX_FLUVIRIN_CAUCASIAN_MALE_AGE_18_40YO_HIGH_RESPONDERS_1DY_3DY_POSITIVE_PREDICTIVE_OF_TITER`$aucThr$selected)


#画出所有直方图
par(mfrow=c(3,3))
for(i in 1:length(cells_assignment)){  
  temp.name <- gsub(' (.*)','',names(cells_assignment)[i])  
  temp.cellname <- rownames(cells_AUC)[grep(temp.name, rownames(cells_AUC))]  
  AUCell_plotHist(cells_AUC[temp.cellname,], aucThr=cells_assignment[[i]]$aucThr$selected)  
  abline(v=cells_assignment[[i]]$aucThr$selected)
}




####AUcell##
#得到每个通路划分的阈值
selectedThresholds <- getThresholdSelected(cells_assignment)

# （表格和热图）
# 提取所有基因集的这些细胞，并将其转化为表格
cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)
# 转换为入射矩阵并绘制为直方图：
assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
assignmentMat[,1:2]



##进行on/off
selectedThresholds <- getThresholdSelected(cells_assignment)

dat = data.frame(sce@meta.data, sce@reductions$umap@cell.embeddings)
umap = data.frame(sce@reductions$umap@cell.embeddings)
umap = as.matrix(umap)
par(mfrow=c(1,1))
plot(umap, cex=.3)

assignmentTable[1:4,]

par(mfrow=c(3,3))
for(geneSetName in names(selectedThresholds)){
  colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(5)#处于off状态下的物种颜色
  colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(5)#处于on状态下的物种颜色
  # 将细胞按照on和off拆开
  passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
  if(sum(passThreshold) >0 )  {
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=5)], names(aucSplit[[1]])),
                   setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=5)], names(aucSplit[[2]])))
    plot(umap, main=geneSetName,
         sub="Pink/red cells pass the threshold",
         col=cellColor[rownames(umap)], pch=16)
  }
}





##挑选通路可视化on/off
par(mfrow=c(2,3))
###这里的tSNE是广义上的降维结果
AUCell_plotTSNE(tSNE=umap, exprMat=exprMatrix, 
                cellsAUC=cells_AUC[c(1,9),], thresholds=selectedThresholds[c(1,9)])

##挑选通路可视化其他方法ggplot
geneSet <- "ANDERSON_BLOOD_CN54GP140_ADJUVANTED_WITH_GLA_AF_AGE_18_45YO_LOW_IGM_RESPONDERS_6HY_1DY_UP"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
sce$AUC  <- aucs
#library(ggplot2)
#library(viridis)
par(mfrow=c(1,1))
library(ggraph)
ggplot(data.frame(sce@meta.data, sce@reductions$umap@cell.embeddings), aes(umap_1, umap_2, color=AUC)
) + geom_point( size=1.5
) + scale_color_viridis(option="A")  + theme_light(base_size = 15)+labs(title = "ANDERSON_BLOOD_CN54GP140_ADJUVANTED_WITH_GLA_AF_AGE_18_45YO_LOW_IGM_RESPONDERS_6HY_1DY_UP")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.2))


VlnPlot(sce,features = 'AUC', pt.size = 0, adjust = 2)

