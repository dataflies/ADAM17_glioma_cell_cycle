rm(list = ls())

####16.CytoTRACE推断细胞分化程度、细胞干性####
# 安装CytoTRACE
# 包下载地址：https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz

#本地安装CytoTRACE
install.packages("CytoTRACE_0.3.3.tar.gz", repos = NULL, type="source")
# ERROR: dependencies 'HiClimR', 'ccaPP', 'nnls', 'egg' are not available for package 'CytoTRACE'
# * removing 'D:/R Language/R-4.2.1/library/CytoTRACE'
# Warning in install.packages :
#   安装程序包‘D:/桌面/CytoTRACE_0.3.3.tar.gz’时退出狀態的值不是0

# 报错说缺依赖包，缺啥包咱就就安装啥包，如果依赖包还报错就查看报错一般会出现包的链接，下载下了继续本地安装
install.packages("HiClimR")
install.packages("ccaPP")
install.packages("nnls")
install.packages("egg")
BiocManager::install("sva")

# 安装好了咱再来安装这个CytoTRACE
install.packages("CytoTRACE_0.3.3.tar.gz", repos = NULL, type="source")
#依赖包
library(CytoTRACE)
# 第一次library CytoTRACE 需要 Python 环境来运行一些依赖的功能。你可以选择安装 Miniconda
# 咱选择N就行，没必要安装python环境

# 加载上皮细胞数据
setwd("D:\\R language\\5.ADAM17\\NO2.scRNA\\CytoTRACE推断细胞分化程度和干性")
load("Epi_已注释.Rdata")
DimPlot(sce, reduction = "umap", group.by = "celltype",label = T,pt.size = 1,repel = T) 


# 将 RNA count 数据转换为矩阵
exp1 <- as.matrix(GetAssayData(sce, assay = "RNA", layer = "counts"))
# 过滤掉在少于 5 个细胞中表达的基因
exp1 <- exp1[apply(exp1 > 0, 1, sum) >= 5,]
# 使用 CytoTRACE 进行分析，设置 ncores = 1 表示使用一个 CPU 核心进行计算
results <- CytoTRACE(exp1, ncores = 1)
# 提取细胞类型注释信息，并将其转换为字符向量
phenot <- sce$celltype
phenot <- as.character(phenot)
# 将细胞类型注释的名字设置为元数据中的行名
names(phenot) <- rownames(sce@meta.data)
# 提取 UMAP 降维后的细胞嵌入信息
emb <- sce@reductions[["umap"]]@cell.embeddings

# 使用 CytoTRACE 结果、细胞类型注释和 UMAP 嵌入信息绘制 CytoTRACE 图，并将结果保存到指定目录
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir ="D:\\R language\\5.ADAM17\\NO2.scRNA\\CytoTRACE推断细胞分化程度和干性")
# 绘制 CytoTRACE 基因表达图，显示前 30 个基因，并将结果保存到指定目录
plotCytoGenes(results, numOfGenes = 30, outputDir = "D:\\R language\\5.ADAM17\\NO2.scRNA\\CytoTRACE推断细胞分化程度和干性")

