rm(list = ls())


## 拟时序分析 ------------------------------
library(Seurat)

#没有monocle要先安装 BiocManager::install()
#BiocManager::install('monocle',update = F,ask = F)

#没有的包先安装
library(BiocGenerics)
library(monocle)
library(tidyverse)
library(patchwork)

setwd("D:\\R language\\7.zhongshan\\空转")
scRNA=readRDS('./scRNA_freash_anno_1v1_scoring.RDS')

scRNAsub=subset(scRNA,celltype=='Epithelial_cells')
#scRNAsub=subset(scRNAsub,tissue_type=='Tumor')

#V5提取data的格式
data=GetAssayData(object=scRNAsub, slot="data", assay=DefaultAssay(scRNAsub))

data <- as(data, 'sparseMatrix')  #转化为稀松矩阵
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
## 以下代码一律不得修改 ！
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())



mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

##使用monocle选择的高变基因，不修改
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
plot_ordering_genes(mycds)

#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序，报错请用4.1.3的R，并重装monocle
mycds <- orderCells(mycds)

## 查看关键差异基因在拟时序中所处的时间
load('G_degs_fil.Rdata')

df2 = G_degs_fil
my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[rownames(df2),],
                                                 # num_clusters = 2, # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)

my_pseudotime_cluster 


#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
plot1

plot4 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
plot4

##合并出图
plotc <- plot1|plot4
plotc

DimPlot(scRNAsub,group.by = 'seurat_clusters',label = T)




###流程二####
rm(list = ls())

# 细胞拟时序分析
library(Seurat)
#没有安装包的需要先安装这些包
#install.packages(c("magrittr", "RColorBrewer", "reshape2", "ggsci", "ggpubr", "data.table"))
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("Biobase")
#BiocManager::install("monocle")
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(reshape2)
library(Biobase)
library(ggsci)
library(ggpubr)
library(data.table)
#.libPaths()
library(igraph)
library(monocle)
#install.packages("https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz",repos=NULL)
#remove.packages('igraph')删除电脑里固有的包以便安装新的包


# 加载Biobase包
#unloadNamespace("igraph")
set.seed(12345)
#saveRDS(cds,'cds.rds')
setwd("D:\\R language\\5.ADAM17\\NO2.scRNA\\GSE117891\\拟时序")
#load("...\\sce")

Idents(sce) <- sce$celltype#将celltype赋值给Idents属性；已经注释好的细胞类型
DimPlot(sce,label = T)
dev.off()
####创建monocle对象####
unique(sce$celltype)
# 01选择需要构建细胞分化轨迹的细胞类型（subset提取感兴趣的名字，上面已经准备好名字了）
seurat=subset(sce,idents = c("G2", "G3", "G4"))
#saveRDS(seurat,'seurat.rds')
unique(seurat$celltype)
#[1] "G1" "G3" "G2" "G6" "G5"
table(seurat$celltype)
# G1  G2  G3  G4  G5  G6 
# 900 654 562 482 226 210 

# 02构建monocle独有的celldata对象
data=GetAssayData(object=seurat, slot="data", assay=DefaultAssay(seurat))

expr_matrix=data#使用counts表达值(第一个准备文件：基因count矩阵)
sample_sheet<-seurat@meta.data#将实验信息赋值新变量（第二个准备文件：细胞表型信息-metadata）
gene_annotation=data.frame(gene_short_name=rownames(seurat))#构建一个含有基因名字的数据框
rownames(gene_annotation)=rownames(seurat)#（第三个文件：基因名为行名的基因名数据框）
pd <- new("AnnotatedDataFrame", data = sample_sheet)#将实验信息变量转化为monocel可以接收的对象
fd <- new("AnnotatedDataFrame", data = gene_annotation)#将基因注释变量转化为monocle可以接收的对象

class(expr_matrix)


cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd,expressionFamily=negbinomial.size())#创建一个monocle的对象
cds #cellData对象；monocle独有


# 03相当于归一化Size factors帮助标准化基因在细胞间的差异, "dispersion"值帮助后面的差异表达分析执行。
cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds)


# 04数据降维，特征基因的选择
# 过滤低表达基因；差异分析获得top1000排序基因（可视化用来排序的基因）
###注意--differentialGeneTest这一步有些久
cds <- detectGenes(cds, min_expr = 0.1)#统计，过滤低表达基因（表达量过滤低表达基因）
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))#（数量过滤低表达基因）

######monocle差异分析选择拟时序相关基因####
# 选择一
diff_celltype <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~celltype",cores=4)#差异分析
head(diff_celltype)#查看前几行数据
write.csv(diff_celltype,'DegForCellOrdering_5.csv')


# #选择二、使用monocle选择的高变基因
# disp_table <- dispersionTable(cds)
# disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
# cds <- setOrderingFilter(cds, disp.genes)
# plot_ordering_genes(cds)#以上是选择monocle选择的高变基因进行轨迹构建，上述为高变基因计算方法，monocle在选择做轨迹的基因上共提供四种方法
# #图中黑点就代表选出的高变基因(即离散度高的基因)


# 对细胞排序，使用q值排序前1000个基因
diff_celltype<- diff_celltype[order(diff_celltype$qval),]#将差异基因按照q值进行升序排列
ordering_genes <- row.names(diff_celltype[1:1100,]) #选择top1000个显著基因用做后续排序，这里可修改
ordering_genes <- row.names(diff_celltype[diff_celltype$qval < 0.01,]) 
cds <- setOrderingFilter(cds,ordering_genes = ordering_genes) #设置排序基因
plot_ordering_genes(cds)#可视化用来排序的基因（有着更高的表达与离散）
ggsave(filename = 'monocle2_ordering_gene.pdf') 


library(igraph)
package.version()
# 降维；绘制细胞轨迹（绘制state与细胞类型的轨迹图）
cds <- reduceDimension(cds, method = 'DDRTree')#降维
cds <- orderCells(cds)#排序
#packageVersion('igraph')
##如果这里报错，先保存数据
saveRDS(cds,'cds_G2+G5.rds')

cds<-readRDS('cds.rds')
####运行完orderCells，保存数据####
#save.image(file = "orderCells_cds_before.RData")  # 保存所有对象到my_workspace.RData
#load("orderCells_cds_done.RData")



p1=plot_cell_trajectory(cds, color_by = "State") +#可视化鉴定的State
  theme(text = element_text(size = 18))  # 设置字体大小为 18
p1
ggsave(p1,filename = 'NoG4_monocle2_state_trajectory.pdf', width = 12, height = 9)
p2=plot_cell_trajectory(cds, color_by = "celltype") +#可视化细胞类型
  theme(text = element_text(size = 18))  # 设置字体大小为 18

p2
ggsave(p2,filename = 'NoG4_monocle2_celltype_trajectory.pdf', width = 10, height = 6)
# 保存排序所需要的拟时间值与细胞state值

#4.按发育时间分类画轨迹图，可以与state的图对应上
p3=plot_cell_trajectory(cds, color_by = "Pseudotime")+
  theme(text = element_text(size = 18))  # 设置字体大小为 18
p3
ggsave(p3,filename = 'NoG2+G3_monocle2_Pseudotime.pdf', width = 11, height = 6) 


library(patchwork)
# 水平排列图形（并排）
plotc <- p1 + p2 + p3
plotc
ggsave(plotc,filename = 'G3+G4_monocle2_ALL.pdf', width = 10, height = 6)

#指定细胞分化轨迹的起始点，即设置root cell。
#自定义了一个名为GM_state的函数，函数主要计算包含特定类型细胞最多的state，并返回state编号
#即特定细胞类型最多的state为指定起点
GM_state <- function(cds){
if (length(unique(pData(cds)$State)) > 1){
 T0_counts <- table(pData(cds)$State, pData(cds)$celltype)[,"G4"]
  return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
} else {return (1)}
}
#有时间点采样，就将imm_cell_type改成hours；0时刻为起始状态
cds <- orderCells(cds, root_state = GM_state(cds)) #指定分化起始点，排序细胞

#指定起点后，绘制细胞伪时间轨迹图与分支图
p4=plot_cell_trajectory(cds, color_by = "Pseudotime")#按pseudotime（伪时间）来可视化细胞分化轨迹
p4
ggsave(p2,filename = 'monocle2_Pseudotime.pdf', width = 12, height = 9) 

p5=plot_cell_trajectory(cds, color_by = "celltype") +#可视化细胞类型
  theme(text = element_text(size = 18))  # 设置字体大小为 18

p5

p2+p3
ggsave('monocle2_celltype_Pseudotime.pdf', width = 20, height = 9)


plot_cell_trajectory(cds, color_by = "celltype")+
  facet_wrap(~celltype, nrow = 2)#facet_wrap函数可以把每个state单独highlight
ggsave(filename = 'G2+G3+G5_monocle2_facet_celltype.pdf', width = 20, height = 14) 

####按celltype分类画轨迹图，并按分支点分类，每个celltype的细胞分别标出####
p4=plot_cell_trajectory(cds, color_by = "celltype")+
  facet_wrap(~celltype, nrow = 5) +
  theme(text = element_text(size = 25))  # 设置字体大小为 18
p4
ggsave(p4,file='monocle2_celltype分开.pdf', width = 12, height = 18)

#树形图
p5 <- plot_complex_cell_trajectory(cds,x = 1, y = 2,
                                   color_by = "celltype")

p5

#细胞密度图
p6 <- ggplot(pData(cds),aes(Pseudotime,colour = celltype,fill = celltype))+
  geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic()

p6
ggsave("G3+G4_Trajectory.Density.pdf",plot=p6,width=10,height=6.5)




#### 可视化特定基因在各个state的表达(可辅助判断细胞起点)####
#注意：可以可视化一些增殖基因在哪些state高表达（一般干性细胞高表达或细胞周期marker）
#然后可以假设其为分化起点。（当然也可以看自身感兴趣的）
blast_genes <- row.names(subset(fData(cds),gene_short_name %in% c('ADAM17','EGFR','HRAS','IMPDH','TP53','MAPK8')))#获取特定基因向量
plot_genes_jitter(cds[blast_genes,],grouping = "State",min_expr = 0.1)
ggsave(filename = 'G2+G3+G4_monocle2_gene_state.pdf',width = 12, height = 18)

blast_genes <- row.names(subset(fData(cds),gene_short_name %in% c("ACTA2","TAGLN","S100A4")))#获取特定基因向量
plot_genes_jitter(cds[blast_genes,],grouping = "State",min_expr = 0.1)#可视化基因表达，散点图
ggsave(filename = 'monocle2_specific_gene_state.pdf',width = 12, height = 18) 
#检查每个阶段celltype的基因含量
plot_genes_jitter(cds[blast_genes,],grouping = "celltype",min_expr = 0.1)#可视化基因表达，散点图

ggsave(filename = 'monocle2_specific_gene_celltype.pdf',width = 20, height = 12) 




#### 3.鉴定伪时间相关的基因，即分化过程相关 ----------------------------------------------------
#鉴定伪时间相关的基因，即分化过程（state）相关的差异基因（可以指定计算个数，20可为100）
diff_test_res <- differentialGeneTest(cds,cores=4,fullModelFormulaStr = "~sm.ns(Pseudotime)")#鉴定伪时间相关的基因，即分化过程相关
write.csv(diff_test_res,file='diff_test_res.csv') 
## 注意可以将top3g改成任何自身感兴趣的基因，进行可视化
top5g=rownames(diff_test_res[order(diff_test_res$qval),])[1:5] # top 5个显著基因
interest_gene =c('ADAM17','EGFR','HRAS','TP53','MAPK8')


plot_genes_in_pseudotime(cds[interest_gene,], color_by="celltype", ncol = 1) +# top 5个显著基因，可视化，pseudotime作为横坐标
  theme(text = element_text(size = 15))  # 设置字体大小为 18
ggsave(filename = 'G2+G3+G4_monocle2_significant_gene.pdf',width = 12, height = 18) 
write.csv(top5g,file='monocle_deg_top5.csv') #可以top10、top100按需改动数字即可


####伪时间相关基因热图####
library(viridis)
#伪时间基因排序
diff_test_res = diff_test_res[order(diff_test_res$qval),]
diff_test_res = diff_test_res[1:40,]

plot_pseudotime_heatmap(cds[interest_gene,], #挑选的20个计算了分化差异的基因
                        num_clusters = 8,
                        cores = 4, #指定线程数
                        show_rownames = T)#对伪时间相关的基因绘制热图，查看具备相似表达模式的基因簇
ggsave(filename = 'NoG2+G4_monocle2_gene_pheatmap.pdf',width = 18, height = 15) 

#颜色
plot_pseudotime_heatmap(cds[rownames(diff_test_res),],
                        num_clusters = 8,
                        cores = 4,
                        show_rownames = T,
                        hmcols = colorRampPalette(viridis(4))(1000))#对伪时间相关的基因绘制热图，查看具备相似表达模式的基因簇



# 4.鉴定分支依赖基因--------------------------------------------------------------
### 分支比较后，分支依赖基因的鉴定（热图重现）
# 基因簇按照自身需求，设置；图形上侧显示两个分支（朝不同方向分化）

remove(BEAM_res)
BEAM_res <- BEAM(cds[interest_gene,], branch_point = 1, cores = 4, progenitor_method = "duplicate")#分支依赖基因鉴定，这里只选择了10个基因进行演示
BEAM_res = BEAM_res[order(diff_test_res$qval),]
BEAM_res = BEAM_res[1:100,]

#?BEAM #鉴定分支依赖的基因
plot_genes_branched_heatmap(cds[interest_gene,],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T)


# 确保 expressed_genes 是 cds 中的有效基因
expressed_genes <- rownames(cds)[rowSums(exprs(cds)) > 0]

# 运行 BEAM 分析
BEAM_res <- BEAM(cds[interest_gene, ], branch_point = 1, cores = 4, progenitor_method = "duplicate")

# 按 qval 排序并取前 100 个基因
BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
BEAM_res <- BEAM_res[1:100, ]

# 检查基因名称是否在 cds 中
valid_genes <- rownames(BEAM_res)[rownames(BEAM_res) %in% fData(cds)$gene_short_name]

# 绘制热图
plot_genes_branched_heatmap(cds[interest_gene, ],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 4,
                            use_gene_short_name = TRUE,
                            show_rownames = TRUE)
pp

http://127.0.0.1:19135/graphics/plot_zoom_png?width=1188&height=865
ggsave(filename = 'top100.pdf',width = 10, height = 20) 


# 提取目标基因的表达数据
genes_to_plot <- c("JUN", "RUNX3", "CSF3R", "PDE4B")

# 绘制基因表达随拟时序变化的趋势图
p_genes_pseudotime <- plot_genes_in_pseudotime(cds[interest_gene, ], 
                                               color_by = "State") +  # 按 State 着色
  theme(text = element_text(size = 18))  # 设置字体大小为 18

# 显示图像
p_genes_pseudotime

# 保存图像
ggsave(p_genes_pseudotime, filename = 'monocle2_genes_pseudotime_trajectory.pdf', width = 12, height = 9)

pData(cds)$ADAM17=log2(exprs(cds)["ADAM17",]+1)
p11 = plot_cell_trajectory(cds,color_by = "ADAM17")+
  scale_color_continuous(type = "viridis")

p11

pData(cds)$EGFR=log2(exprs(cds)["EGFR",]+1)
p12 = plot_cell_trajectory(cds,color_by = "EGFR")+
  scale_color_continuous(type = "viridis")
p12


pData(cds)$PIK3CA=log2(exprs(cds)["",]+1)
p13 = plot_cell_trajectory(cds,color_by = "PIK3CA")+
  scale_color_continuous(type = "viridis")
p13

pData(cds)$AKT1=log2(exprs(cds)["AKT1",]+1)
p14 = plot_cell_trajectory(cds,color_by = "AKT1")+
  scale_color_continuous(type = "viridis")
p14

ggsave("Trajectory.expression.PDF",plot = p11|p12|p13|p14,width = 14,height = 6.5)


pData(cds)$COL16A1=log2(exprs(cds)["COL16A1",]+1)
p15 = plot_cell_trajectory(cds,color_by = "COL16A1")+
  scale_color_continuous(type = "viridis")

pData(cds)$ECE1=log2(exprs(cds)["ECE1",]+1)
p16 = plot_cell_trajectory(cds,color_by = " ECE1")+
  scale_color_continuous(type = "viridis")

pData(cds)$ CDKN2C=log2(exprs(cds)["CDKN2C",]+1)
p17 = plot_cell_trajectory(cds,color_by = " CDKN2C")+
  scale_color_continuous(type = "viridis")

pData(cds)$VCAM1=log2(exprs(cds)["VCAM1",]+1)
p18 = plot_cell_trajectory(cds,color_by = "VCAM1")+
  scale_color_continuous(type = "viridis")

ggsave("Trajectory.expression.fate1.PDF",plot = p15|p16|p17|p18,width = 14,height = 6.5)


rownames(cds)
rownames(BEAM_res)



pdf("SSS")      
dev.off()
#绘制分支热图
ggsave(filename = 'monocle2_beam_heatmap.pdf',width = 15, height = 12) 

### 选取BEAM的三个基因做例证(感兴趣基因可以替换rownames(BEAM_res)[]可视化)
plot_genes_branched_pseudotime(cds[interest_gene,],
                               branch_point = 1,
                               color_by = "celltype",ncol = 1)#绘制分支pseudotime散点图
ggsave(filename = 'NoG4+G3_monocle2_beam_gene.pdf',width = 15, height = 12) 
plot_genes_branched_pseudotime(cds[interest_gene,],
                               branch_point = 1,
                               color_by = "State",ncol = 1)#绘制分支pseudotime散点图
ggsave(filename = 'monocle2_beam2_gene.pdf',width = 15, height = 12) 





























