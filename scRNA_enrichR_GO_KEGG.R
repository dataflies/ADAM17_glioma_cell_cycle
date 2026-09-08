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
G_degs = FindMarkers( sce_pca, 
                      logfc.threshold = 0.25,
                      min.pct = 0.1, # 表达比例的阈值设置，小一点找出更多差异基因
                      only.pos = FALSE, #false既保留上调的基因，也保留下调的基因
                      ident.1 = c("G2", "G3", "G6"), ident.2 =c("G1", "G4", "G5") ) %>% # ident 1 vs 2
  mutate( gene = rownames(.) ) # 将行名新增为列，方便检索基因



# 自定义筛选
G_degs_fil = G_degs %>% 
  filter( pct.1 > 0.2 & p_val_adj < 0.05 ) %>% # &：和；|：或。 
  filter( abs( avg_log2FC ) > 0.3 ) #%>% # 按 logFC 绝对值筛选
# filter( avg_log2FC > 0 ) # 按 logFC 正负筛选
write.csv(G_degs_fil,file = "G_degs_fil")
save(G_degs_fil,file ='G_degs_fil.Rdata')



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
dotplot(cluster3_ego, showCategory=10,title="ADAM17_H vs ADAM17_L GO")
barplot(cluster3_ego, split ="ONTOLOGY", font.size =10,  showCategory=10) +facet_grid(ONTOLOGY ~ ., space ="free_y",scales ="free_y")

#KEGG富集
cluster3_ekg<- enrichKEGG(gene =cluster3_de,
                          organism ='hsa',#物种Homo sapiens
                          pvalueCutoff =0.05,
                          qvalueCutoff =0.05,
                          minGSSize =10,
                          maxGSSize =500)
cluster3_ekg <- enrichKEGG(gene= cluster3_de, organism = "hsa",pvalueCutoff = 0.05)
head(cluster3_ekg)
#气泡图
barplot(cluster3_ekg)
dotplot(cluster3_ekg, showCategory=10,title="ADAM17_H vs ADAM17_L KEGG")


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





rm(list = ls())
