Cellchat=readRDS("Cellchat.rds")

## 加载必要的包
library(CellChat)
library(ggplot2)
library(reshape2)
library(patchwork)

# 假设你已经构建好 CellChat 对象，命名为 cellchat_obj
# 并已运行 computeCommunProb、computeCommunProbPathway、aggregateNet

# 1. 提取所有信号通路和细胞类型
pathways <- Cellchat@netP$pathways
cells <- levels(Cellchat@idents)   # 按顺序的细胞类型名称

# 2. 初始化矩阵（通路 × 细胞类型）
outgoing_mat <- matrix(0, nrow = length(pathways), ncol = length(cells),
                       dimnames = list(pathways, cells))
incoming_mat <- matrix(0, nrow = length(pathways), ncol = length(cells),
                       dimnames = list(pathways, cells))

# 3. 填充矩阵
for (p in pathways) {
  # netP$prob 是四维数组（发送细胞 × 接收细胞 × 通路 × 配受体对）
  # 这里取该通路下所有配受体对的概率之和（已由 CellChat 汇总）
  prob_mat <- Cellchat@netP$prob[,, p]   # 矩阵：发送细胞 × 接收细胞
  
  outgoing <- rowSums(prob_mat, na.rm = TRUE)   # 外发强度：每个发送细胞的总发送量
  incoming <- colSums(prob_mat, na.rm = TRUE)   # 内向强度：每个接收细胞的总接收量
  
  # 确保顺序与细胞类型一致
  outgoing_mat[p, cells] <- outgoing[cells]
  incoming_mat[p, cells] <- incoming[cells]
}

# 总体信号强度 = 外发 + 内向
total_mat <- outgoing_mat + incoming_mat


####归一化
#  total_mat 是的总体信号强度矩阵，行为通路，列为细胞类型
# 按行（margin = 1）进行 min-max 归一化
normalize_row_minmax <- function(mat) {
  # 对每一行应用归一化
  mat_norm <- t(apply(mat, 1, function(row) {
    min_val <- min(row, na.rm = TRUE)
    max_val <- max(row, na.rm = TRUE)
    if (max_val == min_val) {
      # 如果该行所有值相等，则归一化为 0（或保持 0）
      return(rep(0, length(row)))
    } else {
      return((row - min_val) / (max_val - min_val))
    }
  }))
  # 恢复行名和列名
  rownames(mat_norm) <- rownames(mat)
  colnames(mat_norm) <- colnames(mat)
  return(mat_norm)
}

total_mat<- normalize_row_minmax(total_mat)
outgoing_mat<- normalize_row_minmax(outgoing_mat)
incoming_mat<- normalize_row_minmax(incoming_mat)


# 4. 绘制热图的函数（使用 ggplot2）
plot_heatmap <- function(mat, title) {
  df <- melt(mat, varnames = c("Pathway", "CellType"), value.name = "Strength")
  df$Pathway <- factor(df$Pathway, levels = rev(rownames(mat)))   # 倒序使热图顶部为第一个通路
  df$CellType <- factor(df$CellType, levels = colnames(mat))
  
  # 生成从紫到金的连续渐变色（100级）
  my_colors <- colorRampPalette(c("#2166ACFE",'#4393C3FF',"#92C5DEFF",'#D1E5F0FF',
                                  '#FDDBC7FF','#F4A582FF','#D6604DFF','#B2182BFF'))(100)
  
  p <- ggplot(df, aes(x = CellType, y = Pathway, fill = Strength)) +
    geom_tile(color = "White", size = 0.5) +
    scale_fill_gradientn(colors = my_colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = title, x = "", y = "")
  return(p)
}

# 5. 生成三张图
p_out <- plot_heatmap(outgoing_mat, "Outgoing signaling")
p_in  <- plot_heatmap(incoming_mat, "Incoming signaling")
p_total <- plot_heatmap(total_mat, "Total signaling")

# 6. 并列显示
combined <- p_out | p_in | p_total
combined

# 保存图片
ggsave("signaling_patterns.pdf", combined, width = 12, height = 8)
# 或者保存为文件
ggsave("cellchat_comparison.pdf", combined_plot, width = 15, height = 12)






