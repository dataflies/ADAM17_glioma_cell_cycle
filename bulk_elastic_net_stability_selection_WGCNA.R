# ==================== 1. 加载必要的包 ====================
library(tidyverse)
library(GSVA)         # 提供 ssgsea 算法
library(GSEABase)     # 读取 GMT 文件
library(ggplot2)
library(ggpubr)       # 方便添加统计标签

# ==================== 2. 读取表达矩阵并预处理 ====================
# 读取数据（注意：CGGA 文件是制表符分隔，第一列为基因名）
expr <- read.delim("E:/ADAM17 data/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",
                   row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

cat("原始表达矩阵维度：", dim(expr), "\n")
# 去除重复列（样本名重复）
expr <- expr[, !duplicated(colnames(expr))]
# 转换为数值矩阵（强制转换，非数值变 NA）
expr <- as.matrix(expr)
mode(expr) <- "numeric"
# 去除全零的基因（行和为零）
expr <- expr[rowSums(expr, na.rm = TRUE) > 0, ]
# log2(x+1) 转换
expr_log <- log2(expr + 1)
cat("log转换后维度：", dim(expr_log), "\n")

# ==================== 3. 读取 GMT 文件并提取细胞周期相关基因 ====================
# 自定义读取 GMT 的函数（返回 GeneSetCollection）
read_gmt <- function(gmt_path) {
  gsc <- getGmt(gmt_path, geneIdType = SymbolIdentifier())
  return(gsc)
}

# 读取两个 GMT 文件
g2m_gsc <- read_gmt("E:/ADAM17 data/HALLMARK_G2M_CHECKPOINT.v2025.1.Hs.gmt")
e2f_gsc <- read_gmt("E:/ADAM17 data/HALLMARK_E2F_TARGETS.v2025.1.Hs.gmt")

# 提取基因名（取第一个基因集的基因列表）
g2m_genes <- geneIds(g2m_gsc[[1]])
e2f_genes <- geneIds(e2f_gsc[[1]])

# 合并并去重
cell_cycle_genes <- unique(c(g2m_genes, e2f_genes))
cat("细胞周期基因集总基因数：", length(cell_cycle_genes), "\n")

# 与表达矩阵中的基因取交集
expr_genes <- rownames(expr_log)
common_genes <- intersect(cell_cycle_genes, expr_genes)
cat("用于 ssGSEA 的基因数：", length(common_genes), "\n")

# 构建用于 GSVA 的基因列表（一个基因集）
gene_list <- list(Cell_Cycle = common_genes)

# ==================== 4. 自定义 ssGSEA 函数（基于秩次） ====================
# 输入：样本×基因矩阵（行=样本，列=基因），基因集（字符向量）
# 输出：每个样本的 ssGSEA 分数（标准化后）
custom_ssgsea <- function(expr_matrix, gene_set) {
  # expr_matrix: 样本 × 基因（数值矩阵，已 log 转换）
  # gene_set: 基因名向量
  gene_set <- intersect(gene_set, colnames(expr_matrix))
  if (length(gene_set) < 5) stop("基因集过小，无法计算")
  
  # 对每个基因（列）计算秩次（表达量高的秩次小）
  ranks <- apply(expr_matrix, 2, function(x) rank(-x, ties.method = "average"))
  
  # 对每个样本（行），计算基因集中基因的秩次之和 / 基因集大小
  scores <- apply(ranks, 1, function(row) {
    sum(row[gene_set]) / length(gene_set)
  })
  
  # Z-score 标准化
  scores <- scale(scores)
  return(as.numeric(scores))
}

# 准备数据：将 expr_log 转置为 样本 × 基因
expr_sample_gene <- t(expr_log)   # 行=样本，列=基因
cat("转置后矩阵维度：", dim(expr_sample_gene), "\n")  # 应为 693 × 23987

# 计算细胞周期分数
cellcycle_scores <- custom_ssgsea(expr_sample_gene, common_genes)
names(cellcycle_scores) <- rownames(expr_sample_gene)

# 构建评分数据框
cellcycle_df <- data.frame(
  Sample = names(cellcycle_scores),
  CellCycle_score = cellcycle_scores,
  row.names = NULL
)
head(cellcycle_df)

# ==================== 5. 按中位数分组 ====================
median_score <- median(cellcycle_df$CellCycle_score)
cellcycle_df$Group_median <- ifelse(cellcycle_df$CellCycle_score >= median_score,
                                    "High", "Low")
table(cellcycle_df$Group_median)

# ==================== 6. 提取 ADAM17 表达量 ====================
# 从原始 log 矩阵（基因×样本）中提取
if ("ADAM17" %in% rownames(expr_log)) {
  adam17_expr <- expr_log["ADAM17", ]
} else if ("Adam17" %in% rownames(expr_log)) {
  adam17_expr <- expr_log["Adam17", ]
} else if ("TACE" %in% rownames(expr_log)) {
  adam17_expr <- expr_log["TACE", ]
} else {
  stop("ADAM17 基因未在表达矩阵中找到，请检查基因名")
}

# 合并数据框
df_corr <- data.frame(
  Sample = cellcycle_df$Sample,
  CellCycle_score = cellcycle_df$CellCycle_score,
  ADAM17_expr = adam17_expr[match(cellcycle_df$Sample, names(adam17_expr))],
  Group = cellcycle_df$Group_median
)
df_corr <- na.omit(df_corr)
cat("有效样本数：", nrow(df_corr), "\n")

# ==================== 7. 相关性分析（Pearson & Spearman） ====================
pearson_res <- cor.test(df_corr$CellCycle_score, df_corr$ADAM17_expr, method = "pearson")
spearman_res <- cor.test(df_corr$CellCycle_score, df_corr$ADAM17_expr, method = "spearman")

cat("\n========== 相关性结果 ==========\n")
cat(sprintf("Pearson r = %.4f, P = %.4e\n", pearson_res$estimate, pearson_res$p.value))
cat(sprintf("Spearman rho = %.4f, P = %.4e\n", spearman_res$estimate, spearman_res$p.value))

# ==================== 8. 散点图（带回归线） ====================
p1 <- ggplot(df_corr, aes(x = CellCycle_score, y = ADAM17_expr)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  labs(title = "Cell Cycle Score vs ADAM17 expression",
       x = "Cell Cycle Score (ssGSEA, scaled)",
       y = "ADAM17 expression (log2(TPM+1))") +
  theme_minimal()
print(p1)

# ==================== 9. 组间比较（High vs Low） ====================
wilcox_res <- wilcox.test(ADAM17_expr ~ Group, data = df_corr)
cat("\nMann-Whitney U 检验 (High vs Low): P =", format(wilcox_res$p.value, scientific = TRUE), "\n")

p2 <- ggplot(df_corr, aes(x = Group, y = ADAM17_expr, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c("High" = "#F8766D", "Low" = "#00BFC4")) +
  labs(title = "ADAM17 expression by Cell Cycle group",
       y = "ADAM17 expression (log2(TPM+1))") +
  theme_minimal()
print(p2)

# ==================== 10. 保存结果（可选） ====================
write.csv(df_corr, "cellcycle_ADAM17_correlation.csv", row.names = FALSE)
ggsave("scatter_cellcycle_ADAM17.png", p1, width = 6, height = 5, dpi = 300)
ggsave("boxplot_ADAM17_by_group.png", p2, width = 5, height = 4, dpi = 300)

cat("\n分析完成！结果已保存。\n")

# ==================== 11. 遍历所有基因与细胞周期评分的相关性 ====================
# 假设已有：
#   expr_log : 基因 × 样本 矩阵（log2(TPM+1)）
#   cellcycle_df : 包含 Sample 和 CellCycle_score 列
# 确保样本顺序一致
common_samples <- intersect(colnames(expr_log), cellcycle_df$Sample)
expr_log_sub <- expr_log[, common_samples, drop = FALSE]
cellcycle_sub <- cellcycle_df$CellCycle_score[match(common_samples, cellcycle_df$Sample)]

# 转置为 样本 × 基因
expr_sample_gene <- t(expr_log_sub)   # 行=样本，列=基因

# 预定义结果存储
gene_names <- colnames(expr_sample_gene)
n_genes <- length(gene_names)

# 初始化结果矩阵
pearson_r <- numeric(n_genes)
pearson_p <- numeric(n_genes)
spearman_rho <- numeric(n_genes)
spearman_p <- numeric(n_genes)

# 循环计算（基因数较多时可能需要几分钟，可加进度条）
cat("开始计算所有基因与细胞周期评分的相关性...\n")
pb <- txtProgressBar(min = 0, max = n_genes, style = 3)

for (i in seq_len(n_genes)) {
  gene_expr <- expr_sample_gene[, i]
  # Pearson
  pt <- cor.test(cellcycle_sub, gene_expr, method = "pearson")
  pearson_r[i] <- pt$estimate
  pearson_p[i] <- pt$p.value
  # Spearman
  st <- cor.test(cellcycle_sub, gene_expr, method = "spearman")
  spearman_rho[i] <- st$estimate
  spearman_p[i] <- st$p.value
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 构建结果数据框
corr_results <- data.frame(
  Gene = gene_names,
  Pearson_r = pearson_r,
  Pearson_p = pearson_p,
  Spearman_rho = spearman_rho,
  Spearman_p = spearman_p,
  stringsAsFactors = FALSE
)

# 可选：添加 FDR 校正
corr_results$Pearson_FDR <- p.adjust(corr_results$Pearson_p, method = "fdr")
corr_results$Spearman_FDR <- p.adjust(corr_results$Spearman_p, method = "fdr")

# ==================== 12. 按 |Pearson r| > 0.4 筛选基因 ====================
selected_genes <- corr_results[abs(corr_results$Pearson_r) > 0.4, ]
cat(sprintf("\n满足 |Pearson r| > 0.4 的基因数：%d\n", nrow(selected_genes)))

# 也可以按 Spearman |rho| > 0.4 筛选（可选）
selected_genes_spearman <- corr_results[abs(corr_results$Spearman_rho) > 0.4, ]
cat(sprintf("满足 |Spearman rho| > 0.4 的基因数：%d\n", nrow(selected_genes_spearman)))

# ==================== 13. 保存结果 ====================
write.csv(corr_results, "all_genes_correlation_with_cellcycle.csv", row.names = FALSE)
write.csv(selected_genes, "selected_genes_pearson_r_gt_0.4.csv", row.names = FALSE)

# 可选：仅输出基因名列表（一行一个）
writeLines(selected_genes$Gene, "selected_gene_list_pearson_r_gt_0.4.txt")

cat("\n分析完成！结果已保存为：\n",
    "- all_genes_correlation_with_cellcycle.csv (所有基因)\n",
    "- selected_genes_pearson_r_gt_0.4.csv (筛选后详细表格)\n",
    "- selected_gene_list_pearson_r_gt_0.4.txt (仅基因名)\n")

# ==================== 修正版：弹性网络回归（避免过拟合 + 稳定性评估） ====================
library(glmnet)

# 假设已存在：
#   selected_genes : 满足 |Pearson r|>0.4 的基因数据框（至少包含 Gene 列）
#   expr_log : 基因×样本矩阵
#   cellcycle_df : 包含 Sample 和 CellCycle_score

# ---------- 1. 准备数据：过滤近零方差基因 ----------
genes_input <- selected_genes$Gene
expr_sub <- expr_log[genes_input, ]  # 基因 × 样本
# 计算每个基因的表达方差（基于 log 表达值）
gene_vars <- apply(expr_sub, 1, var, na.rm = TRUE)
# 保留方差大于 0.1 的基因（阈值可根据数据分布调整，例如方差 > 0.1 或 > 0.05）
min_var <- 0.1
high_var_genes <- names(gene_vars[gene_vars > min_var])
cat(sprintf("过滤低方差基因：输入 %d 个，保留 %d 个（方差 > %.2f）\n", 
            length(genes_input), length(high_var_genes), min_var))

if (length(high_var_genes) < 2) stop("高方差基因不足，无法进行弹性网络")

# 重新准备 X, Y
common_samples <- intersect(colnames(expr_log), cellcycle_df$Sample)
X <- t(expr_log[high_var_genes, common_samples])  # 样本 × 基因
Y <- cellcycle_df$CellCycle_score[match(common_samples, cellcycle_df$Sample)]

# ---------- 2. 确定最佳 alpha（基于 lambda.1se 的误差）----------
alpha_candidates <- seq(0.1, 1, by = 0.1)  # 也可包含 0，但 0 是岭回归，通常弹性网络 alpha>0
set.seed(123)  # 用于固定初次的 alpha 选择，但后续稳定性测试会更换种子

# 函数：给定 alpha，返回交叉验证的 lambda.1se 对应的平均误差
get_cv_error_1se <- function(alpha, X, Y, nfolds = 10) {
  cv_fit <- cv.glmnet(X, Y, alpha = alpha, nfolds = nfolds, standardize = TRUE)
  # 取 lambda.1se 对应的误差
  idx_1se <- which(cv_fit$lambda == cv_fit$lambda.1se)
  return(cv_fit$cvm[idx_1se])
}

cv_errors_1se <- sapply(alpha_candidates, function(a) {
  get_cv_error_1se(a, X, Y)
})
best_alpha <- alpha_candidates[which.min(cv_errors_1se)]
cat(sprintf("最优 alpha (基于 lambda.1se 误差)：%.1f\n", best_alpha))

# ---------- 3. 多次重复交叉验证（稳定性评估）----------
n_repeats <- 50   # 推荐 50-100 次
set.seed(2025)    # 全局种子，但每次重复使用不同随机种子

# 存储每次重复中非零系数的基因名
selected_genes_list <- vector("list", length = n_repeats)

for (rep in 1:n_repeats) {
  # 每次使用不同的随机种子（基于 rep）
  seed_rep <- 123 + rep * 10
  set.seed(seed_rep)
  
  # 交叉验证（固定 alpha = best_alpha）
  cv_fit_rep <- cv.glmnet(X, Y, alpha = best_alpha, nfolds = 10, standardize = TRUE)
  # 使用 lambda.1se 拟合最终模型
  final_model_rep <- glmnet(X, Y, alpha = best_alpha, lambda = cv_fit_rep$lambda.1se, standardize = TRUE)
  coef_rep <- as.matrix(coef(final_model_rep))
  # 提取非零系数的基因（去掉截距项）
  nonzero_genes_rep <- rownames(coef_rep)[which(coef_rep != 0)][-1]
  selected_genes_list[[rep]] <- nonzero_genes_rep
  
  if (rep %% 10 == 0) cat(sprintf("完成 %d/%d 次重复\n", rep, n_repeats))
}

# 计算每个基因被选中的频率
all_selected_genes <- unique(unlist(selected_genes_list))
frequency <- sapply(all_selected_genes, function(g) {
  sum(sapply(selected_genes_list, function(x) g %in% x)) / n_repeats
})
freq_df <- data.frame(Gene = all_selected_genes, Frequency = frequency)
freq_df <- freq_df[order(-freq_df$Frequency), ]

# 筛选稳定基因（例如频率 >= 0.8，即 50 次中至少出现 40 次）
stable_threshold <- 0.8
stable_genes <- freq_df$Gene[freq_df$Frequency >= stable_threshold]
cat(sprintf("\n稳定性筛选：%d 个基因出现频率 >= %.1f\n", length(stable_genes), stable_threshold))

# ---------- 4. 最终模型（使用全部数据 + 稳定基因，再次拟合以获取系数）----------
if (length(stable_genes) > 0) {
  X_final <- X[, stable_genes, drop = FALSE]
  final_cv <- cv.glmnet(X_final, Y, alpha = best_alpha, nfolds = 10, standardize = TRUE)
  final_model <- glmnet(X_final, Y, alpha = best_alpha, lambda = final_cv$lambda.1se, standardize = TRUE)
  coef_final <- as.matrix(coef(final_model))
  coef_df <- data.frame(
    Gene = rownames(coef_final)[-1],
    Coefficient = as.numeric(coef_final[-1])
  )
  coef_df <- coef_df[order(-abs(coef_df$Coefficient)), ]
  
  # 保存结果
  write.csv(freq_df, "elasticnet_stability_frequencies.csv", row.names = FALSE)
  writeLines(stable_genes, "elasticnet_stable_genes.txt")
  write.csv(coef_df, "elasticnet_final_coefficients.csv", row.names = FALSE)
  
  cat("\n最终稳定基因列表已保存。\n")
} else {
  cat("警告：没有基因通过稳定性阈值，请降低阈值或检查数据。\n")
}

# ---------- 5. 可视化（可选）----------
# 频率柱状图
library(ggplot2)
p_freq <- ggplot(freq_df[1:min(30, nrow(freq_df)), ], aes(x = reorder(Gene, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Elastic Net 特征选择稳定性",
       x = "Gene", y = "Selection Frequency (50 repeats)") +
  theme_minimal()
ggsave("elasticnet_stability_plot.png", p_freq, width = 8, height = 6)

# 最终模型系数图
if (nrow(coef_df) > 0) {
  p_coef <- ggplot(coef_df[1:min(20, nrow(coef_df)), ], aes(x = reorder(Gene, Coefficient), y = Coefficient)) +
    geom_bar(stat = "identity", fill = "coral") +
    coord_flip() +
    labs(title = "Final Elastic Net Coefficients (lambda.1se, stable genes)",
         x = "Gene", y = "Coefficient") +
    theme_minimal()
  ggsave("elasticnet_final_coefficients.png", p_coef, width = 8, height = 5)
}

cat("\n分析完成！结果已保存。\n")

# 在确定 best_alpha 之后，可以绘制最佳 alpha 下的 CV 曲线
cv_best <- cv.glmnet(X, Y, alpha = best_alpha, nfolds = 10, standardize = TRUE)
plot(cv_best)
# 或者用 ggplot 美化
library(ggplot2)
cv_df <- data.frame(loglambda = log(cv_best$lambda), 
                    cvm = cv_best$cvm, 
                    cvup = cv_best$cvup, 
                    cvlo = cv_best$cvlo)
ggplot(cv_df, aes(x = loglambda, y = cvm)) +
  geom_ribbon(aes(ymin = cvlo, ymax = cvup), alpha = 0.2) +
  geom_line(color = "blue") +
  geom_vline(xintercept = log(cv_best$lambda.1se), linetype = "dashed") +
  geom_vline(xintercept = log(cv_best$lambda.min), linetype = "dotted") +
  labs(title = "交叉验证误差 vs log(λ)", x = "log(λ)", y = "MSE") +
  theme_minimal()
ggsave("cv_error_curve.png")

# 使用最佳 alpha 拟合完整路径
fit_path <- glmnet(X, Y, alpha = best_alpha, standardize = TRUE)
plot(fit_path, xvar = "lambda", label = FALSE)
# 或使用 ggplot（需提取系数）
coef_path <- as.matrix(fit_path$beta)
lambda_seq <- fit_path$lambda
# 这里可自行构造长格式数据绘图

# ========== Elastic Net 系数路径图（论文美化版） ==========

library(glmnet)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(viridisLite)
library(grid)

# 1. 使用最佳 alpha 拟合完整路径
fit_path <- glmnet(
  X, 
  Y, 
  alpha = best_alpha,
  standardize = TRUE
)


# 2. 提取系数矩阵
coef_path <- as.matrix(fit_path$beta)
lambda_seq <- fit_path$lambda


# 3. 选择贡献最大的基因
max_abs <- apply(abs(coef_path), 1, max)

top_genes <- names(
  sort(max_abs, decreasing = TRUE)
)[1:min(20, nrow(coef_path))]


coef_sub <- coef_path[top_genes, , drop = FALSE]


# 4. 转成长格式
df_long <- data.frame(
  log_lambda = rep(log(lambda_seq), each = length(top_genes)),
  Coefficient = as.vector(coef_sub),
  Gene = rep(top_genes, times = length(lambda_seq))
)



# 5. 获取每个基因最后一个lambda位置用于标签
label_df <- df_long %>%
  group_by(Gene) %>%
  filter(log_lambda == min(log_lambda)) %>%
  ungroup()



# 6. 绘图
p_path <- ggplot(
  df_long,
  aes(
    x = log_lambda,
    y = Coefficient,
    group = Gene,
    color = Gene
  )
) +
  
  geom_line(
    linewidth = 0.7,
    alpha = 0.85
  ) +
  
  # 零线
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey50"
  ) +
  
  # 末端基因标签
  geom_text_repel(
    data = label_df,
    aes(label = Gene),
    size = 3.2,
    direction = "y",
    hjust = 0,
    nudge_x = 0.15,
    segment.size = 0.2,
    show.legend = FALSE
  ) +
  
  
  scale_color_viridis_d(
    option = "D",
    begin = 0.1,
    end = 0.9
  ) +
  
  labs(
    title = "Elastic Net Coefficient Path",
    subtitle = paste0(
      "Optimal α = ",
      best_alpha,
      " | Top 20 genes ranked by maximum coefficient magnitude"
    ),
    x = expression(log(lambda)),
    y = "Coefficient"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    
    # 删除legend
    legend.position = "none",
    
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5,
      color = "grey40"
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    axis.text = element_text(
      color = "black"
    ),
    
    plot.margin = margin(
      10, 
      80, 
      10, 
      10
    )
  )



# 7. 输出PDF
pdf(
  "ElasticNet_coefficient_path_beautified.pdf",
  width = 10,
  height = 7
)

print(p_path)

dev.off()


cat(
  "Elastic Net 系数路径图已保存\n"
)

# ==================== 弹性网络稳定性详细报告 ====================
# 前提：已运行修正版弹性网络，存在：
#   selected_genes_list (list), X (样本×基因矩阵), Y (向量), best_alpha
#   corr_results (数据框，至少包含 Gene, Pearson_r, Pearson_p, Spearman_rho, Spearman_p)

library(glmnet)
library(ggplot2)
library(reshape2)

# ---------- 1. 获取所有进入弹性网络的基因（高方差基因）----------
genes_in_elastic <- colnames(X)
cat("弹性网络输入基因数：", length(genes_in_elastic), "\n")

# ---------- 2. 对每个基因，计算稳定性统计 ----------
n_repeats <- length(selected_genes_list)

# 存储系数（每次重复中，如果基因被选中则记录系数，否则 NA）
coef_matrix <- matrix(NA, nrow = length(genes_in_elastic), ncol = n_repeats,
                      dimnames = list(genes_in_elastic, paste0("rep", 1:n_repeats)))

# 重新运行每次重复，获取系数（需要重新拟合模型，因为之前只保存了基因名）
cat("正在重新拟合每次重复以提取系数（可能需要几分钟）...\n")
for (rep in 1:n_repeats) {
  set.seed(123 + rep * 10)  # 与之前相同的种子
  cv_rep <- cv.glmnet(X, Y, alpha = best_alpha, nfolds = 10, standardize = TRUE)
  model_rep <- glmnet(X, Y, alpha = best_alpha, lambda = cv_rep$lambda.1se, standardize = TRUE)
  coef_rep <- as.matrix(coef(model_rep))[-1, , drop = FALSE]  # 去掉截距，行名=基因
  # 填充系数矩阵
  coef_matrix[rownames(coef_rep), rep] <- coef_rep[, 1]
  if (rep %% 10 == 0) cat(sprintf("已完成 %d/%d 次重复\n", rep, n_repeats))
}

# 计算每个基因的统计量
stability_stats <- data.frame(
  Gene = genes_in_elastic,
  Selection_Frequency = rowMeans(!is.na(coef_matrix)),  # 非NA比例 = 被选中频率
  Mean_Coeff = rowMeans(coef_matrix, na.rm = TRUE),
  SD_Coeff = apply(coef_matrix, 1, sd, na.rm = TRUE),
  Median_Coeff = apply(coef_matrix, 1, median, na.rm = TRUE),
  Coeff_NonZero_Ratio = rowMeans(coef_matrix != 0 & !is.na(coef_matrix), na.rm = TRUE)
)
# 将 NaN（从未被选中的基因）替换为 0
stability_stats[is.na(stability_stats$Mean_Coeff), "Mean_Coeff"] <- 0
stability_stats[is.na(stability_stats$SD_Coeff), "SD_Coeff"] <- 0

# ---------- 3. 合并相关性结果 ----------
# 确保 corr_results 包含所有基因（包括未进入弹性网络的）
final_report <- merge(corr_results, stability_stats, by = "Gene", all.x = TRUE)
# 对于未进入弹性网络的基因，稳定性统计设为 NA 或 0
final_report$Selection_Frequency[is.na(final_report$Selection_Frequency)] <- 0
final_report$Mean_Coeff[is.na(final_report$Mean_Coeff)] <- 0
final_report$SD_Coeff[is.na(final_report$SD_Coeff)] <- 0

# 按 |Pearson_r| 降序排列
final_report <- final_report[order(-abs(final_report$Pearson_r)), ]

# ---------- 4. 模型性能评估（每次重复的预测 R² 和 MSE）----------
# 计算每次重复的预测性能（使用交叉验证的 CVM 作为估计，或重新计算）
rep_performance <- data.frame(Repetition = 1:n_repeats, R2 = NA, MSE = NA)

for (rep in 1:n_repeats) {
  set.seed(123 + rep * 10)
  cv_rep <- cv.glmnet(X, Y, alpha = best_alpha, nfolds = 10, standardize = TRUE)
  # 使用 lambda.1se 的模型预测
  pred <- predict(cv_rep, newx = X, s = "lambda.1se")
  # 计算 R² 和 MSE
  ss_res <- sum((Y - pred)^2)
  ss_tot <- sum((Y - mean(Y))^2)
  r2 <- 1 - ss_res / ss_tot
  mse <- mean((Y - pred)^2)
  rep_performance$R2[rep] <- r2
  rep_performance$MSE[rep] <- mse
}

performance_summary <- data.frame(
  Metric = c("R2", "MSE"),
  Mean = c(mean(rep_performance$R2), mean(rep_performance$MSE)),
  SD = c(sd(rep_performance$R2), sd(rep_performance$MSE)),
  Min = c(min(rep_performance$R2), min(rep_performance$MSE)),
  Max = c(max(rep_performance$R2), max(rep_performance$MSE))
)

# ---------- 5. 保存结果 ----------
write.csv(final_report, "gene_stability_full_report.csv", row.names = FALSE)
write.csv(rep_performance, "elasticnet_repeat_performance.csv", row.names = FALSE)
write.csv(performance_summary, "elasticnet_performance_summary.csv", row.names = FALSE)

# 筛选出稳定基因（频率 >= 0.8）并单独保存
stable_genes_report <- final_report[final_report$Selection_Frequency >= 0.8, ]
write.csv(stable_genes_report, "stable_genes_with_stats.csv", row.names = FALSE)

cat("\n稳定性报告已保存：\n",
    "- gene_stability_full_report.csv: 所有基因的相关性 + 弹性网络稳定性统计\n",
    "- stable_genes_with_stats.csv: 频率 >= 0.8 的稳定基因详细表格\n",
    "- elasticnet_repeat_performance.csv: 每次重复的 R² 和 MSE\n",
    "- elasticnet_performance_summary.csv: 性能汇总\n")

# ---------- 6. 可视化 ----------
# (1) 选择频率分布直方图
p1 <- ggplot(final_report, aes(x = Selection_Frequency)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  labs(title = "Distribution of Selection Frequency (50 repeats)",
       x = "Selection Frequency", y = "Number of Genes") +
  theme_minimal()
ggsave("selection_frequency_histogram.png", p1, width = 6, height = 4)

# (2) 频率 vs 相关性散点图（标记稳定基因）
p2 <- ggplot(final_report, aes(x = Pearson_r, y = Selection_Frequency)) +
  geom_point(alpha = 0.5, color = "gray40") +
  geom_point(data = stable_genes_report, aes(x = Pearson_r, y = Selection_Frequency),
             color = "red", size = 2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  labs(title = "Selection Frequency vs Pearson Correlation",
       x = "Pearson r", y = "Elastic Net Selection Frequency") +
  theme_minimal()
ggsave("frequency_vs_correlation.png", p2, width = 6, height = 5)

# (3) 系数箱线图（仅展示频率 > 0.5 的基因，最多 30 个）
coeff_data <- coef_matrix[, ]  # 基因 × 重复
coeff_melt <- melt(coeff_data, varnames = c("Gene", "Repetition"), value.name = "Coefficient")
coeff_melt <- coeff_melt[!is.na(coeff_melt$Coefficient), ]
# 计算每个基因的频率并筛选
gene_freq <- rowMeans(!is.na(coef_matrix))
freq_genes <- names(gene_freq[gene_freq > 0.5])
if (length(freq_genes) > 30) freq_genes <- freq_genes[order(gene_freq[freq_genes], decreasing = TRUE)[1:30]]
coeff_subset <- coeff_melt[coeff_melt$Gene %in% freq_genes, ]
if (nrow(coeff_subset) > 0) {
  p3 <- ggplot(coeff_subset, aes(x = reorder(Gene, Coefficient, FUN = median), y = Coefficient)) +
    geom_boxplot(fill = "lightblue", outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    coord_flip() +
    labs(title = "Elastic Net Coefficients (genes selected >50% repeats)",
         x = "Gene", y = "Coefficient") +
    theme_minimal()
  ggsave("coefficients_boxplot.png", p3, width = 8, height = 6)
}

# (4) 性能稳定性：R² 和 MSE 的折线图
rep_perf_long <- reshape2::melt(rep_performance, id.vars = "Repetition", 
                                variable.name = "Metric", value.name = "Value")
p4 <- ggplot(rep_perf_long, aes(x = Repetition, y = Value, color = Metric)) +
  geom_line() + geom_point(size = 1) +
  facet_wrap(~Metric, scales = "free_y", ncol = 1) +
  labs(title = "Model Performance Across 50 Repeats (lambda.1se)",
       x = "Repetition", y = "Value") +
  theme_minimal()
ggsave("performance_across_repeats.png", p4, width = 6, height = 6)

cat("\n所有图表已保存。\n")


# ==================== 1. 加载必要的包 ====================
# 若未安装，请先运行以下命令（去掉注释）：
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2"))

library(clusterProfiler)   # 核心富集分析包
library(org.Hs.eg.db)      # 人类基因注释数据库
library(enrichplot)        # 富集结果可视化
library(ggplot2)           # 绘图

# ==================== 2. 准备基因列表（假设已获得稳定基因） ====================
# 假设我们上一步的稳定基因列表存放在向量 'stable_genes' 中
# 这里提供一个示例，实际使用时请替换为您自己的基因名列表
# stable_genes <- c("GENE1", "GENE2", "GENE3", ...)

# 检查是否有基因
cat("稳定基因数量：", length(stable_genes), "\n")
if (length(stable_genes) == 0) stop("稳定基因列表为空，无法进行富集分析。")

# ==================== 3. 将基因符号转换为 Entrez ID ====================
# 注意：clusterProfiler 的 enrichGO 默认使用 Entrez ID，因此需进行转换
gene_entrez <- bitr(stable_genes, 
                    fromType = "SYMBOL",          # 输入类型：基因符号
                    toType = "ENTREZID",          # 输出类型：Entrez ID
                    OrgDb = org.Hs.eg.db)         # 使用人类注释库

# 提取 Entrez ID 向量
entrez_ids <- gene_entrez$ENTREZID
cat("成功转换为 Entrez ID 的基因数：", length(entrez_ids), "\n")

# 若有未成功转换的基因，可查看：
if (length(entrez_ids) < length(stable_genes)) {
  cat("警告：部分基因未能转换，请检查基因名是否正确。\n")
}

# ==================== 4. 执行 GOBP 富集分析 ====================
# 对 Biological Process (BP) 进行富集分析
go_bp <- enrichGO(
  gene = entrez_ids,                # 输入基因的 Entrez ID
  OrgDb = org.Hs.eg.db,            # 注释数据库
  ont = "BP",                       # 指定分析生物过程 (BP)
  pAdjustMethod = "BH",             # 多重检验校正方法 (Benjamini-Hochberg)
  pvalueCutoff = 0.05,              # p 值阈值
  qvalueCutoff = 0.05,              # q 值阈值（FDR 校正后）
  readable = TRUE                   # 将 Entrez ID 转换为基因符号，便于阅读
)

# ==================== 5. 查看富集结果 ====================
# 显示前 10 个显著富集的 GO 条目
cat("\n===== GOBP 富集分析结果（前 10 行）=====\n")
print(head(go_bp@result, 10))

# 检查是否有显著富集的条目
if (is.null(go_bp) || nrow(go_bp@result) == 0) {
  cat("\n警告：未找到显著富集的 GO 条目。可尝试调整 pvalueCutoff 或 qvalueCutoff。\n")
} else {
  cat("\n找到显著富集的 GO 条目数量：", nrow(go_bp@result), "\n")
}

# ==================== 6. 可视化富集结果 ====================
# 设置展示前 15 个最显著的条目
show_num <- min(15, nrow(go_bp@result))

# 6.1 气泡图（最常用）
p1 <- dotplot(go_bp, showCategory = show_num, title = "GOBP Enrichment - Top 15")
print(p1)

# 6.2 条形图
p2 <- barplot(go_bp, showCategory = show_num, title = "GOBP Enrichment - Top 15")
print(p2)

# 6.3 网络图（展示条目间的基因重叠关系）
# 注意：若显著条目过少，此图可能不适用
if (nrow(go_bp@result) >= 5) {
  p3 <- cnetplot(go_bp, showCategory = show_num, 
                 categorySize = "pvalue", 
                 foldChange = NULL,
                 title = "GOBP Enrichment Network")
  print(p3)
} else {
  cat("显著条目不足5个，跳过网络图绘制。\n")
}

# 6.4 富集图（enrichment map，显示条目间语义相似性）
if (nrow(go_bp@result) >= 5) {
  p4 <- emapplot(go_bp, showCategory = show_num, 
                 layout = "kk", 
                 title = "GOBP Enrichment Map")
  print(p4)
}

# ==================== 7. 保存结果 ====================
# 7.1 保存完整富集结果表
write.csv(go_bp@result, "GOBP_enrichment_results.csv", row.names = FALSE)

# 7.2 保存图形
ggsave("GOBP_dotplot.png", p1, width = 10, height = 8, dpi = 300)
ggsave("GOBP_barplot.png", p2, width = 10, height = 8, dpi = 300)
if (exists("p3")) ggsave("GOBP_cnetplot.png", p3, width = 12, height = 10, dpi = 300)
if (exists("p4")) ggsave("GOBP_emapplot.png", p4, width = 12, height = 10, dpi = 300)

# 7.3 保存会话信息（便于复现）
sink("GOBP_analysis_session_info.txt")
sessionInfo()
sink()

cat("\n分析完成！结果已保存至当前工作目录。\n")


# ==================== Residual-driven driver discovery（使用 GOBP 稳定基因集） ====================
# 目标：使用经过稳定性筛选的细胞周期核心基因集预测细胞周期评分，得到残差后在全基因组中寻找上游驱动因子
# 前提：已存在 expr_log (基因×样本)，cellcycle_df (Sample, CellCycle_score)
#       stable_genes_for_GOBP (用于 GOBP 的基因集，弹性网络频率≥0.8)

library(glmnet)
library(ggplot2)
library(ggpubr)

# ---------- 1. 准备数据 ----------
common_samples <- intersect(colnames(expr_log), cellcycle_df$Sample)
Y <- cellcycle_df$CellCycle_score[match(common_samples, cellcycle_df$Sample)]

# 使用 GOBP 用的稳定基因集（不是原始全部细胞周期基因）
# 确保这些基因在表达矩阵中
genes_for_prediction <- intersect(stable_genes, rownames(expr_log))
cat("用于预测的基因数（稳定核心细胞周期基因）：", length(genes_for_prediction), "\n")
if (length(genes_for_prediction) < 5) stop("基因数太少，无法进行弹性网络预测")

X_core <- t(expr_log[genes_for_prediction, common_samples])  # 样本 × 核心基因

# ---------- 2. 弹性网络预测（使用 lambda.1se 避免过拟合） ----------
set.seed(123)
cv_fit <- cv.glmnet(x = X_core, y = Y, alpha = 0.5, nfolds = 10, standardize = TRUE)
lambda_1se <- cv_fit$lambda.1se
final_model <- glmnet(x = X_core, y = Y, alpha = 0.5, lambda = lambda_1se, standardize = TRUE)
Y_pred <- predict(final_model, newx = X_core, s = lambda_1se)[, 1]
residuals <- Y - Y_pred
cat("模型 R² =", 1 - var(residuals)/var(Y), "\n")

# ---------- 3. 全基因组（排除用于预测的核心基因）与残差的相关性 ----------
all_genes <- rownames(expr_log)
non_core_genes <- setdiff(all_genes, genes_for_prediction)   # 排除核心基因集
cat("全基因组其他基因数：", length(non_core_genes), "\n")

# 循环计算相关性
gene_resid_corr <- data.frame(Gene = non_core_genes, 
                              Pearson_r = NA, 
                              Pearson_p = NA)
pb <- txtProgressBar(min = 0, max = length(non_core_genes), style = 3)
for (i in seq_along(non_core_genes)) {
  g <- non_core_genes[i]
  expr_g <- expr_log[g, common_samples]
  ct <- cor.test(expr_g, residuals, method = "pearson")
  gene_resid_corr$Pearson_r[i] <- ct$estimate
  gene_resid_corr$Pearson_p[i] <- ct$p.value
  setTxtProgressBar(pb, i)
}
close(pb)

gene_resid_corr$FDR <- p.adjust(gene_resid_corr$Pearson_p, method = "fdr")
gene_resid_corr <- gene_resid_corr[order(-abs(gene_resid_corr$Pearson_r)), ]

# ---------- 4. 检查 ADAM17 ----------
adam17_entry <- gene_resid_corr[gene_resid_corr$Gene %in% c("ADAM17", "Adam17", "TACE"), ]
if (nrow(adam17_entry) > 0) {
  cat("\n===== ADAM17 与残差的相关性 =====\n")
  print(adam17_entry)
  # 散点图
  adam17_expr <- expr_log[adam17_entry$Gene[1], common_samples]
  df_plot <- data.frame(Residual = residuals, ADAM17 = adam17_expr)
  p <- ggplot(df_plot, aes(x = Residual, y = ADAM17)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    labs(title = "ADAM17 vs Cell Cycle Score Residual (using core cycle genes)",
         x = "Residual (observed - predicted by core cycle genes)",
         y = "ADAM17 expression (log2(TPM+1))") +
    theme_minimal()
  print(p)
  ggsave("ADAM17_vs_residual_using_core_genes.png", p, width = 6, height = 5, dpi = 300)
} else {
  cat("\nADAM17 未出现在非核心基因列表中（或未找到）。\n")
  # 检查 ADAM17 是否在核心基因集里
  if ("ADAM17" %in% genes_for_prediction) cat("注意：ADAM17 属于核心细胞周期基因集，已从残差分析中排除。\n")
}

# ---------- 5. 筛选显著驱动基因（例如 |r|>0.3, FDR<0.05） ----------
r_thresh <- 0.3
fdr_thresh <- 0.05
drivers <- gene_resid_corr[abs(gene_resid_corr$Pearson_r) > r_thresh & 
                             gene_resid_corr$FDR < fdr_thresh, ]
cat(sprintf("\n显著驱动基因数（|r|>%.1f, FDR<%.2f）：%d\n", r_thresh, fdr_thresh, nrow(drivers)))

# ---------- 6. （可选）保存结果 ----------
write.csv(gene_resid_corr, "residual_correlation_all_non_core_genes.csv", row.names = FALSE)
write.csv(drivers, "residual_drivers_significant.csv", row.names = FALSE)

# 若需要，保存残差值供后续分析
residual_df <- data.frame(Sample = common_samples, Residual = residuals)
write.csv(residual_df, "cellcycle_score_residuals.csv", row.names = FALSE)

cat("\n分析完成！结果已保存。\n")