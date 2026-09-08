# ============================================================================
# 完整流程：基于CGGA数据的JNK/ADAM17分组药敏预测分析 (oncoPredict策略)
# ============================================================================
# 作者：Bioinformatics R Script
# 日期：自动生成
# 说明：本脚本包含数据读取、分组、GDSC2药物敏感性预测、两组差异分析及可视化
# ============================================================================

# ----------------------------- 第1节：安装与加载包 -----------------------------
# 如果以下包未安装，请取消对应行的注释并运行一次
install.packages("oncoPredict")      # 主包
# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("ggpubr")
# install.packages("pheatmap")
# install.packages("RColorBrewer")

library(oncoPredict)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)

set.seed(12345)   # 保证可重复性

# ----------------------------- 第2节：设置路径（请修改为您的实际路径） -----------------------------
data_path <- "E:/ADAM17 data/CGGA.mRNAseq_693.RSEM-genes.20200506.txt"   # CGGA表达数据
train_dir <- "E:/ADAM17 data/Training Data/"          # 存放GDSC2训练数据的文件夹

# 训练数据文件名（请确保这两个文件存在于 train_dir 中）
train_expr_file <- file.path(train_dir, "GDSC2_Expr (RMA Normalized and Log Transformed).rds")
train_res_file  <- file.path(train_dir, "GDSC2_Res.rds")

# ----------------------------- 第3节：读取CGGA数据 -----------------------------
cat("===== 读取CGGA数据 =====\n")
if (!file.exists(data_path)) {
  stop(paste("错误：找不到CGGA数据文件:", data_path))
}
cgga_data <- read.table(data_path,
                        header = TRUE,
                        row.names = 1,
                        sep = "\t",
                        check.names = FALSE)
cat(paste0("CGGA数据维度: ", nrow(cgga_data), " 个基因, ", ncol(cgga_data), " 个样本\n"))

# ----------------------------- 第4节：确定JNK和ADAM17基因名 -----------------------------
cat("\n===== 查找JNK和ADAM17基因名 =====\n")
cat("数据中与JNK相关的基因：\n")
print(grep("MAPK8|MAPK9|MAPK10|JNK", rownames(cgga_data), value = TRUE, ignore.case = TRUE))
cat("\n数据中与ADAM17相关的基因：\n")
print(grep("ADAM17", rownames(cgga_data), value = TRUE))

# ===== 请根据上面输出，手动修改下面两个变量 =====
jnk_gene <- "MAPK8"      # 若输出为"JNK"则改为"JNK"
adam_gene <- "ADAM17"    # 若输出不同则相应修改

# 检查基因是否存在于数据中
if (!jnk_gene %in% rownames(cgga_data)) {
  stop(paste("错误：基因", jnk_gene, "不在数据中，请重新设置jnk_gene"))
}
if (!adam_gene %in% rownames(cgga_data)) {
  stop(paste("错误：基因", adam_gene, "不在数据中，请重新设置adam_gene"))
}

# ----------------------------- 第5节：四分组并提取目标两组 -----------------------------
cat("\n===== 进行四分组 =====\n")
jnk_expr <- as.numeric(cgga_data[jnk_gene, ])
adam_expr <- as.numeric(cgga_data[adam_gene, ])

# 按中位数划分高/低
jnk_high <- jnk_expr > median(jnk_expr, na.rm = TRUE)
adam_high <- adam_expr > median(adam_expr, na.rm = TRUE)

group_4 <- case_when(
  jnk_high & adam_high ~ "High_High",
  jnk_high & !adam_high ~ "High_Low",
  !jnk_high & adam_high ~ "Low_High",
  !jnk_high & !adam_high ~ "Low_Low"
)
names(group_4) <- colnames(cgga_data)
cat("四分组样本数量：\n")
print(table(group_4))

# 筛选目标两组
target_groups <- c("High_Low", "Low_High")
samples_interest <- names(group_4)[group_4 %in% target_groups]
cat(paste0("\n目标两组样本数: ", length(samples_interest), "\n"))

if (length(samples_interest) < 6) {
  warning("目标两组总样本数少于6，统计检验可能不可靠")
}

expr_for_pred <- cgga_data[, samples_interest, drop = FALSE]
group_factor <- factor(group_4[samples_interest], levels = c("High_Low", "Low_High"))

# ----------------------------- 第6节：加载GDSC2训练数据 -----------------------------
cat("\n===== 加载GDSC2训练数据 =====\n")
if (!file.exists(train_expr_file)) {
  stop(paste("错误：找不到训练表达文件", train_expr_file, 
             "\n请从 https://osf.io/c6tfx/ 下载 GDSC2_Expr (RMA Normalized and Log Transformed).rds 并放入", train_dir))
}
if (!file.exists(train_res_file)) {
  stop(paste("错误：找不到训练响应文件", train_res_file,
             "\n请从 https://osf.io/c6tfx/ 下载 GDSC2_Res.rds 并放入", train_dir))
}

trainingExprData <- readRDS(train_expr_file)
trainingPtype <- readRDS(train_res_file)

# GDSC2的IC50存储为log值，需转回原始值（oncoPredict默认期望原始IC50）
trainingPtype <- exp(trainingPtype)

cat(paste0("训练表达矩阵: ", nrow(trainingExprData), " 基因, ", ncol(trainingExprData), " 细胞系\n"))
cat(paste0("训练响应矩阵: ", nrow(trainingPtype), " 细胞系, ", ncol(trainingPtype), " 药物\n"))

# ----------------------------- 第7节：药物敏感性预测 -----------------------------
cat("\n===== 开始药物敏感性预测 (此步骤可能需要较长时间) =====\n")
testExprData <- as.matrix(expr_for_pred)

predicted_ic50 <- calcPhenotype(
  trainingExprData = trainingExprData,
  trainingPtype = trainingPtype,
  testExprData = testExprData,
  batchCorrect = "eb",
  powerTransformPhenotype = TRUE,
  removeLowVaryingGenes = 0.2,
  minNumSamples = 10,
  printOutput = TRUE,
  removeLowVaringGenesFrom = 'homogenizeData'  # <--- 添加这一行
)

cat("\n预测完成！\n")
cat(paste0("预测结果维度: ", nrow(predicted_ic50), " 个样本, ", ncol(predicted_ic50), " 种药物\n"))

# ----------------------------- 第8节：整理预测结果并保存 -----------------------------
ic50_df <- as.data.frame(predicted_ic50)
ic50_df$Sample <- rownames(ic50_df)
ic50_df$Group <- group_factor[rownames(ic50_df)]

# 保存完整结果
write.csv(ic50_df, file = "Predicted_IC50_all_drugs.csv", row.names = FALSE)
cat("\n完整预测结果已保存至: Predicted_IC50_all_drugs.csv\n")

# ----------------------------- 第9节：两组间药物敏感性差异分析 -----------------------------
cat("\n===== 进行两组间差异分析 =====\n")
drug_names <- setdiff(colnames(ic50_df), c("Sample", "Group"))

diff_results <- data.frame(
  Drug = character(),
  P_value = numeric(),
  Log2FC = numeric(),
  stringsAsFactors = FALSE
)

for (drug in drug_names) {
  high_low_ic50 <- ic50_df[ic50_df$Group == "High_Low", drug]
  low_high_ic50 <- ic50_df[ic50_df$Group == "Low_High", drug]
  
  if (length(high_low_ic50) >= 3 && length(low_high_ic50) >= 3) {
    test_res <- t.test(low_high_ic50, high_low_ic50)
    log2fc <- log2(mean(low_high_ic50, na.rm = TRUE) / mean(high_low_ic50, na.rm = TRUE))
    diff_results <- rbind(diff_results, data.frame(
      Drug = drug,
      P_value = test_res$p.value,
      Log2FC = log2fc
    ))
  }
}

diff_results$Adj_P_value <- p.adjust(diff_results$P_value, method = "BH")
diff_results$Significant <- ifelse(diff_results$Adj_P_value < 0.05 & abs(diff_results$Log2FC) > 0.5,
                                   "Yes", "No")
diff_results <- diff_results[order(diff_results$P_value), ]

cat(paste0("显著差异药物数量 (adj.P<0.05 且 |log2FC|>0.5): ", sum(diff_results$Significant == "Yes"), "\n"))
write.csv(diff_results, file = "Drug_sensitivity_diff_results.csv", row.names = FALSE)

# ----------------------------- 第10节：可视化 ---------------------------------
cat("\n===== 生成可视化图表 =====\n")

# 10.1 火山图
p_volcano <- ggplot(diff_results, aes(x = Log2FC, y = -log10(Adj_P_value), color = Significant)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("No" = "grey70", "Yes" = "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(title = "Drug Sensitivity: Low_High vs High_Low",
       x = "Log2 Fold Change (Low_High / High_Low)",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("Volcano_Plot_Drug_Sensitivity.png", p_volcano, width = 8, height = 6, dpi = 300)

# 10.2 Top显著药物箱线图
sig_drugs <- diff_results$Drug[diff_results$Significant == "Yes"]
if (length(sig_drugs) > 0) {
  top_n <- min(6, length(sig_drugs))
  top_drugs <- sig_drugs[1:top_n]
  plot_data <- ic50_df %>%
    select(Sample, Group, all_of(top_drugs)) %>%
    pivot_longer(cols = -c(Sample, Group), names_to = "Drug", values_to = "IC50")
  
  p_box <- ggplot(plot_data, aes(x = Group, y = IC50, fill = Group)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.8) +
    facet_wrap(~ Drug, scales = "free_y", ncol = 3) +
    stat_compare_means(method = "t.test", label = "p.format",
                       comparisons = list(c("High_Low", "Low_High"))) +
    labs(title = "Top Significant Drugs: IC50 Comparison") +
    theme_minimal() +
    theme(legend.position = "bottom")
  ggsave("Top_Drugs_Boxplot.png", p_box, width = 12, height = 8, dpi = 300)
} else {
  cat("没有显著差异药物，跳过箱线图生成。\n")
}

# 10.3 热图（显著差异药物）
if (length(sig_drugs) > 0) {
  heatmap_data <- ic50_df %>%
    select(all_of(sig_drugs)) %>%
    as.matrix()
  rownames(heatmap_data) <- ic50_df$Sample
  
  ann_col <- data.frame(Group = ic50_df$Group)
  rownames(ann_col) <- ic50_df$Sample
  
  pheatmap(heatmap_data,
           scale = "row",
           annotation_row = ann_col,
           main = paste0("Drug Sensitivity Heatmap (", length(sig_drugs), " significant drugs)"),
           fontsize = 8,
           filename = "Drug_Sensitivity_Heatmap.png",
           width = 10, height = 8)
} else {
  cat("没有显著差异药物，跳过热图生成。\n")
}

# ----------------------------- 第11节：特定药物分析（替莫唑胺示例） -----------------------------
cat("\n===== 特定药物分析（替莫唑胺） =====\n")
tmz_related <- grep("Temozolomide|TMZ", drug_names, value = TRUE, ignore.case = TRUE)
if (length(tmz_related) > 0) {
  tmz_drug <- tmz_related[1]
  p_tmz <- ggplot(ic50_df, aes(x = Group, y = .data[[tmz_drug]], fill = Group)) +
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = paste("Predicted IC50 for", tmz_drug),
         x = "Group", y = "Predicted IC50") +
    theme_minimal() +
    theme(legend.position = "bottom")
  ggsave("TMZ_sensitivity_comparison.png", p_tmz, width = 6, height = 5, dpi = 300)
  cat("替莫唑胺分析图已保存: TMZ_sensitivity_comparison.png\n")
} else {
  cat("未在药物列表中找到替莫唑胺相关药物名。\n")
}

# ----------------------------- 第12节：完成提示 -----------------------------
cat("\n===== 所有分析完成！生成的文件如下： =====\n")
cat("  1. Predicted_IC50_all_drugs.csv\n")
cat("  2. Drug_sensitivity_diff_results.csv\n")
cat("  3. Volcano_Plot_Drug_Sensitivity.png\n")
cat("  4. Top_Drugs_Boxplot.png (若有显著药物)\n")
cat("  5. Drug_Sensitivity_Heatmap.png (若有显著药物)\n")
cat("  6. TMZ_sensitivity_comparison.png (若存在替莫唑胺)\n")
cat("\n请检查输出文件并解读结果。\n")
sessionInfo()

# ============================================================================
# 重写可视化部分：全部输出为 PDF + 生成“所有趋势药物合并大图”
# ============================================================================

library(ggplot2)
library(tidyr)
library(ggpubr)
library(patchwork)   # 用于拼图（如需组合趋势图）

cat("\n===== 开始生成 PDF 图表 =====\n")

# ------------------- 1. 确保数据已整理且分组已重命名 -------------------
if (!exists("ic50_df")) {
  stop("错误：ic50_df 不存在，请先成功运行 calcPhenotype() 并构建 ic50_df")
}

# 确认分组名称是否正确
cat("当前分组名称：", paste(levels(ic50_df$Group), collapse = ", "), "\n")

# ------------------- 2. 差异统计（若之前没跑过则重跑） -------------------
drug_names <- setdiff(colnames(ic50_df), c("Sample", "Group"))

diff_results <- data.frame(
  Drug = character(),
  P_value = numeric(),
  Log2FC = numeric(),
  stringsAsFactors = FALSE
)

for (drug in drug_names) {
  g1 <- ic50_df[ic50_df$Group == "ADAM17_low_JNK_high", drug]
  g2 <- ic50_df[ic50_df$Group == "ADAM17_high_JNK_low", drug]
  if (length(g1) >= 3 && length(g2) >= 3) {
    test_res <- t.test(g2, g1)
    log2fc <- log2(mean(g2, na.rm = TRUE) / mean(g1, na.rm = TRUE))
    diff_results <- rbind(diff_results, data.frame(
      Drug = drug,
      P_value = test_res$p.value,
      Log2FC = log2fc
    ))
  }
}
diff_results$Adj_P_value <- p.adjust(diff_results$P_value, method = "BH")
diff_results$Significant <- ifelse(diff_results$Adj_P_value < 0.05 & abs(diff_results$Log2FC) > 0.5,
                                   "Yes", "No")
diff_results <- diff_results[order(diff_results$P_value), ]

# 保存 Excel 表格
write.csv(diff_results, file = "All_Drugs_Sensitivity_Results.csv", row.names = FALSE)
cat("✅ 1. Excel 表格已生成: All_Drugs_Sensitivity_Results.csv\n")

# ------------------- 3. 图1：Top 10 药物箱线图 (PDF) -------------------
cat("\n生成 Top 10 药物箱线图...\n")
top10_drugs <- head(diff_results$Drug, 10)

plot_top10 <- ic50_df %>%
  select(Sample, Group, all_of(top10_drugs)) %>%
  pivot_longer(cols = -c(Sample, Group), names_to = "Drug", values_to = "IC50")

plot_top10$Drug <- factor(plot_top10$Drug, levels = top10_drugs)

p_top10 <- ggplot(plot_top10, aes(x = Group, y = IC50, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.size = 0.8) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
  facet_wrap(~ Drug, scales = "free_y", ncol = 5) +
  stat_compare_means(method = "t.test", label = "p.format", size = 3,
                     comparisons = list(c("ADAM17_low_JNK_high", "ADAM17_high_JNK_low"))) +
  labs(title = "Top 10 Differential Drugs (by P-value)",
       x = NULL, y = "Predicted IC50") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(face = "bold"))

ggsave("Top10_Drugs_Boxplot.pdf", p_top10, width = 14, height = 10)
cat("✅ 2. Top10 图已保存: Top10_Drugs_Boxplot.pdf\n")

# ------------------- 4. 图2：【您要的总图】所有有趋势的药物合并大图 (PDF) -------------------
# 筛选“有趋势”的药物：adj.P < 0.05 且 |log2FC| > 0.5
trend_drugs <- diff_results$Drug[diff_results$Significant == "Yes"]

# 如果没有显著药物，则自动取 Top 20 作为“趋势药物”
if (length(trend_drugs) == 0) {
  cat("⚠️ 没有发现 adj.P<0.05 且 |log2FC|>0.5 的显著药物，自动改用 Top 20 药物绘制大图。\n")
  trend_drugs <- head(diff_results$Drug, 20)
} else {
  cat(paste0("共筛选出 ", length(trend_drugs), " 种有趋势的药物，准备绘制合并大图。\n"))
}

# 准备绘图数据
plot_all <- ic50_df %>%
  select(Sample, Group, all_of(trend_drugs)) %>%
  pivot_longer(cols = -c(Sample, Group), names_to = "Drug", values_to = "IC50")

# 按显著性排序因子水平（P值小的在前，方便阅读）
plot_all$Drug <- factor(plot_all$Drug, levels = trend_drugs)

# 动态计算分面列数（若药物数量 <= 6 则 2 列，<= 12 则 3 列，<= 30 则 4 列，否则 5 列）
n_drugs <- length(trend_drugs)
if (n_drugs <= 6) {
  ncol_facet <- 2
} else if (n_drugs <= 12) {
  ncol_facet <- 3
} else if (n_drugs <= 30) {
  ncol_facet <- 4
} else {
  ncol_facet <- 5
}
# 动态计算 PDF 高度（每个小图约 3.5 英寸高，让箱线图清晰可见）
plot_height <- ceiling(n_drugs / ncol_facet) * 3.5 + 1

# 绘制大图（所有趋势药物）
p_all <- ggplot(plot_all, aes(x = Group, y = IC50, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.size = 0.6) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
  facet_wrap(~ Drug, scales = "free_y", ncol = ncol_facet) +
  stat_compare_means(method = "t.test", label = "p.format", size = 2.5,
                     comparisons = list(c("ADAM17_low_JNK_high", "ADAM17_high_JNK_low"))) +
  labs(title = paste0("All Trend Drugs (n = ", n_drugs, ", adj.P < 0.05 & |log2FC| > 0.5)"),
       x = NULL, y = "Predicted IC50") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        strip.text = element_text(face = "bold", size = 9),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("All_Trend_Drugs_Boxplot.pdf", p_all, width = 16, height = plot_height, limitsize = FALSE)
cat(paste0("✅ 3. 所有趋势药物合并大图已保存: All_Trend_Drugs_Boxplot.pdf (", n_drugs, " 种药物, 高度 ", round(plot_height, 1), " 英寸)\n"))

# ------------------- 5. 图3（补充）：两张综合趋势箱线图 (PDF) -------------------
cat("\n生成综合趋势图（平均 Z-score + 显著药物占比）...\n")

# 计算每个样本的平均 Z-score（反映整体敏感性）
ic50_matrix <- as.matrix(ic50_df[, drug_names])
ic50_zscore <- scale(ic50_matrix)
ic50_df$Mean_Zscore <- rowMeans(ic50_zscore, na.rm = TRUE)

# 计算每个样本中“显著敏感/耐药”药物的数量（用于评估个体趋势）
# 这里我们统计每个样本中 IC50 低于中位数的药物占比（简化版：对每个药物以整体中位数为界计算敏感度评分）
# 更稳健：计算每个样本相对于所有样本的耐药指数
sample_median <- apply(ic50_matrix, 1, median, na.rm = TRUE)
ic50_df$Median_IC50 <- sample_median

# 绘制两个补充子图：平均 Z-score 和 中位数 IC50
p_trend1 <- ggplot(ic50_df, aes(x = Group, y = Mean_Zscore, fill = Group)) +
  geom_boxplot(width = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.5) +
  stat_compare_means(method = "t.test", label = "p.format", size = 4,
                     comparisons = list(c("ADAM17_low_JNK_high", "ADAM17_high_JNK_low"))) +
  labs(title = "Overall Sensitivity (Mean Z-score across all drugs)",
       x = NULL, y = "Mean Standardized IC50 (lower = more sensitive)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

p_trend2 <- ggplot(ic50_df, aes(x = Group, y = Median_IC50, fill = Group)) +
  geom_boxplot(width = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.5) +
  stat_compare_means(method = "t.test", label = "p.format", size = 4,
                     comparisons = list(c("ADAM17_low_JNK_high", "ADAM17_high_JNK_low"))) +
  labs(title = "Overall Sensitivity (Median IC50 across all drugs)",
       x = NULL, y = "Median Predicted IC50 (lower = more sensitive)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# 使用 patchwork 拼成一张图
p_combined <- (p_trend1 | p_trend2) + 
  plot_annotation(title = "Overall Drug Sensitivity Trend Comparison",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

ggsave("Overall_Trend_Summary.pdf", p_combined, width = 12, height = 6)
cat("✅ 4. 综合趋势总结图已保存: Overall_Trend_Summary.pdf\n")

# ------------------- 6. 完成信息 -------------------
cat("\n===== 🎉 所有 PDF 图表生成完毕 =====\n")
cat("输出文件清单：\n")
cat("  📊 1. All_Drugs_Sensitivity_Results.csv       （所有药物差异统计表格）\n")
cat("  📈 2. Top10_Drugs_Boxplot.pdf                 （Top 10 分面图）\n")
cat("  📈 3. All_Trend_Drugs_Boxplot.pdf             （【您要的总图】所有趋势药物合并大图）\n")
cat("  📈 4. Overall_Trend_Summary.pdf               （综合趋势补充图：平均Z-score + 中位数IC50）\n")
cat("\n提示：若 PDF 中文字或图形过小，可用 Adobe Acrobat 或浏览器放大查看。\n")

# ============================================================================
# Drug sensitivity volcano plot with labels
# ============================================================================


library(ggplot2)
library(ggrepel)


# ---------------------------
# 数据准备
# ---------------------------

# diff_results需要包含:
# Drug
# Log2FC
# FDR


diff_results$negLogFDR <-
  -log10(diff_results$FDR)



# ---------------------------
# 定义显著性
# ---------------------------

diff_results$Significance <- "NS"


diff_results$Significance[
  diff_results$FDR < 0.05 &
    diff_results$Log2FC > 0.5
] <- "Resistance in ADAM17-high"


diff_results$Significance[
  diff_results$FDR < 0.05 &
    diff_results$Log2FC < -0.5
] <- "Sensitivity in ADAM17-high"



# ---------------------------
# 选择需要标注的药物
# ---------------------------

# 综合评分:
# 越靠上 + 越偏离中心越优先


diff_results$Label_score <-
  abs(diff_results$Log2FC) *
  diff_results$negLogFDR



# Top 10
label_drugs <-
  diff_results %>%
  arrange(desc(Label_score)) %>%
  head(10)



# ---------------------------
# 绘图
# ---------------------------


p_volcano <- ggplot(
  diff_results,
  aes(
    x = Log2FC,
    y = negLogFDR,
    color = Significance
  )
)+
  
  
  geom_point(
    size=2.2,
    alpha=0.8
  )+
  
  
  
  # 阈值线
  
  geom_vline(
    xintercept=c(-0.5,0.5),
    linetype="dashed",
    color="grey50"
  )+
  
  
  geom_hline(
    yintercept=-log10(0.05),
    linetype="dashed",
    color="grey50"
  )+
  
  
  
  # 药物名称
  
  geom_text_repel(
    data = label_drugs,
    aes(
      label = Drug
    ),
    size=3.5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color="grey50"
  )+
  
  
  
  scale_color_manual(
    values=c(
      "NS"="grey75",
      "Resistance in ADAM17-high"="#E64B35",
      "Sensitivity in ADAM17-high"="#4DBBD5"
    )
  )+
  
  
  labs(
    title="Differential Drug Sensitivity",
    x="log2(IC50 Fold Change)",
    y="-log10(FDR)",
    color=NULL
  )+
  
  
  theme_classic(base_size=14)+
  
  theme(
    axis.text=element_text(color="black"),
    axis.title=element_text(face="bold"),
    plot.title=
      element_text(
        face="bold",
        hjust=0.5
      ),
    legend.position="top"
  )



# 输出PDF

ggsave(
  "Drug_sensitivity_volcano_labeled.pdf",
  p_volcano,
  width=8,
  height=7
)


