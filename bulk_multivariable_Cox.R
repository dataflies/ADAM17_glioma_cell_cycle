# ================================================================
# 多因素 Cox 回归分析及森林图（修正版）
# 目标人群：ADAM17 高表达（>中位数）且 JNK（MAPK8）低表达（<中位数）
# ================================================================

# 1. 加载必要包 ----------------------------------------------------
library(survival)
library(survminer)
library(dplyr)
library(tidyr)

# 2. 设置文件路径 --------------------------------------------------
clinical_file <- "E:/ADAM17 data/CGGA.mRNAseq_693_clinical.20200506.txt"
expr_file     <- "E:/ADAM17 data/CGGA.mRNAseq_693.RSEM-genes.20200506.txt"

# 3. 读取数据 ------------------------------------------------------
clinical <- read.delim(clinical_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
expr     <- read.delim(expr_file,     header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# 查看临床数据列名（调试用）
cat("临床数据列名：\n"); print(colnames(clinical))

# 4. 提取目标基因（ADAM17 和 MAPK8 / JNK）------------------------
gene1 <- "ADAM17"
gene2 <- "MAPK8"   # JNK 基因的官方符号为 MAPK8，若数据中为 "JNK" 请改回

if (!(gene1 %in% expr$Gene_Name)) stop(paste("基因", gene1, "不存在！"))
if (!(gene2 %in% expr$Gene_Name)) stop(paste("基因", gene2, "不存在！"))

expr_genes <- expr %>%
  filter(Gene_Name %in% c(gene1, gene2)) %>%
  column_to_rownames(var = "Gene_Name") %>%
  t() %>% as.data.frame() %>%
  rename(ADAM17 = !!gene1, JNK = !!gene2)

expr_genes$CGGA_ID <- rownames(expr_genes)

# 5. 计算中位数并筛选 ----------------------------------------------
median_ADAM17 <- median(expr_genes$ADAM17, na.rm = TRUE)
median_JNK    <- median(expr_genes$JNK,    na.rm = TRUE)
cat("ADAM17 中位数 =", median_ADAM17, "\n")
cat("JNK 中位数    =", median_JNK, "\n")

selected_samples <- expr_genes %>%
  filter(ADAM17 > median_ADAM17 & JNK < median_JNK)
cat("符合筛选条件的样本数：", nrow(selected_samples), "\n")

# 6. 合并临床数据 --------------------------------------------------
merged <- selected_samples %>%
  left_join(clinical, by = "CGGA_ID")

if (nrow(merged) == 0) stop("合并后无数据，请检查CGGA_ID匹配！")

# 7. 预处理临床变量（列名含有空格和特殊字符，使用反引号）----------
# 将需要使用的列名存入变量，便于引用
col_OS      <- "OS"
col_censor  <- "Censor (alive=0; dead=1)"
col_radio   <- "Radio_status (treated=1;un-treated=0)"
col_chemo   <- "Chemo_status (TMZ treated=1;un-treated=0)"
col_grade   <- "Grade"
col_age     <- "Age"
col_gender  <- "Gender"
col_idh     <- "IDH_mutation_status"
col_1p19q   <- "1p19q_codeletion_status"
col_mgmt    <- "MGMTp_methylation_status"

# 检查这些列是否存在于合并数据中
required_cols <- c(col_OS, col_censor, col_radio, col_chemo, col_grade,
                   col_age, col_gender, col_idh, col_1p19q, col_mgmt)
missing_cols <- required_cols[!required_cols %in% colnames(merged)]
if (length(missing_cols) > 0) {
  warning("以下列在合并数据中不存在，请检查临床数据：", 
          paste(missing_cols, collapse = ", "))
  # 可创建缺省列或退出
}

# 重命名生存相关列为简洁名称（避免后续引用困难）
merged <- merged %>%
  rename(
    OS      = !!col_OS,
    Censor  = !!col_censor,
    Radio   = !!col_radio,
    Chemo   = !!col_chemo
  )

# 因子转换（动态获取实际水平）
merged$Grade <- factor(merged$Grade)
merged$Gender <- factor(merged$Gender)
merged$IDH_mutation_status <- factor(merged$IDH_mutation_status)
# 处理 1p19q 和 MGMTp（可能含有缺失值）
merged$`1p19q_codeletion_status` <- factor(merged[[col_1p19q]])  # 注意列名含数字开头，用反引号或[[
merged$MGMTp_methylation_status <- factor(merged[[col_mgmt]])

# 放化疗转为因子（0/1 转 No/Yes）
merged$Radio <- factor(merged$Radio, levels = c(0,1), labels = c("No", "Yes"))
merged$Chemo <- factor(merged$Chemo, levels = c(0,1), labels = c("No", "Yes"))

# 年龄作为连续变量（可考虑分类，此处保持数值）
merged$Age <- as.numeric(merged$Age)

# 检查是否有缺失值（关键变量）
cat("\n生存变量及协变量缺失情况：\n")
summary(merged[, c("OS", "Censor", "Radio", "Chemo", "Grade", "Age", "Gender",
                   "IDH_mutation_status", "1p19q_codeletion_status", "MGMTp_methylation_status")])

# 8. 多因素 Cox 回归 ------------------------------------------------
# 使用更新后的列名
cox_model <- coxph(Surv(OS, Censor) ~ Radio + Chemo + Grade + Age + Gender +
                     IDH_mutation_status + `1p19q_codeletion_status` + MGMTp_methylation_status,
                   data = merged)

# 9. 输出结果 ------------------------------------------------------
summary(cox_model)

# 提取 HR、CI、P 值
res <- summary(cox_model)
hr <- round(res$coefficients[, "exp(coef)"], 3)
ci_lower <- round(res$conf.int[, "lower .95"], 3)
ci_upper <- round(res$conf.int[, "upper .95"], 3)
p <- round(res$coefficients[, "Pr(>|z|)"], 4)
results <- data.frame(
  Variable = rownames(res$coefficients),
  HR = hr,
  CI_lower = ci_lower,
  CI_upper = ci_upper,
  P = p
)
print(results)

# 10. 森林图 -------------------------------------------------------
forest_plot <- ggforest(cox_model, data = merged,
                        main = "多因素 Cox 回归森林图 (ADAM17高 & JNK低亚组)",
                        cpositions = c(0.10, 0.22, 0.40),
                        fontsize = 0.7,
                        refLabel = "reference",
                        noDigits = 3)

# 11. 保存结果 ----------------------------------------------------
# 保存模型摘要
sink("Cox_model_summary.txt")
print(summary(cox_model))
sink()

# 保存结果表格
write.csv(results, "Cox_results.csv", row.names = FALSE)

# 保存森林图
ggsave("Forest_plot.pdf", forest_plot, width = 8, height = 6)

cat("\n分析完成！结果文件已保存：Cox_model_summary.txt, Cox_results.csv, Forest_plot.pdf\n")
