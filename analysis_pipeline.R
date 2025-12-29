# 1. 环境准备 (请确保安装了 cowplot)
# ==========================================
# install.packages("cowplot") 
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(patchwork)
library(cowplot) # 引入 cowplot 处理网格对齐

# 读取与清洗数据部分保持不变
raw_df <- read.delim("E:/GSE102485gene.txt/GSE102485gene.txt", header = TRUE, check.names = FALSE)
rownames(raw_df) <- raw_df[,1]
raw_df <- raw_df[,-1]

target_ids <- c("EGFR" = "ENSG00000146648", "HERC4" = "ENSG00000101446", 
                "CTGF" = "ENSG00000118523", "ACTA2" = "ENSG00000107796")

clean_rows <- gsub("\\..*", "", rownames(raw_df))
plot_matrix <- raw_df[clean_rows %in% target_ids, ]
id_map <- setNames(names(target_ids), target_ids)
rownames(plot_matrix) <- id_map[clean_rows[clean_rows %in% target_ids]]

df_clean <- as.data.frame(t(plot_matrix)) %>%
  mutate(across(everything(), ~as.numeric(as.character(.)))) %>%
  na.omit()

df_clean$SampleID <- rownames(df_clean)
df_clean$Group <- factor(ifelse(1:nrow(df_clean) <= (ncol(raw_df)/2), "Control", "High_Glucose"), 
                         levels = c("Control", "High_Glucose"))

df_long <- df_clean %>%
  pivot_longer(cols = -c(Group, SampleID), names_to = "Gene", values_to = "Exp") %>%
  mutate(Exp_log = log2(Exp + 1))

# ==========================================
# 4. 绘图模块 A：差异表达箱线图 (p1)
# ==========================================
p1 <- ggplot(df_long, aes(x = Group, y = Exp_log, fill = Group)) +
  geom_violin(alpha = 0.2, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4, aes(color = Group)) +
  facet_wrap(~Gene, scales = "free_y", nrow = 1) +
  stat_compare_means(method = "t.test", label = "p.signif") + 
  scale_fill_manual(values = c("Control" = "#74add1", "High_Glucose" = "#f46d43")) +
  scale_color_manual(values = c("Control" = "#4575b4", "High_Glucose" = "#d73027")) +
  theme_bw() +
  labs(y = "log2(FPKM + 1)", x = "", title = "A. Differential Expression Analysis") +
  theme(strip.text = element_text(face = "italic"), legend.position = "none")

# ==========================================
# 5. 绘图模块 B：轴相关性散点图 (p2)
# ==========================================
p2 <- ggplot(df_clean, aes(x = log2(EGFR + 1), y = log2(HERC4 + 1))) +
  geom_point(aes(color = Group), size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  stat_cor(method = "pearson", digits = 3) +
  scale_color_manual(values = c("Control" = "#74add1", "High_Glucose" = "#f46d43")) +
  theme_classic() +
  labs(x = "EGFR (log2)", y = "HERC4 (log2)", title = "B. Axis Correlation")

# ==========================================
# 6. 绘图模块 C：相关性热图 (解决覆盖的关键)
# ==========================================
# 步骤1: 创建一个绘制热图的函数
draw_corr <- function() {
  cor_res <- cor(df_clean[, names(target_ids)], method = "pearson")
  # 增加 margins 并调整 tl.pos 确保文字不被截断
  corrplot(cor_res, method = "color", addCoef.col = "black", type = "upper",
           tl.col = "black", tl.srt = 45, diag = FALSE,
           col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200),
           mar = c(0,0,2,0)) # 顶部预留空间给标题
}

# 步骤2: 将热图转化为一个符合 cowplot 规范的对象
p3_heatmap <- as_grob(draw_corr)

# ==========================================
# 7. 终极排版 (T字型布局，防止覆盖)
# ==========================================
# 下方行：将 B 和 C 并排，设置标签
bottom_row <- plot_grid(p2, p3_heatmap, 
                        labels = c("", "C. Co-expression Matrix"), 
                        label_size = 12, rel_widths = c(1, 1))

# 垂直组合：上方放 A，下方放组合好的 bottom_row
# 通过 rel_heights 增加下方的空间，防止热图文字重叠
final_plot <- plot_grid(p1, bottom_row, 
                        ncol = 1, 
                        rel_heights = c(1, 1.2)) 

# ==========================================
# 8. 高清导出 (防止小窗口预览时的重叠)
# ==========================================
print(final_plot)

# 建议导出为 12x10 英寸，此时 R 渲染引擎会自动调整文字间距，防止重叠
ggsave("Analysis_Result_Final.png", final_plot, width = 12, height = 10, dpi = 300)

