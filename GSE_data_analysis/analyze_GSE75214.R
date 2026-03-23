# -*- coding: utf-8 -*-
#
# GSE75214 基因表达数据分析（R 版本）
# 1) 读取/下载 GEO 数据（series matrix）
# 2) 根据样本 title 推断分组（Control / UC_active / CD_active / CD ...）
# 3) 每个探针做 Welch t-test 得到 log2FC、pvalue，并做 BH 校正
# 4) 探针到 gene symbol 注释（GPL6244）
# 5) 火山图、热图
# 6) GO/KEGG 富集分析（clusterProfiler）
#
# 依赖包：
#   GEOquery, ggplot2, dplyr, pheatmap, clusterProfiler, org.Hs.eg.db
#
# 运行示例：
#   Rscript analyze_GSE75214.R

rm(list = ls())

safe_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

safe_require("GEOquery")
safe_require("ggplot2")
safe_require("dplyr")
safe_require("pheatmap")
safe_require("clusterProfiler")
safe_require("org.Hs.eg.db")

library(GEOquery)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# -------- 基础路径 --------
script_dir <- tryCatch({
  normalizePath(dirname(sys.frame(1)$ofile))
}, error = function(e) {
  getwd()
})
setwd(script_dir)

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

cat(strrep("=", 60), "\n", sep = "")
cat("GSE75214 分析开始\n")
cat(strrep("=", 60), "\n", sep = "")

infer_group <- function(title) {
  t <- tolower(ifelse(is.na(title), "", title))
  if (grepl("control", t) || grepl("non-ibd", t) || grepl("healthy", t)) {
    return("Control")
  }
  if (grepl("cd_colon_active", t) || (grepl("crohn", t) && grepl("active", t))) {
    return("CD_active")
  }
  if (grepl("cd_", t) || grepl("crohn", t)) {
    return("CD")
  }
  if (grepl("uc_colon_active", t) || (grepl("uc", t) && grepl("active", t))) {
    return("UC_active")
  }
  if (grepl("uc_", t) || grepl("ulcerative", t)) {
    return("UC_inactive")
  }
  return("Other")
}

extract_titles_from_pdata <- function(expr_set) {
  pdata <- pData(expr_set)
  # 常见列：title、source_name_ch1、characteristics_ch1...
  if ("title" %in% colnames(pdata)) return(as.character(pdata$title))
  if ("source_name_ch1" %in% colnames(pdata)) return(as.character(pdata$source_name_ch1))
  # fallback：用样本名
  return(rownames(pdata))
}

load_expression_matrix <- function(gse_id = "GSE75214") {
  # GEOquery 会自动下载并解析 series matrix
  gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE, AnnotGPL = FALSE)
  # GSE 可能是列表（多个平台），这里取第一个 ExpressionSet
  if (is.list(gse)) expr_set <- gse[[1]] else expr_set <- gse
  expr <- exprs(expr_set)
  expr <- as.matrix(expr)
  return(list(expr = expr, expr_set = expr_set))
}

load_gpl_annotation <- function(gpl_id = "GPL6244") {
  gpl <- getGEO(gpl_id, AnnotGPL = FALSE, getGPL = TRUE)
  gpl_tab <- Table(gpl)

  # python 里优先 gene_assignment，再解析取第一个/第二段里的 gene symbol
  # GEOquery 表字段名可能随版本不同，这里做一些兼容：
  # 典型列：ID、gene_assignment 或 Gene Symbol
  colnames_lower <- tolower(colnames(gpl_tab))

  if ("gene_assignment" %in% colnames_lower) {
    # 映射列名：找到实际列名
    ga_col <- colnames(gpl_tab)[which(colnames_lower == "gene_assignment")][1]
    id_col <- ifelse("id" %in% colnames_lower, colnames(gpl_tab)[which(colnames_lower == "id")][1], "ID")

    extract_symbol <- function(s) {
      if (is.na(s) || s == "---") return("")
      s <- as.character(s)
      # 复刻 python：parts = str(s).split(" /// ")[0].split(" // ")
      parts <- strsplit(s, " /// ", fixed = TRUE)[[1]][1]
      parts2 <- strsplit(parts, " // ", fixed = TRUE)[[1]]
      if (length(parts2) >= 2) return(trimws(parts2[2])) else return("")
    }

    annot <- data.frame(
      probe = as.character(gpl_tab[[id_col]]),
      gene_symbol = vapply(gpl_tab[[ga_col]], extract_symbol, character(1)),
      stringsAsFactors = FALSE
    )
    annot <- annot %>%
      filter(!is.na(gene_symbol), gene_symbol != "", gene_symbol != "nan")
    return(annot)
  }

  if ("gene symbol" %in% colnames(gpl_tab)) {
    annot <- data.frame(
      probe = as.character(gpl_tab$ID),
      gene_symbol = as.character(gpl_tab$`Gene Symbol`),
      stringsAsFactors = FALSE
    )
    annot <- annot %>%
      filter(!is.na(gene_symbol), gene_symbol != "", gene_symbol != "---")
    colnames(annot)[1] <- "probe"
    return(annot)
  }

  # fallback：空注释
  return(data.frame(probe = rownames(gpl_tab), gene_symbol = NA_character_, stringsAsFactors = FALSE)[0, ])
}

dl <- load_expression_matrix("GSE75214")
expr <- dl$expr
expr_set <- dl$expr_set

cat("[1/7] 表达矩阵维度：", nrow(expr), "探针 x ", ncol(expr), "样本\n", sep = "")

titles <- extract_titles_from_pdata(expr_set)
sample_names <- colnames(expr)
groups <- sapply(titles, infer_group)

meta <- data.frame(
  sample = sample_names,
  title = titles,
  group = groups,
  stringsAsFactors = FALSE
)
meta <- meta %>% filter(sample %in% colnames(expr))

cat(meta$group %>% table() %>% as.character(), "\n", sep = "")

case_label <- "UC_active"
ctrl_label <- "Control"
case_samples <- meta %>% filter(group == case_label) %>% pull(sample) %>% intersect(colnames(expr))
ctrl_samples <- meta %>% filter(group == ctrl_label) %>% pull(sample) %>% intersect(colnames(expr))

if (length(case_samples) < 2 || length(ctrl_samples) < 2) {
  case_label <- "CD_active"
  case_samples <- meta %>% filter(group == case_label) %>% pull(sample) %>% intersect(colnames(expr))
}
if (length(case_samples) < 2 || length(ctrl_samples) < 2) {
  case_label <- "CD"
  case_samples <- meta %>% filter(group == case_label) %>% pull(sample) %>% intersect(colnames(expr))
}

if (length(case_samples) < 2 || length(ctrl_samples) < 2) {
  stop("无法找到足够的 case/ctrl 样本（至少各 2 个）。")
}

cat("选择分组：", case_label, " vs ", ctrl_label, "\n", sep = "")
cat("case 样本数：", length(case_samples), "\n", sep = "")
cat("ctrl 样本数：", length(ctrl_samples), "\n", sep = "")

expr <- expr[, c(case_samples, ctrl_samples), drop = FALSE]
case_samples <- c(case_samples)
ctrl_samples <- c(ctrl_samples)

cat("[2/7] 差异表达分析：Welch t-test...\n")
expr_num <- apply(expr, 2, as.numeric)
rownames(expr_num) <- rownames(expr)
colnames(expr_num) <- colnames(expr)

# python 判断：min>0 && max>100 则 log2(expr+1)
if (min(expr_num, na.rm = TRUE) > 0 && max(expr_num, na.rm = TRUE) > 100) {
  expr_num <- log2(expr_num + 1)
  cat("已进行 log2 转换\n")
}

probe_ids <- rownames(expr_num)
results <- vector("list", length(probe_ids))

for (i in seq_along(probe_ids)) {
  probe <- probe_ids[i]
  case_vals <- expr_num[probe, case_samples, drop = TRUE]
  ctrl_vals <- expr_num[probe, ctrl_samples, drop = TRUE]

  # Welch t-test
  tt <- tryCatch(t.test(case_vals, ctrl_vals, var.equal = FALSE), error = function(e) NULL)
  if (is.null(tt)) {
    results[[i]] <- data.frame(
      probe = probe,
      log2FC = NA_real_,
      pvalue = NA_real_,
      mean_case = NA_real_,
      mean_ctrl = NA_real_
    )
  } else {
    results[[i]] <- data.frame(
      probe = probe,
      log2FC = mean(case_vals, na.rm = TRUE) - mean(ctrl_vals, na.rm = TRUE),
      pvalue = tt$p.value,
      mean_case = mean(case_vals, na.rm = TRUE),
      mean_ctrl = mean(ctrl_vals, na.rm = TRUE)
    )
  }
}

deg <- bind_rows(results) %>% filter(!is.na(pvalue))
deg$padj <- p.adjust(deg$pvalue, method = "BH")

padj_cut <- 0.05
lfc_cut <- 1.0
deg$sig <- (deg$padj < padj_cut) & (abs(deg$log2FC) >= lfc_cut)
cat("显著差异探针数：", sum(deg$sig), "\n", sep = "")

cat("[3/7] 探针注释：GPL6244...\n")
annot <- load_gpl_annotation("GPL6244")

if (nrow(annot) == 0) {
  annot <- data.frame(probe = probe_ids, gene_symbol = probe_ids, stringsAsFactors = FALSE)
}

deg_annot <- deg %>%
  mutate(probe = as.character(probe)) %>%
  left_join(annot %>% mutate(probe = as.character(probe)), by = "probe") %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

sig_deg <- deg_annot %>% filter(sig)
up_genes <- sig_deg %>% filter(log2FC > 0) %>% pull(gene_symbol) %>% unique() %>% as.character()
down_genes <- sig_deg %>% filter(log2FC < 0) %>% pull(gene_symbol) %>% unique() %>% as.character()

cat("上调基因：", length(up_genes), " 下调基因：", length(down_genes), "\n", sep = "")

write.csv(deg_annot, file = file.path("results", "GSE75214_DEG_gene_level.csv"), row.names = FALSE)

cat("[4/7] 火山图...\n")
volcano <- deg_annot %>%
  mutate(
    neglog10_padj = -log10(pmax(padj, 1e-300)),
    sig_label = ifelse(sig, "Significant", "Not significant")
  )

vol_plot <- ggplot(volcano, aes(x = log2FC, y = neglog10_padj, color = sig_label)) +
  geom_point(size = 1.3, alpha = 0.8) +
  scale_color_manual(values = c("Not significant" = "grey70", "Significant" = "red")) +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "blue", alpha = 0.6) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black", alpha = 0.6) +
  labs(x = "log2FC", y = "-log10(padj)", title = paste0("GSE75214 Volcano: ", case_label, " vs ", ctrl_label)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = file.path("figures", "GSE75214_volcano.png"), plot = vol_plot, dpi = 300, width = 8, height = 6)
cat("已保存 figures/GSE75214_volcano.png\n")

cat("[5/7] 热图...\n")
topN <- min(50, length(unique(sig_deg$gene_symbol)))

if (topN > 0) {
  top_genes <- sig_deg %>%
    arrange(padj) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    slice_head(n = topN) %>%
    pull(gene_symbol)

  probe_to_gene <- annot %>% select(probe, gene_symbol)
  probe_gene_map <- probe_to_gene$gene_symbol[match(rownames(expr_num), probe_to_gene$probe)]

    valid <- !is.na(probe_gene_map) & probe_gene_map != ""
    expr_valid <- expr_num[valid, , drop = FALSE]              # probes x samples
    gene_valid <- probe_gene_map[valid]                        # length = n_probes_valid

    # 按 gene_symbol 对探针做均值合并，得到 gene x samples
    expr_valid_df <- data.frame(
      Gene = gene_valid,
      expr_valid,
      check.names = FALSE
    )
    expr_gene_df <- expr_valid_df %>%
      group_by(Gene) %>%
      summarise(across(where(is.numeric), mean))

    expr_gene_mat <- as.matrix(expr_gene_df[, -1, drop = FALSE])
    rownames(expr_gene_mat) <- expr_gene_df$Gene
    colnames(expr_gene_mat) <- colnames(expr_valid)

  available <- intersect(top_genes, rownames(expr_gene_mat))
  if (length(available) > 0) {
    heatmat <- expr_gene_mat[available, c(case_samples, ctrl_samples), drop = FALSE]

    # Z-score：每行标准化
    heatmat_z <- t(scale(t(heatmat), center = TRUE, scale = TRUE))
    heatmat_z[is.na(heatmat_z)] <- 0

    group_vec <- c(rep("Case", length(case_samples)), rep("Ctrl", length(ctrl_samples)))
    names(group_vec) <- c(case_samples, ctrl_samples)

    ann_col <- data.frame(Group = group_vec)
    rownames(ann_col) <- colnames(heatmat_z)
    ann_colors <- list(Group = c("Case" = "#e74c3c", "Ctrl" = "#3498db"))

    png(file.path("figures", "GSE75214_heatmap.png"), width = 1000, height = 1000, res = 300)
    pheatmap(
      heatmat_z,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_colnames = FALSE,
      fontsize_row = 8,
      annotation_col = ann_col,
      annotation_colors = ann_colors,
      border_color = NA
    )
    dev.off()
    cat("已保存 figures/GSE75214_heatmap.png\n")
  } else {
    cat("热图：可用 top genes 为 0，跳过热图\n")
  }
}

cat("[6/7] GO/KEGG 富集分析...\n")
max_genes <- 200
all_sig_genes <- unique(c(up_genes, down_genes))
genes_for_enrich <- head(all_sig_genes, max_genes)

if (length(genes_for_enrich) >= 5) {
  # symbol -> entrez
  gene_map <- bitr(
    genes_for_enrich,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  entrez_ids <- unique(gene_map$ENTREZID)
  if (length(entrez_ids) >= 5) {
    go_res <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
      write.csv(as.data.frame(go_res), file = file.path("results", "GSE75214_GO_enrichment.csv"), row.names = FALSE)
      cat("GO 富集：", nrow(as.data.frame(go_res)), " 条结果 -> results/GSE75214_GO_enrichment.csv\n", sep = "")
    }

    kegg_res <- enrichKEGG(
      gene = entrez_ids,
      organism = "hsa",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    )
    if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
      write.csv(as.data.frame(kegg_res), file = file.path("results", "GSE75214_KEGG_enrichment.csv"), row.names = FALSE)
      cat("KEGG 富集：", nrow(as.data.frame(kegg_res)), " 条结果 -> results/GSE75214_KEGG_enrichment.csv\n", sep = "")
    }
  } else {
    cat("富集分析：基因转 Entrez 后数量不足，跳过\n")
  }
} else {
  cat("显著基因数不足（<5），跳过富集分析\n")
}

cat(strrep("=", 60), "\n", sep = "")
cat("分析完成\n")
cat(strrep("=", 60), "\n", sep = "")

