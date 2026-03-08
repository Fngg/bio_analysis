# -*- coding: utf-8 -*-
"""
GSE75214 基因表达数据分析
- 差异基因分析
- 火山图、热图
- GO/KEGG 富集分析
"""

import os
import sys
import gzip

os.makedirs("data", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


def download_series_matrix():
    """下载 series matrix 并解析"""
    import urllib.request
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75214/matrix/GSE75214_series_matrix.txt.gz"
    path = "data/GSE75214_series_matrix.txt.gz"
    if not os.path.exists(path):
        print("  下载 series matrix...")
        urllib.request.urlretrieve(url, path)
    return path


def parse_series_matrix(path):
    """解析 series matrix 文件"""
    meta_lines = []
    table_lines = []
    in_table = False
    with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if line.startswith('!'):
                if not in_table:
                    meta_lines.append(line)
                if 'series_matrix_table_begin' in line:
                    in_table = True
            elif 'series_matrix_table_end' in line:
                break
            elif in_table:
                table_lines.append(line)

    def strip_quotes(s):
        return str(s).strip().strip('"')

    # 解析 metadata：Sample_geo_accession 和 Sample_title
    samples = []
    titles = []
    for line in meta_lines:
        if line.startswith('!Sample_geo_accession'):
            parts = line.split('\t')[1:]
            samples = [strip_quotes(p) for p in parts]
        elif line.startswith('!Sample_title'):
            parts = line.split('\t')[1:]
            titles = [strip_quotes(p) for p in parts]

    sample_meta = {}
    for i, s in enumerate(samples):
        sample_meta[s] = {'group': 'Unknown'}
        text = (titles[i] if i < len(titles) else '') + ' '
        # 从 title 推断分组，如 CD_colon_active_1, UC_colon_inactive_1, Control_1
        if 'control' in text.lower() or 'non-ibd' in text.lower() or 'healthy' in text.lower():
            sample_meta[s]['group'] = 'Control'
        elif 'cd_colon_active' in text.lower() or ('crohn' in text.lower() and 'active' in text.lower()):
            sample_meta[s]['group'] = 'CD_active'
        elif 'cd_' in text.lower() or 'crohn' in text.lower():
            sample_meta[s]['group'] = 'CD'
        elif 'uc_colon_active' in text.lower() or ('uc' in text.lower() and 'active' in text.lower()):
            sample_meta[s]['group'] = 'UC_active'
        elif 'uc_' in text.lower() or 'ulcerative' in text.lower():
            sample_meta[s]['group'] = 'UC_inactive'
        elif sample_meta[s]['group'] == 'Unknown':
            sample_meta[s]['group'] = 'Other'

    # 解析表达矩阵
    if not table_lines:
        raise ValueError("No expression table found")
    header = [strip_quotes(x) for x in table_lines[0].split('\t')]
    rows = []
    for line in table_lines[1:]:
        parts = line.split('\t')
        if len(parts) >= len(header):
            rows.append([strip_quotes(p) for p in parts[:len(header)]])
    df = pd.DataFrame(rows, columns=header)
    df = df.set_index(header[0])
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df = df.dropna(how='all')
    return df, sample_meta


def load_gpl_annotation(gpl_id="GPL6244"):
    """加载平台注释（使用 GEOparse）"""
    try:
        import GEOparse
        gpl = GEOparse.get_GEO(geo=gpl_id, destdir="data")
        if "gene_assignment" in gpl.table.columns:
            def extract_symbol(s):
                if pd.isna(s) or s == "---":
                    return ""
                parts = str(s).split(" /// ")[0].split(" // ")
                return parts[1].strip() if len(parts) >= 2 else ""
            annot = gpl.table[["ID"]].copy()
            annot["gene_symbol"] = gpl.table["gene_assignment"].apply(extract_symbol)
        elif "Gene Symbol" in gpl.table.columns:
            annot = gpl.table[["ID", "Gene Symbol"]].rename(columns={"Gene Symbol": "gene_symbol", "ID": "probe"})
        else:
            return pd.DataFrame()
        annot = annot.rename(columns={"ID": "probe"})
        annot = annot[(annot["gene_symbol"].notna()) & (annot["gene_symbol"] != "") & (annot["gene_symbol"] != "---")]
        return annot
    except Exception as e:
        print(f"  GPL 注释加载失败: {e}")
        return pd.DataFrame()


def main():
    print("=" * 60)
    print("GSE75214 分析开始")
    print("=" * 60)

    # 1. 下载并解析
    print("\n[1/7] 下载 GSE75214 数据...")
    path = download_series_matrix()
    expr, sample_meta = parse_series_matrix(path)
    print(f"  表达矩阵: {expr.shape[0]} 探针 x {expr.shape[1]} 样本")

    meta = pd.DataFrame.from_dict(sample_meta, orient='index')
    meta = meta[meta.index.isin(expr.columns)]
    print(meta['group'].value_counts())

    # 选择分组
    case_label, ctrl_label = "UC_active", "Control"
    case_samples = meta.index[meta["group"] == case_label].intersection(expr.columns)
    ctrl_samples = meta.index[meta["group"] == ctrl_label].intersection(expr.columns)

    if len(case_samples) < 2 or len(ctrl_samples) < 2:
        case_label, ctrl_label = "CD_active", "Control"
        case_samples = meta.index[meta["group"] == case_label].intersection(expr.columns)
        ctrl_samples = meta.index[meta["group"] == ctrl_label].intersection(expr.columns)

    if len(case_samples) < 2 or len(ctrl_samples) < 2:
        case_label, ctrl_label = "CD", "Control"
        case_samples = meta.index[meta["group"] == case_label].intersection(expr.columns)
        ctrl_samples = meta.index[meta["group"] == ctrl_label].intersection(expr.columns)

    print(f"  {case_label}: {len(case_samples)} 样本")
    print(f"  {ctrl_label}: {len(ctrl_samples)} 样本")

    expr = expr[case_samples.tolist() + ctrl_samples.tolist()]

    # 2. 差异分析
    print("\n[3/7] 差异表达分析...")
    if expr.min().min() > 0 and expr.max().max() > 100:
        expr = np.log2(expr + 1)
        print("  已进行 log2 转换")

    results = []
    for probe in expr.index:
        case_vals = expr.loc[probe, case_samples]
        ctrl_vals = expr.loc[probe, ctrl_samples]
        t_stat, p_val = stats.ttest_ind(case_vals, ctrl_vals, equal_var=False, nan_policy="omit")
        log2fc = case_vals.mean() - ctrl_vals.mean()
        results.append({
            "probe": probe,
            "log2FC": log2fc,
            "pvalue": p_val,
            "mean_case": case_vals.mean(),
            "mean_ctrl": ctrl_vals.mean(),
        })

    deg = pd.DataFrame(results).dropna(subset=["pvalue"])
    deg["padj"] = multipletests(deg["pvalue"], method="fdr_bh")[1]

    padj_cut, lfc_cut = 0.05, 1.0
    deg["sig"] = (deg["padj"] < padj_cut) & (deg["log2FC"].abs() >= lfc_cut)
    print(f"  显著差异探针数: {deg['sig'].sum()}")

    # 3. 探针注释
    print("\n[4/7] 探针注释...")
    annot = load_gpl_annotation()
    if annot.empty:
        annot = pd.DataFrame({'probe': expr.index, 'gene_symbol': expr.index})
    else:
        annot['gene_symbol'] = annot['gene_symbol'].astype(str).str.split(' /// ').str[0].str.strip()
        annot = annot[(annot['gene_symbol'].notna()) & (annot['gene_symbol'] != '') & (annot['gene_symbol'] != 'nan')]

    deg["probe"] = deg["probe"].astype(str)
    annot["probe"] = annot["probe"].astype(str)
    deg_annot = deg.merge(annot, on="probe", how="left")
    deg_annot = deg_annot[(deg_annot["gene_symbol"].notna()) & (deg_annot["gene_symbol"] != "")]

    sig_deg = deg_annot[deg_annot["sig"]].copy()
    up_genes = sig_deg[sig_deg["log2FC"] > 0]["gene_symbol"].drop_duplicates().tolist()
    down_genes = sig_deg[sig_deg["log2FC"] < 0]["gene_symbol"].drop_duplicates().tolist()
    print(f"  上调基因: {len(up_genes)}, 下调基因: {len(down_genes)}")

    deg_annot.to_csv("results/GSE75214_DEG_gene_level.csv", index=False)

    # 4. 火山图
    print("\n[5/7] 绘制火山图...")
    volcano = deg_annot.copy()
    volcano["-log10_padj"] = -np.log10(np.maximum(volcano["padj"].values, 1e-300))

    fig, ax = plt.subplots(figsize=(8, 6))
    is_sig = volcano["sig"]
    ax.scatter(volcano.loc[~is_sig, "log2FC"], volcano.loc[~is_sig, "-log10_padj"],
               s=8, c="grey", alpha=0.4, label="Not significant")
    ax.scatter(volcano.loc[is_sig, "log2FC"], volcano.loc[is_sig, "-log10_padj"],
               s=8, c="red", alpha=0.7, label="Significant")
    ax.axhline(-np.log10(padj_cut), color="blue", linestyle="--", alpha=0.6)
    ax.axvline(lfc_cut, color="black", linestyle="--", alpha=0.6)
    ax.axvline(-lfc_cut, color="black", linestyle="--", alpha=0.6)
    ax.set_xlabel("log2FC")
    ax.set_ylabel("-log10(padj)")
    ax.set_title(f"GSE75214 Volcano: {case_label} vs {ctrl_label}")
    ax.legend()
    plt.tight_layout()
    plt.savefig("figures/GSE75214_volcano.png", dpi=300)
    plt.close()
    print("  已保存 figures/GSE75214_volcano.png")

    # 5. 热图
    print("\n[6/7] 绘制热图...")
    topN = min(50, len(sig_deg))
    if topN > 0:
        top_genes = sig_deg.sort_values("padj").drop_duplicates("gene_symbol").head(topN)["gene_symbol"].tolist()
        probe_to_gene = annot.set_index("probe")["gene_symbol"].to_dict() if not annot.empty else {}
        expr_gene = expr.copy()
        expr_gene["gene_symbol"] = expr_gene.index.map(lambda x: probe_to_gene.get(x, x))
        expr_gene = expr_gene.dropna(subset=["gene_symbol"])
        expr_gene = expr_gene.groupby("gene_symbol").mean()

        available = [g for g in top_genes if g in expr_gene.index]
        if available:
            heatmat = expr_gene.loc[available, case_samples.tolist() + ctrl_samples.tolist()]
            std = heatmat.std(axis=1).replace(0, 1e-8)
            heatmat_z = (heatmat - heatmat.mean(axis=1).values.reshape(-1, 1)) / std.values.reshape(-1, 1)
            col_colors = pd.Series(index=heatmat_z.columns, dtype=str)
            col_colors[case_samples] = "#e74c3c"
            col_colors[ctrl_samples] = "#3498db"

            g = sns.clustermap(heatmat_z, cmap="RdBu_r", center=0, col_colors=col_colors,
                               xticklabels=False, yticklabels=True, figsize=(10, 10))
            g.fig.suptitle(f"GSE75214 Top{len(available)} DEG Heatmap", y=1.02)
            plt.savefig("figures/GSE75214_heatmap.png", dpi=300, bbox_inches="tight")
            plt.close()
            print("  已保存 figures/GSE75214_heatmap.png")

    # 6. GO/KEGG 富集
    print("\n[7/7] GO/KEGG 富集分析...")
    try:
        import gget
        max_genes = 200
        all_sig_genes = up_genes + down_genes
        genes_for_enrich = list(dict.fromkeys(all_sig_genes))[:max_genes]

        if len(genes_for_enrich) >= 5:
            go_res = gget.enrichr(genes_for_enrich, database="ontology")
            kegg_res = gget.enrichr(genes_for_enrich, database="pathway")

            if go_res is not None and len(go_res) > 0:
                go_res.to_csv("results/GSE75214_GO_enrichment.csv", index=False)
                print(f"  GO 富集: {len(go_res)} 条结果 -> results/GSE75214_GO_enrichment.csv")
            if kegg_res is not None and len(kegg_res) > 0:
                kegg_res.to_csv("results/GSE75214_KEGG_enrichment.csv", index=False)
                print(f"  KEGG 富集: {len(kegg_res)} 条结果 -> results/GSE75214_KEGG_enrichment.csv")
        else:
            print("  显著基因数不足，跳过富集分析")
    except Exception as e:
        print(f"  富集分析跳过: {e}")

    print("\n" + "=" * 60)
    print("分析完成")
    print("=" * 60)
    return deg_annot, sig_deg


if __name__ == "__main__":
    main()
