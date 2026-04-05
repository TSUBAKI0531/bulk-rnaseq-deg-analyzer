"""
RNA-seq DEG Analysis & Visualization Tool
==========================================
pyDESeq2 を用いた差次発現遺伝子（DEG）解析パイプライン
ボルカノプロット・ヒートマップ・PCA・GO エンリッチメント解析 を自動出力

必要ライブラリ:
  pip install streamlit pydeseq2 seaborn matplotlib pandas scikit-learn gseapy
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO

# ============================================================
# ページ設定
# ============================================================
st.set_page_config(page_title="RNA-seq DEG Analyzer", layout="wide")
st.title("🧬 RNA-seq DEG Analysis Pipeline")
st.caption("pyDESeq2 による差次発現解析 → ボルカノプロット / ヒートマップ / PCA / GO エンリッチメント 自動出力")

# ============================================================
# サイドバー: パラメータ設定
# ============================================================
st.sidebar.header("⚙️ 解析パラメータ")
padj_threshold = st.sidebar.slider("adjusted p-value 閾値", 0.001, 0.1, 0.05, 0.005)
lfc_threshold = st.sidebar.slider("|log2FC| 閾値", 0.5, 3.0, 1.0, 0.25)
top_n_genes = st.sidebar.slider("ヒートマップ表示遺伝子数", 10, 100, 30, 5)
min_count_filter = st.sidebar.number_input("最小カウントフィルタ（全サンプル合計）", 0, 100, 10)


# ============================================================
# デモデータ生成関数
# ============================================================
def generate_demo_data(n_genes=2000, n_samples=6):
    """
    テスト用のカウントデータとメタデータを生成
    - 前半3サンプル: Control, 後半3サンプル: Treatment
    - 一部遺伝子に意図的な発現差を付与
    """
    np.random.seed(42)
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    sample_names = [f"Sample_{i+1}" for i in range(n_samples)]

    # ベースカウント（負の二項分布でシミュレーション）
    base_mean = np.random.lognormal(mean=5, sigma=1.5, size=n_genes)
    counts = np.zeros((n_genes, n_samples), dtype=int)

    for j in range(n_samples):
        for i in range(n_genes):
            mu = base_mean[i]
            # Treatment群（後半）の一部遺伝子に発現変動を付与
            if j >= n_samples // 2:
                if i < 50:       # Up-regulated genes
                    mu *= np.random.uniform(3, 8)
                elif i < 100:    # Down-regulated genes
                    mu *= np.random.uniform(0.1, 0.3)
            # 負の二項分布からサンプリング
            r = 5  # dispersion
            p = r / (r + mu)
            counts[i, j] = np.random.negative_binomial(r, p)

    count_df = pd.DataFrame(counts, index=gene_names, columns=sample_names)

    metadata_df = pd.DataFrame({
        "sample": sample_names,
        "condition": ["Control"] * (n_samples // 2) + ["Treatment"] * (n_samples // 2)
    })

    return count_df, metadata_df


# ============================================================
# データ入力セクション
# ============================================================
st.header("📂 Step 1: データ入力")

use_demo = st.checkbox("✅ デモデータを使用する（初回はこちらで動作確認）", value=True)

if use_demo:
    count_df, metadata_df = generate_demo_data()
    st.success("デモデータを生成しました（2000遺伝子 × 6サンプル）")

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("カウントマトリクス（先頭5行）")
        st.dataframe(count_df.head(), use_container_width=True)
    with col2:
        st.subheader("メタデータ")
        st.dataframe(metadata_df, use_container_width=True)
else:
    st.info(
        "**カウントマトリクス CSV**: 行=遺伝子, 列=サンプル（インデックスが遺伝子名）\n\n"
        "**メタデータ CSV**: 'sample' 列と 'condition' 列が必要"
    )
    col1, col2 = st.columns(2)
    with col1:
        count_file = st.file_uploader("カウントマトリクス CSV", type=["csv", "tsv"])
    with col2:
        meta_file = st.file_uploader("メタデータ CSV", type=["csv", "tsv"])

    if count_file and meta_file:
        sep = "\t" if count_file.name.endswith(".tsv") else ","
        count_df = pd.read_csv(count_file, index_col=0, sep=sep)
        metadata_df = pd.read_csv(meta_file, sep=sep)
        st.success(f"読込完了: {count_df.shape[0]} 遺伝子 × {count_df.shape[1]} サンプル")
    else:
        st.stop()

# ============================================================
# 条件選択
# ============================================================
conditions = metadata_df["condition"].unique().tolist()
if len(conditions) < 2:
    st.error("メタデータに2群以上の条件が必要です")
    st.stop()

col1, col2 = st.columns(2)
with col1:
    ref_condition = st.selectbox("Reference（対照群）", conditions, index=0)
with col2:
    remaining = [c for c in conditions if c != ref_condition]
    test_condition = st.selectbox("Test（処理群）", remaining, index=0)

# ============================================================
# DEG解析実行
# ============================================================
st.header("🔬 Step 2: DEG 解析")

if st.button("▶ 解析を実行", type="primary", use_container_width=True):

    with st.spinner("pyDESeq2 で差次発現解析を実行中..."):

        # --- 前処理: 低発現遺伝子フィルタ ---
        gene_total = count_df.sum(axis=1)
        count_filtered = count_df[gene_total >= min_count_filter]
        n_filtered = count_df.shape[0] - count_filtered.shape[0]
        st.info(f"低発現フィルタ: {n_filtered} 遺伝子を除外 → {count_filtered.shape[0]} 遺伝子で解析")

        # --- 2群に絞り込み ---
        target_samples = metadata_df[
            metadata_df["condition"].isin([ref_condition, test_condition])
        ]["sample"].tolist()
        count_filtered = count_filtered[[c for c in target_samples if c in count_filtered.columns]]
        meta_subset = metadata_df[metadata_df["sample"].isin(count_filtered.columns)].copy()
        meta_subset = meta_subset.set_index("sample")
        # カラム順をメタデータに揃える
        count_filtered = count_filtered[meta_subset.index]

        # --- pyDESeq2 実行 ---
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats

            dds = DeseqDataSet(
                counts=count_filtered.T,  # pyDESeq2はサンプル×遺伝子
                metadata=meta_subset,
                design="~condition",
            )
            dds.deseq2()

            stats = DeseqStats(
                dds,
                contrast=["condition", test_condition, ref_condition],
            )
            stats.summary()

            results_df = stats.results_df.copy()
            results_df = results_df.dropna(subset=["padj", "log2FoldChange"])

        except ImportError:
            st.error(
                "pydeseq2 がインストールされていません。\n\n"
                "`pip install pydeseq2` を実行してください。"
            )
            st.stop()
        except Exception as e:
            st.error(f"DESeq2 解析中にエラーが発生: {e}")
            st.stop()

    # --- DEG分類 ---
    results_df["deg_status"] = "Not Significant"
    results_df.loc[
        (results_df["padj"] < padj_threshold) & (results_df["log2FoldChange"] > lfc_threshold),
        "deg_status"
    ] = "Up-regulated"
    results_df.loc[
        (results_df["padj"] < padj_threshold) & (results_df["log2FoldChange"] < -lfc_threshold),
        "deg_status"
    ] = "Down-regulated"

    n_up = (results_df["deg_status"] == "Up-regulated").sum()
    n_down = (results_df["deg_status"] == "Down-regulated").sum()

    st.success(f"解析完了！ Up: {n_up} 遺伝子 / Down: {n_down} 遺伝子")

    # セッションに保存
    st.session_state["results_df"] = results_df
    st.session_state["count_filtered"] = count_filtered
    st.session_state["meta_subset"] = meta_subset

# ============================================================
# 可視化セクション
# ============================================================
if "results_df" in st.session_state:
    results_df = st.session_state["results_df"]
    count_filtered = st.session_state["count_filtered"]
    meta_subset = st.session_state["meta_subset"]

    st.header("📊 Step 3: 可視化")

    # ----------------------------------------------------------
    # 3-1. ボルカノプロット
    # ----------------------------------------------------------
    st.subheader("🌋 ボルカノプロット")

    fig_volcano, ax = plt.subplots(figsize=(10, 7))

    color_map = {
        "Up-regulated": "#e74c3c",
        "Down-regulated": "#3498db",
        "Not Significant": "#bdc3c7"
    }

    for status, color in color_map.items():
        subset = results_df[results_df["deg_status"] == status]
        ax.scatter(
            subset["log2FoldChange"],
            -np.log10(subset["padj"]),
            c=color, label=status, alpha=0.6, s=15, edgecolors="none"
        )

    ax.axhline(-np.log10(padj_threshold), color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(lfc_threshold, color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(-lfc_threshold, color="gray", linestyle="--", linewidth=0.8)

    # 上位遺伝子のラベル表示
    sig_genes = results_df[results_df["deg_status"] != "Not Significant"].copy()
    sig_genes["rank_score"] = sig_genes["log2FoldChange"].abs() * (-np.log10(sig_genes["padj"]))
    top_labels = sig_genes.nlargest(10, "rank_score")
    for gene_name, row in top_labels.iterrows():
        ax.annotate(
            gene_name,
            (row["log2FoldChange"], -np.log10(row["padj"])),
            fontsize=7, alpha=0.8,
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
            textcoords="offset points", xytext=(8, 5)
        )

    ax.set_xlabel("log₂ Fold Change", fontsize=12)
    ax.set_ylabel("-log₁₀ (adjusted p-value)", fontsize=12)
    ax.set_title(f"Volcano Plot: {test_condition} vs {ref_condition}", fontsize=14)
    ax.legend(loc="upper right")
    plt.tight_layout()
    st.pyplot(fig_volcano)

    # ----------------------------------------------------------
    # 3-2. MA プロット
    # ----------------------------------------------------------
    st.subheader("📈 MA プロット")
    st.caption(
        "X軸 = 平均発現量（baseMean の log₁₀）、Y軸 = log₂ Fold Change。"
        "ボルカノプロットと併せて見ると、発現量と変動の関係を把握できます。"
    )

    fig_ma, ax_ma = plt.subplots(figsize=(10, 7))

    # baseMean が pyDESeq2 の結果に含まれる
    if "baseMean" in results_df.columns:
        ma_x = np.log10(results_df["baseMean"] + 1)
        ma_x_label = "log₁₀(baseMean + 1)"
    else:
        # baseMean がない場合はカウント平均で代用
        gene_means = count_filtered.mean(axis=1)
        ma_x = np.log10(gene_means.reindex(results_df.index) + 1)
        ma_x_label = "log₁₀(mean count + 1)"

    for status, color in color_map.items():
        subset = results_df[results_df["deg_status"] == status]
        ax_ma.scatter(
            ma_x.loc[subset.index],
            subset["log2FoldChange"],
            c=color, label=status, alpha=0.5, s=12, edgecolors="none"
        )

    # 閾値ライン
    ax_ma.axhline(0, color="black", linewidth=0.8)
    ax_ma.axhline(lfc_threshold, color="gray", linestyle="--", linewidth=0.8)
    ax_ma.axhline(-lfc_threshold, color="gray", linestyle="--", linewidth=0.8)

    # 上位DEG ラベル（ボルカノと同じ上位10遺伝子）
    for gene_name, row in top_labels.iterrows():
        ax_ma.annotate(
            gene_name,
            (ma_x.loc[gene_name], row["log2FoldChange"]),
            fontsize=7, alpha=0.8,
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
            textcoords="offset points", xytext=(8, 5)
        )

    ax_ma.set_xlabel(ma_x_label, fontsize=12)
    ax_ma.set_ylabel("log₂ Fold Change", fontsize=12)
    ax_ma.set_title(f"MA Plot: {test_condition} vs {ref_condition}", fontsize=14)
    ax_ma.legend(loc="upper right")
    plt.tight_layout()
    st.pyplot(fig_ma)

    # ----------------------------------------------------------
    # 3-3. ヒートマップ（上位DEG）
    # ----------------------------------------------------------
    st.subheader("🗺️ クラスターヒートマップ")

    sig_genes_sorted = results_df[
        results_df["deg_status"] != "Not Significant"
    ].reindex(results_df["padj"].sort_values().index).head(top_n_genes)

    if len(sig_genes_sorted) > 0:
        heatmap_genes = sig_genes_sorted.index.tolist()
        heatmap_data = count_filtered.loc[
            count_filtered.index.isin(heatmap_genes)
        ]

        # log2(CPM+1) 正規化
        lib_size = heatmap_data.sum(axis=0)
        cpm = heatmap_data.div(lib_size, axis=1) * 1e6
        log_cpm = np.log2(cpm + 1)

        # Z-score（行方向）
        z_scores = log_cpm.subtract(log_cpm.mean(axis=1), axis=0).div(
            log_cpm.std(axis=1), axis=0
        ).dropna()

        # 条件ラベルのカラーバー
        condition_colors = meta_subset["condition"].map({
            ref_condition: "#3498db",
            test_condition: "#e74c3c"
        })

        fig_heatmap = sns.clustermap(
            z_scores,
            cmap="RdBu_r",
            center=0,
            col_colors=condition_colors,
            figsize=(10, max(6, len(z_scores) * 0.25)),
            dendrogram_ratio=(0.15, 0.1),
            cbar_pos=(0.02, 0.8, 0.03, 0.15),
            yticklabels=True,
            xticklabels=True,
            linewidths=0.5,
        )
        fig_heatmap.ax_heatmap.set_ylabel("")
        fig_heatmap.fig.suptitle(
            f"Top {len(z_scores)} DEGs (Z-score of log₂CPM)",
            y=1.02, fontsize=14
        )
        st.pyplot(fig_heatmap.fig)
    else:
        st.warning("有意なDEGが見つかりませんでした。閾値を緩和してください。")

    # ----------------------------------------------------------
    # 3-4. PCA プロット
    # ----------------------------------------------------------
    st.subheader("📐 PCA プロット")

    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    # log2(CPM+1) → 全遺伝子で PCA
    lib_size_all = count_filtered.sum(axis=0)
    cpm_all = count_filtered.div(lib_size_all, axis=1) * 1e6
    log_cpm_all = np.log2(cpm_all + 1).T  # サンプル×遺伝子

    scaler = StandardScaler()
    scaled = scaler.fit_transform(log_cpm_all)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(scaled)

    pca_df = pd.DataFrame(pcs, columns=["PC1", "PC2"], index=log_cpm_all.index)
    pca_df["condition"] = meta_subset.loc[pca_df.index, "condition"]

    fig_pca, ax_pca = plt.subplots(figsize=(8, 6))
    for cond, color in [(ref_condition, "#3498db"), (test_condition, "#e74c3c")]:
        subset = pca_df[pca_df["condition"] == cond]
        ax_pca.scatter(subset["PC1"], subset["PC2"], c=color, label=cond, s=100, edgecolors="black")
        for idx, row in subset.iterrows():
            ax_pca.annotate(idx, (row["PC1"], row["PC2"]), fontsize=8, textcoords="offset points", xytext=(6, 6))

    ax_pca.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)", fontsize=12)
    ax_pca.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)", fontsize=12)
    ax_pca.set_title("PCA Plot (all genes)", fontsize=14)
    ax_pca.legend()
    plt.tight_layout()
    st.pyplot(fig_pca)

    # ----------------------------------------------------------
    # Step 4: GO エンリッチメント解析 (gseapy)
    # ----------------------------------------------------------
    st.header("🧩 Step 4: GO エンリッチメント解析")

    st.info(
        "**gseapy** の Enrichr API を利用し、DEG リストに対して "
        "Over-Representation Analysis (ORA) を実行します。\n\n"
        "デモデータ（Gene_XXXX 形式）の場合はモック結果を表示します。"
        " 実データでは HUGO 遺伝子シンボル（例: TP53, BRCA1）を使用してください。"
    )

    # --- サイドバーに GO 解析パラメータを追加 ---
    st.sidebar.divider()
    st.sidebar.header("🧩 GO 解析パラメータ")
    go_gene_sets = st.sidebar.multiselect(
        "Gene Set ライブラリ",
        options=[
            "GO_Biological_Process_2023",
            "GO_Molecular_Function_2023",
            "GO_Cellular_Component_2023",
            "KEGG_2021_Human",
            "Reactome_2022",
            "WikiPathway_2023_Human",
        ],
        default=["GO_Biological_Process_2023"],
        help="Enrichr に登録されているライブラリから選択"
    )
    go_top_n = st.sidebar.slider("GO 結果表示数", 5, 30, 10, 5)
    go_cutoff = st.sidebar.slider("GO adjusted p-value 閾値", 0.01, 0.1, 0.05, 0.01)

    # --- DEGリスト抽出 ---
    up_genes = results_df[results_df["deg_status"] == "Up-regulated"].index.tolist()
    down_genes = results_df[results_df["deg_status"] == "Down-regulated"].index.tolist()
    all_deg_genes = up_genes + down_genes

    # --- デモデータ判定 ---
    is_demo_genes = all(g.startswith("Gene_") for g in all_deg_genes[:10]) if all_deg_genes else True

    def generate_mock_go_results(gene_list, direction_label, gene_sets):
        """デモ用: 生物学的にリアルなモック GO 結果を生成"""
        np.random.seed(hash(direction_label) % 2**31)
        mock_terms = {
            "GO_Biological_Process_2023": [
                "immune response (GO:0006955)",
                "inflammatory response (GO:0006954)",
                "cytokine-mediated signaling (GO:0019221)",
                "antigen processing and presentation (GO:0019882)",
                "T cell activation (GO:0042110)",
                "B cell receptor signaling (GO:0050853)",
                "regulation of apoptotic process (GO:0042981)",
                "cell proliferation (GO:0008283)",
                "signal transduction (GO:0007165)",
                "response to virus (GO:0009615)",
                "NF-kappaB signaling (GO:0038061)",
                "JAK-STAT cascade (GO:0007259)",
                "complement activation (GO:0006956)",
                "phagocytosis (GO:0006909)",
                "leukocyte migration (GO:0050900)",
            ],
            "GO_Molecular_Function_2023": [
                "cytokine receptor binding (GO:0005126)",
                "protein kinase activity (GO:0004672)",
                "antigen binding (GO:0003823)",
                "MHC class I binding (GO:0042288)",
                "chemokine activity (GO:0008009)",
                "immunoglobulin receptor binding (GO:0034987)",
                "receptor signaling protein activity (GO:0005057)",
                "GTPase activity (GO:0003924)",
                "ubiquitin ligase activity (GO:0061630)",
                "transcription factor binding (GO:0008134)",
            ],
            "GO_Cellular_Component_2023": [
                "plasma membrane (GO:0005886)",
                "extracellular space (GO:0005615)",
                "MHC class I complex (GO:0042612)",
                "immunological synapse (GO:0001772)",
                "cell surface (GO:0009986)",
                "endoplasmic reticulum (GO:0005783)",
                "proteasome complex (GO:0000502)",
                "cytoplasm (GO:0005737)",
                "nucleus (GO:0005634)",
                "mitochondrion (GO:0005739)",
            ],
            "KEGG_2021_Human": [
                "Cytokine-cytokine receptor interaction",
                "NF-kappa B signaling pathway",
                "TNF signaling pathway",
                "Chemokine signaling pathway",
                "Natural killer cell mediated cytotoxicity",
                "Antigen processing and presentation",
                "B cell receptor signaling pathway",
                "T cell receptor signaling pathway",
                "Toll-like receptor signaling pathway",
                "JAK-STAT signaling pathway",
            ],
            "Reactome_2022": [
                "Immune System R-HSA-168256",
                "Adaptive Immune System R-HSA-1280218",
                "Innate Immune System R-HSA-168249",
                "Cytokine Signaling R-HSA-1280215",
                "Interferon Signaling R-HSA-913531",
                "Interleukin-1 signaling R-HSA-9020702",
                "TCR signaling R-HSA-202403",
                "MHC class I presentation R-HSA-983169",
                "Apoptosis R-HSA-109581",
                "Signaling by Interleukins R-HSA-449147",
            ],
            "WikiPathway_2023_Human": [
                "Inflammatory Response Pathway WP453",
                "Type II interferon signaling WP619",
                "B Cell Receptor Signaling WP23",
                "T Cell Receptor Signaling WP69",
                "Apoptosis WP254",
                "Toll-like Receptor Signaling WP75",
                "Cytokines and Inflammatory Response WP530",
                "NF-kB Signaling Pathway WP4562",
                "Complement and Coagulation Cascades WP558",
                "IL-1 signaling pathway WP195",
            ],
        }

        rows = []
        for gs in gene_sets:
            terms = mock_terms.get(gs, mock_terms["GO_Biological_Process_2023"])
            n_terms = min(go_top_n, len(terms))
            for i in range(n_terms):
                pval = np.random.exponential(0.005) * (i + 1) * 0.3
                padj_val = min(pval * n_terms * 0.5, 1.0)
                overlap_n = max(3, int(np.random.uniform(3, min(20, len(gene_list)))))
                bg_n = np.random.randint(50, 500)
                overlap_genes = gene_list[:overlap_n]
                odds = (overlap_n / len(gene_list)) / (bg_n / 20000) if len(gene_list) > 0 else 1
                rows.append({
                    "Gene_set": gs,
                    "Term": terms[i],
                    "Overlap": f"{overlap_n}/{bg_n}",
                    "P-value": pval,
                    "Adjusted P-value": padj_val,
                    "Odds Ratio": round(odds, 2),
                    "Genes": ";".join(overlap_genes[:5]),
                    "Direction": direction_label,
                })
        df = pd.DataFrame(rows)
        return df

    def run_enrichr(gene_list, direction_label, gene_sets):
        """gseapy.enrich で Enrichr API を呼び出す"""
        import gseapy as gp

        all_results = []
        for gs in gene_sets:
            try:
                enr = gp.enrich(
                    gene_list=gene_list,
                    gene_sets=gs,
                    organism="human",
                    outdir=None,        # ファイル出力なし
                    no_plot=True,
                    cutoff=go_cutoff,
                )
                if enr.results is not None and len(enr.results) > 0:
                    res = enr.results.head(go_top_n).copy()
                    res["Direction"] = direction_label
                    all_results.append(res)
            except Exception as e:
                st.warning(f"⚠️ {gs} の解析中にエラー: {e}")

        if all_results:
            return pd.concat(all_results, ignore_index=True)
        else:
            return pd.DataFrame()

    # --- 解析実行 ---
    go_tabs = st.tabs(["⬆️ Up-regulated", "⬇️ Down-regulated", "🔄 All DEGs"])

    for tab, gene_list, label in zip(
        go_tabs,
        [up_genes, down_genes, all_deg_genes],
        ["Up-regulated", "Down-regulated", "All DEGs"]
    ):
        with tab:
            if len(gene_list) == 0:
                st.warning(f"{label}: 該当する DEG がありません。")
                continue

            st.write(f"**入力遺伝子数:** {len(gene_list)}")

            if is_demo_genes:
                st.warning(
                    "🔶 デモデータ (Gene_XXXX) のため、Enrichr API は使用せずモック結果を表示します。\n\n"
                    "実際の遺伝子シンボル（TP53, IL6 等）を使えば Enrichr API で本物の GO 解析が走ります。"
                )
                go_result_df = generate_mock_go_results(gene_list, label, go_gene_sets)
            else:
                with st.spinner(f"{label} の GO エンリッチメント解析を実行中..."):
                    try:
                        go_result_df = run_enrichr(gene_list, label, go_gene_sets)
                    except ImportError:
                        st.error(
                            "gseapy がインストールされていません。\n\n"
                            "`pip install gseapy` を実行してください。"
                        )
                        continue

            if go_result_df.empty:
                st.info("有意なエンリッチメントは見つかりませんでした。")
                continue

            # --- 結果をセッションに保存 ---
            session_key = f"go_result_{label}"
            st.session_state[session_key] = go_result_df

            # --- 可視化: バープロット ---
            st.markdown(f"##### バープロット（{label}）")

            for gs_name in go_result_df["Gene_set"].unique():
                gs_data = go_result_df[go_result_df["Gene_set"] == gs_name].copy()
                gs_data = gs_data.sort_values("Adjusted P-value").head(go_top_n)

                if gs_data.empty:
                    continue

                fig_bar, ax_bar = plt.subplots(figsize=(10, max(3, len(gs_data) * 0.4)))

                # 短縮表示用のTerm名
                gs_data["Term_short"] = gs_data["Term"].apply(
                    lambda x: (x[:50] + "...") if len(x) > 53 else x
                )

                bars = ax_bar.barh(
                    range(len(gs_data)),
                    -np.log10(gs_data["Adjusted P-value"].values),
                    color="#e74c3c" if "Up" in label else "#3498db" if "Down" in label else "#9b59b6",
                    edgecolor="white",
                    linewidth=0.5,
                )

                ax_bar.set_yticks(range(len(gs_data)))
                ax_bar.set_yticklabels(gs_data["Term_short"].values, fontsize=9)
                ax_bar.set_xlabel("-log₁₀(Adjusted P-value)", fontsize=11)
                ax_bar.set_title(f"{gs_name} — {label}", fontsize=12)
                ax_bar.invert_yaxis()

                # p-value 閾値ライン
                ax_bar.axvline(-np.log10(go_cutoff), color="gray", linestyle="--", linewidth=0.8, label=f"p={go_cutoff}")
                ax_bar.legend(fontsize=8)

                plt.tight_layout()
                st.pyplot(fig_bar)

                # バープロットのダウンロード
                buf_bar = BytesIO()
                fig_bar.savefig(buf_bar, format="png", dpi=150, bbox_inches="tight")
                st.download_button(
                    f"⬇️ {gs_name}_{label}_barplot.png",
                    data=buf_bar.getvalue(),
                    file_name=f"GO_{gs_name}_{label}_barplot.png",
                    mime="image/png",
                )

            # --- 可視化: ドットプロット ---
            st.markdown(f"##### ドットプロット（{label}）")

            for gs_name in go_result_df["Gene_set"].unique():
                gs_data = go_result_df[go_result_df["Gene_set"] == gs_name].copy()
                gs_data = gs_data.sort_values("Adjusted P-value").head(go_top_n)

                if gs_data.empty:
                    continue

                # Overlap からカウント抽出
                gs_data["Overlap_count"] = gs_data["Overlap"].apply(
                    lambda x: int(x.split("/")[0])
                )

                gs_data["Term_short"] = gs_data["Term"].apply(
                    lambda x: (x[:50] + "...") if len(x) > 53 else x
                )

                fig_dot, ax_dot = plt.subplots(figsize=(10, max(3, len(gs_data) * 0.4)))

                scatter = ax_dot.scatter(
                    -np.log10(gs_data["Adjusted P-value"].values),
                    range(len(gs_data)),
                    s=gs_data["Overlap_count"].values * 20,    # ドットサイズ = Overlap数
                    c=-np.log10(gs_data["Adjusted P-value"].values),
                    cmap="YlOrRd",
                    edgecolors="black",
                    linewidth=0.5,
                    alpha=0.8,
                )

                ax_dot.set_yticks(range(len(gs_data)))
                ax_dot.set_yticklabels(gs_data["Term_short"].values, fontsize=9)
                ax_dot.set_xlabel("-log₁₀(Adjusted P-value)", fontsize=11)
                ax_dot.set_title(f"{gs_name} — {label} (dot size = gene count)", fontsize=12)
                ax_dot.invert_yaxis()

                cbar = plt.colorbar(scatter, ax=ax_dot, shrink=0.6, pad=0.02)
                cbar.set_label("-log₁₀(padj)", fontsize=9)

                plt.tight_layout()
                st.pyplot(fig_dot)

            # --- 結果テーブル ---
            st.markdown(f"##### 結果テーブル（{label}）")
            display_go_cols = [c for c in ["Gene_set", "Term", "Overlap", "P-value", "Adjusted P-value", "Genes"]
                              if c in go_result_df.columns]
            st.dataframe(
                go_result_df[display_go_cols].sort_values("Adjusted P-value"),
                use_container_width=True,
                height=300,
            )

            # GO結果 CSVダウンロード
            go_csv_buf = BytesIO()
            go_result_df.to_csv(go_csv_buf, index=False, encoding="utf-8-sig")
            st.download_button(
                f"⬇️ GO結果CSV ({label})",
                data=go_csv_buf.getvalue(),
                file_name=f"go_enrichment_{label}.csv",
                mime="text/csv",
            )

    # ----------------------------------------------------------
    # 結果テーブル & ダウンロード
    # ----------------------------------------------------------
    st.header("📋 Step 5: DEG 結果出力")

    st.subheader("DEG 結果テーブル")
    display_cols = ["log2FoldChange", "pvalue", "padj", "deg_status"]
    st.dataframe(
        results_df[display_cols].sort_values("padj"),
        use_container_width=True,
        height=400
    )

    # CSVダウンロード
    csv_buf = BytesIO()
    results_df.to_csv(csv_buf, encoding="utf-8-sig")
    st.download_button(
        "⬇️ DEG結果をCSVダウンロード",
        data=csv_buf.getvalue(),
        file_name="deg_results.csv",
        mime="text/csv",
        use_container_width=True
    )

    # 図のダウンロード
    col1, col2, col3, col4 = st.columns(4)
    for col, fig_obj, name in [
        (col1, fig_volcano, "volcano_plot.png"),
        (col2, fig_ma, "ma_plot.png"),
        (col3, fig_heatmap.fig if len(sig_genes_sorted) > 0 else None, "heatmap.png"),
        (col4, fig_pca, "pca_plot.png")
    ]:
        if fig_obj is not None:
            buf = BytesIO()
            fig_obj.savefig(buf, format="png", dpi=150, bbox_inches="tight")
            col.download_button(
                f"⬇️ {name}",
                data=buf.getvalue(),
                file_name=name,
                mime="image/png",
                use_container_width=True
            )

# ============================================================
# フッター
# ============================================================
st.divider()
st.caption(
    "RNA-seq DEG Analysis Pipeline | "
    "Powered by pyDESeq2 + Streamlit | "
    "GitHub: [your-username]/rnaseq-deg-analyzer"
)
