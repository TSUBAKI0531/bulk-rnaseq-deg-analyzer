# 🧬 RNA-seq DEG Analysis Pipeline

> **Python-only RNA-seq differential expression analysis pipeline with automated visualization and pathway enrichment.**

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://your-app-url.streamlit.app)
![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![Streamlit](https://img.shields.io/badge/Streamlit-1.30%2B-FF4B4B?logo=streamlit&logoColor=white)
![License](https://img.shields.io/badge/License-MIT-green)

RNA-seq カウントデータから差次発現遺伝子（DEG）の統計検定、多角的な可視化、GO/KEGG パスウェイ解析までを **R を使わず Python だけで完結** する Web アプリケーションです。  
デモデータ内蔵のため、実データがなくても即座に全機能を試すことができます。

---

## Demo

<!-- デプロイ後、実際のスクリーンショットに差し替えてください -->

| Volcano Plot | MA Plot |
|:---:|:---:|
| ![Volcano](docs/screenshots/volcano_plot.png) | ![MA](docs/screenshots/ma_plot.png) |

| Cluster Heatmap | PCA Plot |
|:---:|:---:|
| ![Heatmap](docs/screenshots/heatmap.png) | ![PCA](docs/screenshots/pca_plot.png) |

| GO Enrichment Bar Plot | GO Enrichment Dot Plot |
|:---:|:---:|
| ![GO Bar](docs/screenshots/go_barplot.png) | ![GO Dot](docs/screenshots/go_dotplot.png) |

> **💡 スクリーンショットの追加方法:**  
> アプリをローカルで起動 → 各図を PNG ダウンロード → `docs/screenshots/` フォルダに配置してください。

---

## Features

### 🔬 DEG 検定 — pyDESeq2

R の DESeq2 と同等の負の二項分布モデルを Python で実行します。低発現遺伝子の自動フィルタリング、正規化、統計検定を一括処理します。

- 2群比較（Control vs Treatment）
- adjusted p-value / log2 Fold Change の閾値をサイドバーからリアルタイム調整
- 全結果を CSV でダウンロード可能

### 📊 4種の可視化を自動出力

| 図 | 何がわかるか |
|---|---|
| **Volcano Plot** | どの遺伝子が有意に変動したか（統計的有意性 × 変動量） |
| **MA Plot** | 発現量と変動量の関係（低発現遺伝子の偽陽性チェックに有用） |
| **Cluster Heatmap** | 上位 DEG の発現パターン（サンプル・遺伝子の階層クラスタリング付き） |
| **PCA Plot** | サンプル間の全体的な類似性（品質評価・外れ値検出） |

すべての図は上位 DEG の遺伝子名ラベルを自動表示し、PNG でダウンロードできます。

### 🧩 GO エンリッチメント解析 — gseapy

DEG リストに対して Over-Representation Analysis (ORA) を実行し、生物学的な意味づけを行います。

- **6つの Gene Set ライブラリ**に対応:
  GO Biological Process / Molecular Function / Cellular Component, KEGG, Reactome, WikiPathway
- **Up / Down / All DEG をタブで切替**して方向別に解析
- バープロット + ドットプロットで結果を可視化
- Enrichr API 経由で最新のアノテーションデータを使用

### ⚙️ インタラクティブなパラメータ調整

サイドバーから以下をリアルタイムに変更可能です:

- adjusted p-value 閾値 (0.001〜0.1)
- |log2FC| 閾値 (0.5〜3.0)
- ヒートマップ表示遺伝子数 (10〜100)
- 低発現遺伝子フィルタ閾値
- GO 解析: ライブラリ選択、表示数、p-value 閾値

---

## Architecture

```
┌─────────────────────────────────────────────────┐
│                  Streamlit UI                    │
│  ┌───────────┐  ┌──────────┐  ┌──────────────┐  │
│  │ File      │  │ Sidebar  │  │ Results &    │  │
│  │ Upload /  │  │ Params   │  │ Download     │  │
│  │ Demo Data │  │ Control  │  │ (CSV + PNG)  │  │
│  └─────┬─────┘  └────┬─────┘  └──────▲───────┘  │
│        │              │               │          │
│  ┌─────▼──────────────▼───────────────┤          │
│  │         Analysis Engine            │          │
│  │  ┌──────────┐  ┌───────────────┐   │          │
│  │  │ pyDESeq2 │  │ scikit-learn  │   │          │
│  │  │ DEG Test │  │ PCA           │   │          │
│  │  └────┬─────┘  └───────┬───────┘   │          │
│  │       │                │           │          │
│  │  ┌────▼────────────────▼────────┐  │          │
│  │  │   Visualization Layer        │  │          │
│  │  │  matplotlib / seaborn        │  │          │
│  │  │  Volcano, MA, Heatmap, PCA   │  │          │
│  │  └──────────────────────────────┘  │          │
│  │                                    │          │
│  │  ┌──────────────────────────────┐  │          │
│  │  │  gseapy (Enrichr API)        │──┘          │
│  │  │  GO / KEGG / Reactome        │             │
│  │  └──────────────────────────────┘             │
│  └───────────────────────────────────────────────┘
└─────────────────────────────────────────────────┘
```

---

## Quick Start

### ローカル実行

```bash
# 1. リポジトリをクローン
git clone https://github.com/your-username/rnaseq-deg-analyzer.git
cd rnaseq-deg-analyzer

# 2. 仮想環境を作成（Anaconda の場合）
conda create -n rnaseq python=3.10 -y
conda activate rnaseq

# 3. ライブラリインストール
pip install -r requirements.txt

# 4. アプリ起動
streamlit run rnaseq_deg_app.py
```

ブラウザで `http://localhost:8501` が自動で開きます。  
「デモデータを使用する」にチェックが入った状態で **▶ 解析を実行** をクリックすれば、全機能を確認できます。

### Streamlit Cloud で試す

デプロイ済みのアプリはこちらから利用できます:  
👉 **[https://your-app-url.streamlit.app](https://your-app-url.streamlit.app)**

---

## Input Format

本ツールは **2つの CSV ファイル** を入力として受け付けます。

### 1. カウントマトリクス CSV

行 = 遺伝子、列 = サンプルの整数カウント値。先頭列が遺伝子名（インデックス）です。  
HUGO 遺伝子シンボル（TP53, BRCA1 等）を使用すると GO エンリッチメント解析が実データで動作します。

```csv
,Sample_1,Sample_2,Sample_3,Sample_4,Sample_5,Sample_6
TP53,120,135,118,450,520,480
BRCA1,80,92,75,210,245,198
IL6,45,52,48,320,380,355
MYC,200,180,210,95,88,102
```

### 2. メタデータ CSV

`sample` 列（カウントマトリクスのカラム名と一致）と `condition` 列を含みます。

```csv
sample,condition
Sample_1,Control
Sample_2,Control
Sample_3,Control
Sample_4,Treatment
Sample_5,Treatment
Sample_6,Treatment
```

> **📝 Note:** TSV（タブ区切り）にも対応しています。拡張子が `.tsv` であれば自動判定します。

---

## Analysis Pipeline

```
Step 1: データ入力
  │  CSV アップロード or デモデータ（2000 遺伝子 × 6 サンプル）
  ▼
Step 2: DEG 解析 (pyDESeq2)
  │  低発現フィルタ → 正規化 → 負の二項分布モデル → Wald 検定
  │  → log2FoldChange, pvalue, padj を算出
  ▼
Step 3: 可視化（4種の図を自動生成）
  │  ├── 🌋 Volcano Plot   — 有意な DEG を色分け表示
  │  ├── 📈 MA Plot         — 発現量 vs 変動量の関係
  │  ├── 🗺️ Cluster Heatmap — 上位 DEG の Z-score ヒートマップ
  │  └── 📐 PCA Plot        — サンプル間の品質チェック
  ▼
Step 4: GO エンリッチメント解析 (gseapy → Enrichr API)
  │  ├── ⬆️ Up-regulated DEGs
  │  ├── ⬇️ Down-regulated DEGs
  │  └── 🔄 All DEGs
  │  各方向で バープロット + ドットプロット + テーブル を出力
  ▼
Step 5: 結果出力
     ├── DEG 結果テーブル（インタラクティブ表示）
     ├── CSV ダウンロード（DEG 結果 / GO 結果）
     └── PNG ダウンロード（Volcano / MA / Heatmap / PCA / GO 図）
```

---

## Tech Stack

| カテゴリ | ライブラリ | 役割 |
|---------|-----------|------|
| Web UI | **Streamlit** | インタラクティブ Web アプリケーション |
| DEG 検定 | **pyDESeq2** | DESeq2 の Python 実装（負の二項分布モデル） |
| パスウェイ解析 | **gseapy** | Enrichr API 経由の GO / KEGG / Reactome 解析 |
| 可視化 | **matplotlib** | Volcano Plot, MA Plot, PCA Plot |
| 可視化 | **seaborn** | 階層クラスタリング付きヒートマップ (clustermap) |
| 次元削減 | **scikit-learn** | PCA（主成分分析） |
| データ処理 | **pandas / numpy** | カウントデータの操作・正規化・統計量計算 |

---

## Project Structure

```
rnaseq-deg-analyzer/
├── .streamlit/
│   └── config.toml           # Streamlit テーマ・サーバー設定
├── docs/
│   └── screenshots/          # スクリーンショット格納用（任意）
├── .gitignore
├── LICENSE
├── README.md
├── packages.txt              # Streamlit Cloud 用システムライブラリ
├── requirements.txt          # Python 依存パッケージ
└── rnaseq_deg_app.py         # メインアプリケーション（約 770 行）
```

---

## Deployment

### Streamlit Cloud（推奨）

1. このリポジトリを GitHub にプッシュ
2. [share.streamlit.io](https://share.streamlit.io) に GitHub アカウントでログイン
3. 「**New app**」→ リポジトリ / ブランチ (`main`) / メインファイル (`rnaseq_deg_app.py`) を指定
4. 「**Deploy!**」をクリック

初回デプロイは pyDESeq2 のビルドに 3〜5 分かかります。  
`packages.txt` により必要なシステムライブラリが自動インストールされます。

### Docker（オプション）

```dockerfile
FROM python:3.10-slim
WORKDIR /app
RUN apt-get update && apt-get install -y build-essential gfortran libopenblas-dev
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
EXPOSE 8501
CMD ["streamlit", "run", "rnaseq_deg_app.py", "--server.address", "0.0.0.0"]
```

```bash
docker build -t rnaseq-deg-analyzer .
docker run -p 8501:8501 rnaseq-deg-analyzer
```

---

## Roadmap

- [ ] GSEA（Gene Set Enrichment Analysis）のランクベース解析対応
- [ ] 多群比較（3群以上）への拡張
- [ ] バッチ効果補正（ComBat 等）の統合
- [ ] Plotly によるインタラクティブ図への切替
- [ ] シングルセル RNA-seq 対応（scanpy 統合）

---

## Background

### なぜ Python で RNA-seq 解析なのか

RNA-seq の差次発現解析は長年 R（DESeq2, edgeR）が標準でしたが、近年 **pyDESeq2** の登場により Python 単独での解析が実用レベルになりました。Python で統一するメリットは以下の通りです:

- 機械学習・深層学習パイプラインとのシームレスな統合
- Streamlit / Dash 等による Web アプリ化が容易
- バイオインフォマティクスと創薬研究の橋渡し

本ツールは、ウェットラボ研究者がドライ解析を自力で行えるよう、**コードを書かずに GUI で操作できる**ことを設計方針としています。

---

## References

- [pyDESeq2](https://github.com/owkin/PyDESeq2) — Owkin による DESeq2 の Python 実装
- [gseapy](https://github.com/zqfang/GSEApy) — Enrichr / GSEA の Python ラッパー
- [DESeq2 原著論文](https://doi.org/10.1186/s13059-014-0550-8) — Love, Huber & Anders (2014)
- [Enrichr](https://maayanlab.cloud/Enrichr/) — Ma'ayan Lab の遺伝子セット解析ツール

---

## License

MIT License — 詳細は [LICENSE](LICENSE) を参照してください。

## Author

**[Your Name]**  
Bioinformatics / Antibody Drug Development  
GitHub: [@your-username](https://github.com/your-username)
