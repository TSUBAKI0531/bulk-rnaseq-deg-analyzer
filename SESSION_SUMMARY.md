# スレッド要約: RNA-seq DEG Analysis Pipeline 開発記録

**作成日:** 2026年3月26日  
**テーマ:** RNA-seq 発現解析 & ボルカノプロット + ヒートマップ自動生成ツール  
**最終成果物:** Streamlit Web アプリ（774行）+ デプロイ一式 + 転職用言語化ガイド

---

## 全体の流れ（6ターン）

| # | ユーザーの依頼 | 実施内容 | 主な成果物 |
|---|--------------|---------|-----------|
| 1 | テーマ提示 + たたき台コードを作りたい | 処理フロー設計 → MVP コード作成 → 職務経歴書記載案 | `rnaseq_deg_app.py`（初版）, `requirements.txt` |
| 2 | GO エンリッチメント解析を追加したい | gseapy 統合、モック結果生成、バー/ドットプロット追加 | アプリに Step 4 を追加（約200行増） |
| 3 | MA plot を追加したい | ボルカノとヒートマップの間に MA plot を挿入 | アプリに 3-2 セクション追加（約50行増） |
| 4 | Streamlit Cloud へのデプロイを進めたい | デプロイ用設定ファイル一式 + 手順書を作成 | config.toml, packages.txt, .gitignore, DEPLOY_GUIDE.md |
| 5 | README を作成したい | ポートフォリオ品質の README に全面刷新 | README.md（313行）, LICENSE, docs/screenshots/ |
| 6 | 転職用にどう言語化すべきか | 職務経歴書・面接Q&A・スキルシート・ポートフォリオ戦略 | CAREER_GUIDE.md（192行） |

---

## ターン 1: MVP コード作成

### 設計した処理フロー

```
カウントCSV + メタデータCSV
  → 低発現フィルタ（pandas）
  → DEG検定（pyDESeq2: DeseqDataSet → DeseqStats）
  → ボルカノプロット（matplotlib）
  → ヒートマップ（seaborn clustermap, Z-score正規化）
  → PCA（scikit-learn）
  → CSV + PNG ダウンロード
```

### コードの主要設計判断

- **デモデータ生成機能を内蔵**: 負の二項分布で2000遺伝子×6サンプルのカウントデータを生成。先頭50遺伝子をUp、次の50をDownに意図的に設定し、DEGが検出される状態を保証
- **pyDESeq2 の入力形式**: 一般的なカウントマトリクス（遺伝子×サンプル）を `.T` で転置してpyDESeq2 に渡す処理を実装
- **ヒートマップの正規化**: log2(CPM+1) → Z-score（行方向）の2段階正規化。サンプル間のライブラリサイズ差を補正した上で、遺伝子ごとの相対変動を可視化
- **ボルカノプロットの上位ラベル**: |log2FC| × -log10(padj) のランクスコアで上位10遺伝子を自動ラベル表示
- **Streamlit のセッション管理**: `st.session_state` に解析結果を保存し、ボタン再押下なしで可視化セクションを維持

### 使用ライブラリと役割

| ライブラリ | 役割 |
|-----------|------|
| pydeseq2 | DESeq2のPython実装。負の二項分布モデルによるDEG検定 |
| seaborn | clustermap で階層クラスタリング付きヒートマップを1行で描画 |
| matplotlib | ボルカノプロット・PCA等のカスタム描画 |
| scikit-learn | PCA（主成分分析）+ StandardScaler |
| pandas | CSV読込・データ操作の全基盤 |
| streamlit | Web UIフレームワーク。サイドバー・ボタン・ダウンロード等 |

---

## ターン 2: GO エンリッチメント解析の追加

### 追加した処理フロー

```
DEG結果 → Up / Down / All のリスト抽出
  → gseapy.enrich() で Enrichr API に送信
  → 6ライブラリ（GO_BP, GO_MF, GO_CC, KEGG, Reactome, WikiPathway）
  → バープロット + ドットプロット + テーブル出力
```

### 主要な実装判断

- **デモデータ判定ロジック**: 遺伝子名が `Gene_XXXX` 形式かどうかを先頭10遺伝子で判定。デモ時はEnrichr APIを呼ばず、生物学的にリアルなモック結果（immune response, cytokine signaling 等の実在GOタームを使用）を返す
- **Up / Down / All のタブ切替**: `st.tabs()` で3つのタブを生成し、各方向のDEGリストを個別に解析。方向別のパスウェイ解釈を可能にした
- **ドットプロットのドットサイズ**: Overlap数（遺伝子セットとの重複数）に比例させることで、統計的有意性と生物学的関連性の両方を1つの図で表現
- **サイドバーへのパラメータ追加**: `st.sidebar.divider()` でDEG解析パラメータとGO解析パラメータを視覚的に分離

### gseapy の使い方

```python
import gseapy as gp
enr = gp.enrich(
    gene_list=gene_list,       # DEG遺伝子名のリスト
    gene_sets="GO_Biological_Process_2023",  # ライブラリ名
    organism="human",
    outdir=None,               # ファイル出力なし（メモリ上で完結）
    no_plot=True,              # gseapy内蔵プロットは使わず独自描画
    cutoff=0.05,               # adjusted p-value 閾値
)
results_df = enr.results      # pandas DataFrame で結果取得
```

---

## ターン 3: MA plot の追加

### MA plot の位置づけ

- **X軸**: log10(baseMean + 1) — pyDESeq2が算出する正規化済み平均発現量
- **Y軸**: log2 Fold Change — ボルカノプロットと同じ
- **ボルカノとの補完関係**: ボルカノは「何が有意か」、MAは「発現量と変動の関係（低発現遺伝子の偽陽性チェック）」

### 実装のポイント

- `baseMean` カラムが pyDESeq2 結果にあればそれを使用、なければカウント平均で代用するフォールバック処理
- ボルカノプロットと同じ `color_map` と `top_labels`（上位10遺伝子）を共有し、図の一貫性を確保
- ダウンロードセクションを3カラム→4カラムに拡張（Volcano / MA / Heatmap / PCA）

---

## ターン 4: Streamlit Cloud デプロイ準備

### 作成したファイル一覧

| ファイル | 用途 |
|---------|------|
| `.streamlit/config.toml` | テーマ色（primaryColor: #e74c3c）、アップロード上限200MB |
| `packages.txt` | Streamlit Cloud用のシステムライブラリ（build-essential, gfortran, libopenblas-dev） |
| `.gitignore` | __pycache__, .env, venv 等の除外 |
| `DEPLOY_GUIDE.md` | VS Code → GitHub → Streamlit Cloud の手順書（トラブルシューティング付き） |

### packages.txt が必要な理由

pyDESeq2 は内部で numpy/scipy の C 拡張をビルドするため、Streamlit Cloud のLinux環境で `build-essential`（gccコンパイラ）と `libopenblas-dev`（線形代数ライブラリ）が必要。

---

## ターン 5: README 作成

### 旧版からの主な改善点

| 追加セクション | 目的 |
|--------------|------|
| Demo（スクリーンショット6枚テーブル） | GitHub を開いた瞬間の視覚的訴求 |
| Architecture 図（ASCII） | 技術的全体像の一覧化。面接での説明材料 |
| Background（なぜPythonか） | 技術選定の意図を明文化 → 設計判断力のアピール |
| Roadmap（チェックボックス） | 将来の拡張計画 → 成長の方向性を示す |
| Docker デプロイ例 | コンテナ化の知識もアピール |
| References | 原著論文リンクで学術的裏付け |
| バッジ（Python / Streamlit / License） | 一目で技術スタックが分かる |

### スクリーンショットの追加手順

アプリをローカルで起動 → デモデータで解析実行 → 各図の「⬇️ PNGダウンロード」で保存 → `docs/screenshots/` に配置 → commit & push。

---

## ターン 6: 転職用言語化

### CAREER_GUIDE.md の構成

1. **職務経歴書 記載例 3パターン** — 詳細版（バイオインフォ専任向け）/ バランス版（ウェット兼任向け）/ 簡潔版（1〜3行）
2. **面接 Q&A 8問** — プロジェクト説明、なぜPythonか、pyDESeq2とは、GOの解釈、ボルカノ vs MA、PCAの見方、苦労した点、今後の改善
3. **スキルシート記載例** — スキル × レベル × 実績 の対応表
4. **ポートフォリオ全体での位置づけ** — 既存プロジェクト（SDM-Primer, GlycoAntibody, CAR-T等）との関係性。RNA-seqツールが「ウェット→ドライの転換点」として機能
5. **求人要件との対応表** — 典型的な求人キーワードと本プロジェクトの該当箇所
6. **書類通過率を上げる Tips** — 避けるべき表現（×「AIで作った」→ ○「Pythonで設計・実装」）

---

## 最終成果物一覧

```
rnaseq-deg-analyzer/
├── .streamlit/config.toml       # テーマ・サーバー設定
├── .gitignore                   # Git 除外設定
├── LICENSE                      # MIT License
├── README.md                    # ポートフォリオ品質 README（313行）
├── packages.txt                 # Streamlit Cloud 用システムライブラリ
├── requirements.txt             # Python 依存パッケージ（8ライブラリ）
├── rnaseq_deg_app.py            # メインアプリケーション（774行）
├── docs/screenshots/.gitkeep    # スクリーンショット格納先
├── DEPLOY_GUIDE.md              # デプロイ手順書（175行）※リポジトリ非公開用
└── CAREER_GUIDE.md              # 転職用言語化ガイド（192行）※リポジトリ非公開用
```

### アプリの最終機能構成

```
Step 1: データ入力（CSV アップロード / デモデータ）
Step 2: DEG 解析（pyDESeq2）
Step 3: 可視化
  ├── 3-1 ボルカノプロット（上位10遺伝子ラベル付き）
  ├── 3-2 MA プロット（baseMean vs log2FC）
  ├── 3-3 クラスターヒートマップ（Z-score, 階層クラスタリング）
  └── 3-4 PCA プロット（寄与率表示）
Step 4: GO エンリッチメント解析
  ├── Up / Down / All DEG タブ切替
  ├── 6ライブラリ対応（GO_BP/MF/CC, KEGG, Reactome, WikiPathway）
  ├── バープロット + ドットプロット
  └── デモデータ時はモック結果を自動生成
Step 5: 結果出力
  ├── DEG テーブル（インタラクティブ表示）
  ├── CSV ダウンロード（DEG結果 / GO結果）
  └── PNG ダウンロード（4図 + GOプロット）
```

### 使用技術の全体像

```
Python 3.10+ / Streamlit / pyDESeq2 / gseapy / seaborn / matplotlib
scikit-learn / pandas / numpy / Git / GitHub / Streamlit Cloud
```
