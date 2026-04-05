# Streamlit Cloud デプロイ手順書

## 前提条件

- VS Code がインストール済み
- GitHub アカウントを持っている
- VS Code で GitHub にサインイン済み

---

## Step 1: ローカルでの動作確認

ターミナル（VS Code 内の Terminal でOK）で以下を実行：

```bash
# プロジェクトフォルダを作成
mkdir rnaseq-deg-analyzer
cd rnaseq-deg-analyzer

# ダウンロードしたファイルをこのフォルダに配置
# 配置するファイル一覧:
#   rnaseq_deg_app.py
#   requirements.txt
#   packages.txt
#   README.md
#   .gitignore
#   .streamlit/config.toml

# 仮想環境を作成（Anaconda の場合）
conda create -n rnaseq python=3.10 -y
conda activate rnaseq

# または venv の場合
# python -m venv venv
# source venv/bin/activate  (Mac/Linux)
# venv\Scripts\activate     (Windows)

# ライブラリインストール
pip install -r requirements.txt

# ローカルで起動テスト
streamlit run rnaseq_deg_app.py
```

ブラウザで `http://localhost:8501` が開き、デモデータで動作すれば OK。

---

## Step 2: GitHub リポジトリの作成（VS Code から）

### 方法 A: VS Code のソース管理 UI を使う（推奨）

1. VS Code でプロジェクトフォルダを開く
2. 左サイドバーの **ソース管理**（Source Control）アイコンをクリック
3. **「リポジトリを初期化」** をクリック
4. すべてのファイルをステージング（`+` ボタン）
5. コミットメッセージ: `Initial commit: RNA-seq DEG analyzer with GO enrichment`
6. **「Publish Branch」** をクリック → GitHub にパブリックリポジトリとして公開

### 方法 B: ターミナルから

```bash
cd rnaseq-deg-analyzer

git init
git add .
git commit -m "Initial commit: RNA-seq DEG analyzer with GO enrichment"

# GitHub CLI がある場合
gh repo create rnaseq-deg-analyzer --public --push

# GitHub CLI がない場合は GitHub.com で手動でリポジトリ作成後:
# git remote add origin https://github.com/YOUR_USERNAME/rnaseq-deg-analyzer.git
# git branch -M main
# git push -u origin main
```

---

## Step 3: Streamlit Cloud にデプロイ

### 3-1. Streamlit Cloud にサインアップ

1. ブラウザで https://share.streamlit.io にアクセス
2. **「Sign up」** → **「Continue with GitHub」** を選択
3. GitHub アカウントで認証

### 3-2. アプリをデプロイ

1. ダッシュボードで **「New app」** ボタンをクリック
2. 以下を入力：

   | 項目 | 入力値 |
   |------|--------|
   | Repository | `YOUR_USERNAME/rnaseq-deg-analyzer` |
   | Branch | `main` |
   | Main file path | `rnaseq_deg_app.py` |

3. **「Deploy!」** をクリック
4. 初回デプロイは 3〜5 分かかる（pydeseq2 のビルドに時間がかかるため）

### 3-3. デプロイ完了

成功すると以下のような URL が発行される：

```
https://your-username-rnaseq-deg-analyzer-rnaseq-deg-app-xxxxx.streamlit.app
```

この URL を職務経歴書やポートフォリオに記載できる。

---

## Step 4: カスタム URL の設定（任意）

1. Streamlit Cloud ダッシュボードでアプリの **「Settings」** を開く
2. **「General」** タブ → **「App URL」** で短い URL に変更可能

例: `https://rnaseq-deg-analyzer.streamlit.app`

---

## Step 5: README のバッジ URL を更新

`README.md` の以下の行を、実際のデプロイ URL に書き換える：

```markdown
[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://your-actual-app-url.streamlit.app)
```

変更後、VS Code から commit & push すれば GitHub 上の README も更新される。

---

## トラブルシューティング

### デプロイが失敗する場合

1. **Streamlit Cloud のログを確認**
   - ダッシュボード → アプリ → **「Manage app」** → ログを表示
   - エラーメッセージを確認

2. **よくあるエラーと対処法**

   | エラー | 対処法 |
   |--------|--------|
   | `ModuleNotFoundError: pydeseq2` | `requirements.txt` に `pydeseq2>=0.4.0` があるか確認 |
   | ビルドエラー（C拡張） | `packages.txt` に `build-essential` があるか確認 |
   | メモリ不足 | デモデータの遺伝子数を減らす（2000 → 1000） |
   | タイムアウト | Streamlit Cloud 無料枠のリソース上限。遺伝子数を減らして対応 |

3. **Enrichr API 接続エラー**
   - Streamlit Cloud からは外部 API（Enrichr）にアクセス可能
   - 一時的なネットワークエラーの場合はリトライで解消されることが多い

### ローカルでは動くがクラウドで動かない場合

- `requirements.txt` のバージョン指定を確認
- Python バージョンの差異: Streamlit Cloud は Python 3.10+ を使用
- `packages.txt` でシステムライブラリが必要な場合がある

---

## リポジトリの最終的なファイル構成

```
rnaseq-deg-analyzer/
├── .streamlit/
│   └── config.toml          # テーマ・サーバー設定
├── .gitignore                # Git 除外ファイル
├── README.md                 # プロジェクト説明
├── packages.txt              # システムレベルの依存パッケージ
├── requirements.txt          # Python ライブラリ
└── rnaseq_deg_app.py         # メインアプリケーション
```
