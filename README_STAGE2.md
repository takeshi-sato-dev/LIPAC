# Stage 2 Hierarchical Causal Analysis

## 実行方法

```bash
cd .
./RUN_STAGE2.sh
```

## 生成されるファイル

### 1. MAIN_FIGUREファイル（既存機能、更新される）
場所：`stage1_contact_analysis/output/bayesian_analysis/`

- `MAIN_FIGURE_lipac_causal_effects_comparison_contacts.png` / `.svg`
- `MAIN_FIGURE_lipac_causal_effects_comparison_unique_molecules.png` / `.svg`

### 2. Hierarchical Analysis結果（新機能）
場所：`stage1_contact_analysis/output/bayesian_analysis/hierarchical_analysis/`

**各脂質について以下が生成されます：**

- **DIPC_contacts_mu_beta_trace.png / .svg** - μβのトレースプロット
- **DIPC_contacts_forest.png / .svg** - 個別タンパク質効果
- **DIPC_contacts_full_trace.png / .svg** - 全パラメータのトレース
- **CHOL_contacts_mu_beta_trace.png / .svg**
- **DPSM_contacts_mu_beta_trace.png / .svg**
- **DIPC_unique_mu_beta_trace.png / .svg**
- **CHOL_unique_mu_beta_trace.png / .svg**
- **DPSM_unique_mu_beta_trace.png / .svg**

**サマリー：**
- **hierarchical_summary.csv** - 全脂質の収束診断（r̂, ESS）を含む
- **hierarchical_summary.txt** - 人間が読みやすい形式

### 3. Tracesデータ
場所：`stage1_contact_analysis/output/bayesian_analysis/hierarchical_analysis/traces/`

- `DIPC_contacts_hierarchical_trace.nc`
- `CHOL_contacts_hierarchical_trace.nc`
- ...（各脂質について）

## 収束診断

### EphA2の場合（期待される結果）
- **r̂ < 1.05** (収束成功)
- **ESS > 400** (十分な有効サンプル)
- トレースプロット：4本のchainが重なって見える

### Notchの場合（期待される結果）
- **r̂ > 3.0** (収束失敗)
- **ESS < 10** (不十分なサンプル)
- トレースプロット：chainが発散している

## トラブルシューティング

### エラー: ModuleNotFoundError: No module named 'pandas'

LIPACを実行した時と同じPython環境を使用してください。

### MAIN_FIGUREが更新されない

Causal analysisがスキップされています。以下を確認：

1. `stage1_contact_analysis/output/`にdpg3_causal_data_*.csvファイルが存在するか
2. Stage2実行時に "Causal Bayesian analysis" が実行されたか（ログ確認）

### Hierarchical_analysisディレクトリが作成されない

1. Pythonモジュール（pymc, arviz）がインストールされているか確認
2. エラーログを確認（ログの最後の方）
