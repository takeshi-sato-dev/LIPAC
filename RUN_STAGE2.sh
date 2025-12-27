#!/bin/bash
# Stage 2 Hierarchical Causal Analysis実行スクリプト
# 全ての脂質（DIPC, CHOL, DPSM, etc.）についてhierarchical modelを生成

echo "=================================="
echo "STAGE 2 HIERARCHICAL CAUSAL ANALYSIS"
echo "=================================="

cd "$(dirname "$0")/stage2_contact_analysis"

# 実行前の確認
echo ""
echo "実行内容："
echo "  - Individual causal analysis (各タンパク質ごと)"
echo "  - Hierarchical causal analysis (集団レベル)"
echo "  - 全脂質type（DIPC, CHOL, DPSM, etc.）について実行"
echo ""
echo "出力先："
echo "  - MAIN_FIGURE_*: ../stage1_contact_analysis/output/bayesian_analysis/"
echo "  - Hierarchical結果: ../stage1_contact_analysis/output/bayesian_analysis/hierarchical_analysis/"
echo ""

# Python環境確認
echo "Python環境確認..."
python --version || python3 --version

echo ""
echo "必要なモジュール確認..."
python -c "import pandas; import pymc; import arviz; print('✓ All required modules available')" 2>/dev/null || \
python3 -c "import pandas; import pymc; import arviz; print('✓ All required modules available')" 2>/dev/null || \
{
    echo "エラー: pandas, pymc, arvizがインストールされていません"
    echo "LIPACを実行した時と同じPython環境を使用してください"
    exit 1
}

echo ""
echo "Stage 2を実行します..."
echo ""

# Stage 2実行
python main.py || python3 main.py

# 結果確認
echo ""
echo "=================================="
echo "実行完了"
echo "=================================="
echo ""
echo "生成されたファイル："
echo ""
echo "【MAIN_FIGUREファイル】"
ls -lh ../stage1_contact_analysis/output/bayesian_analysis/MAIN_FIGURE* 2>/dev/null || echo "  生成されませんでした"
echo ""
echo "【Hierarchical analysis結果】"
ls -lh ../stage1_contact_analysis/output/bayesian_analysis/hierarchical_analysis/*mu_beta* 2>/dev/null || echo "  生成されませんでした"
echo ""
echo "【詳細】"
ls ../stage1_contact_analysis/output/bayesian_analysis/hierarchical_analysis/ 2>/dev/null || echo "  ディレクトリが作成されませんでした"
