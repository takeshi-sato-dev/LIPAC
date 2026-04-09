# Synthetic Validation for Mixture Causal Model

This directory contains synthetic validation experiments for the Bayesian mixture causal model. The validation tests the model's ability to correctly classify scenarios and recover true parameter values across seven synthetic datasets.

## Scenarios

| Scenario | Type | β_direct | β_coop | π | σ | Description |
|----------|------|----------|--------|-----|-----|-------------|
| S1 | Linear | -5.0 | 0.0 | 0.0 | 4.0 | Negative direct effect only |
| S2 | Cooperative | -5.0 | 20.0 | 0.10 | 4.0 | Negative direct + strong cooperative |
| S3 | Linear | +3.0 | 0.0 | 0.0 | 4.0 | Positive direct effect only |
| S4 | Cooperative | +3.0 | 5.0 | 0.15 | 4.0 | Positive direct + moderate cooperative |
| S5 | Cooperative | 0.0 | 15.0 | 0.20 | 4.0 | No direct effect, cooperative only |
| S6 | Linear | -10.0 | 0.0 | 0.0 | 4.0 | Large negative direct effect |
| S7 | Cooperative | -1.3 | 19.4 | 0.16 | 4.0 | Small direct + strong cooperative |

Each scenario generates data for 4 protein copies × 3000 frames with ~15% bound fraction.

## Usage

Run all validation scenarios:

```bash
cd synthetic_validation
python run_synthetic.py
```

This will:
1. Generate synthetic data for each scenario
2. Fit both linear and mixture models
3. Compute WAIC and classify each scenario
4. Report parameter recovery errors
5. Save results to `synthetic_validation_results/`

## Output

- `validation_summary.csv`: Overall classification accuracy and parameter recovery
- Individual scenario directories with MCMC traces and diagnostics

## Expected Results

The mixture model should:
- Correctly classify all linear scenarios (ΔWAIC ≤ 2)
- Correctly classify all cooperative scenarios (ΔWAIC > 2, β_coop 95% CI excludes zero)
- Recover β_coop and π within ±2.0 and ±0.05, respectively, for cooperative scenarios
