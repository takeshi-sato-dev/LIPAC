"""
Synthetic Validation for Mixture Causal Model

Generates synthetic data for 7 scenarios matching Table 1 from the manuscript
and validates the mixture model's ability to recover true parameters.

Author: Takeshi Sato
Kyoto Pharmaceutical University
"""

import sys
import os
import numpy as np
import pandas as pd

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from stage2_contact_analysis.analysis.mixture_causal_analysis import analyze_lipid_type


def generate_synthetic_data(scenario_name, beta_direct, beta_coop, pi, sigma,
                            n_proteins=4, n_frames=3000, bound_fraction=0.15,
                            sigma_alpha=2.0, sigma_beta_direct=1.0, random_seed=42):
    """Generate synthetic data for one scenario

    Parameters
    ----------
    scenario_name : str
        Name of scenario (e.g., 'S1', 'S2')
    beta_direct : float
        Direct binding effect
    beta_coop : float
        Cooperative effect (0 for linear scenarios)
    pi : float
        Mixing probability (0 for linear scenarios)
    sigma : float
        Observation noise
    n_proteins : int
        Number of protein copies
    n_frames : int
        Number of frames per protein
    bound_fraction : float
        Fraction of frames in bound state
    sigma_alpha : float
        SD of per-protein intercepts
    sigma_beta_direct : float
        SD of per-protein direct effects
    random_seed : int
        Random seed

    Returns
    -------
    pd.DataFrame
        Synthetic data with columns: frame, protein_copy, binding_state, lipid_type, contact_count
    """
    rng = np.random.default_rng(random_seed)

    # Generate per-protein parameters
    mu_alpha = 50.0  # baseline contact count
    alpha_p = rng.normal(mu_alpha, sigma_alpha, n_proteins)

    if sigma_beta_direct > 0:
        beta_direct_p = rng.normal(beta_direct, sigma_beta_direct, n_proteins)
    else:
        beta_direct_p = np.full(n_proteins, beta_direct)

    # Generate data
    data = []

    for p in range(n_proteins):
        # Generate binding states
        n_bound = int(n_frames * bound_fraction)
        n_unbound = n_frames - n_bound

        T = np.concatenate([np.zeros(n_unbound), np.ones(n_bound)])
        rng.shuffle(T)

        # Generate contact counts
        Y = np.zeros(n_frames)

        for i in range(n_frames):
            if T[i] == 0:
                # Unbound
                Y[i] = rng.normal(alpha_p[p], sigma)
            else:
                # Bound
                if pi == 0:
                    # Linear model
                    mu = alpha_p[p] + beta_direct_p[p]
                    Y[i] = rng.normal(mu, sigma)
                else:
                    # Mixture model
                    if rng.random() < pi:
                        # Cooperative component
                        mu = alpha_p[p] + beta_direct_p[p] + beta_coop
                        Y[i] = rng.normal(mu, sigma)
                    else:
                        # Direct-only component
                        mu = alpha_p[p] + beta_direct_p[p]
                        Y[i] = rng.normal(mu, sigma)

        # Store data
        for i in range(n_frames):
            data.append({
                'frame': i,
                'protein_copy': f'protein_{p}',
                'binding_state': int(T[i]),
                'lipid_type': scenario_name,
                'contact_count': Y[i]
            })

    df = pd.DataFrame(data)
    return df


def run_scenario(scenario_name, beta_direct, beta_coop, pi, sigma, output_dir):
    """Run one synthetic validation scenario

    Parameters
    ----------
    scenario_name : str
        Name of scenario
    beta_direct : float
        True direct effect
    beta_coop : float
        True cooperative effect
    pi : float
        True mixing probability
    sigma : float
        True observation noise
    output_dir : str
        Output directory

    Returns
    -------
    dict
        Validation results
    """
    print(f"\n{'='*60}")
    print(f"Scenario {scenario_name}")
    print(f"{'='*60}")
    print(f"True parameters:")
    print(f"  β_direct = {beta_direct:.2f}")
    print(f"  β_coop = {beta_coop:.2f}")
    print(f"  π = {pi:.2f}")
    print(f"  σ = {sigma:.2f}")

    # Generate synthetic data
    df = generate_synthetic_data(
        scenario_name=scenario_name,
        beta_direct=beta_direct,
        beta_coop=beta_coop,
        pi=pi,
        sigma=sigma
    )

    # Analyze with mixture model
    result = analyze_lipid_type(df, scenario_name, output_dir)

    # Compute recovery errors
    beta_coop_error = abs(result['beta_coop_mean'] - beta_coop)
    pi_error = abs(result['pi_mean'] - pi)

    # Check classification
    expected_class = 'cooperative' if beta_coop != 0 else 'linear'
    correct_classification = result['classification'] == expected_class

    print(f"\nRecovery:")
    print(f"  β_coop estimate: {result['beta_coop_mean']:.2f} (error: {beta_coop_error:.2f})")
    print(f"  π estimate: {result['pi_mean']:.2f} (error: {pi_error:.2f})")
    print(f"  Classification: {result['classification']} (expected: {expected_class})")
    print(f"  Correct: {correct_classification}")

    return {
        'scenario': scenario_name,
        'true_beta_direct': beta_direct,
        'true_beta_coop': beta_coop,
        'true_pi': pi,
        'true_sigma': sigma,
        'est_beta_coop': result['beta_coop_mean'],
        'est_pi': result['pi_mean'],
        'est_beta_direct': result['beta_direct_mean'],
        'beta_coop_error': beta_coop_error,
        'pi_error': pi_error,
        'delta_waic': result['delta_waic'],
        'classification': result['classification'],
        'expected_class': expected_class,
        'correct_classification': correct_classification,
        'linear_converged': result['linear_converged'],
        'mixture_converged': result['mixture_converged']
    }


def main():
    """Run all synthetic validation scenarios"""
    # Define scenarios (Table 1)
    scenarios = [
        # Scenario, β_direct, β_coop, π, σ
        ('S1', -5.0, 0.0, 0.0, 4.0),    # Linear, negative
        ('S2', -5.0, 20.0, 0.10, 4.0),  # Cooperative, negative direct
        ('S3', +3.0, 0.0, 0.0, 4.0),    # Linear, positive
        ('S4', +3.0, 5.0, 0.15, 4.0),   # Cooperative, positive direct
        ('S5', 0.0, 15.0, 0.20, 4.0),   # Cooperative, no direct
        ('S6', -10.0, 0.0, 0.0, 4.0),   # Linear, large negative
        ('S7', -1.3, 19.4, 0.16, 4.0)   # Cooperative, small direct
    ]

    # Create output directory
    output_dir = 'synthetic_validation_results'
    os.makedirs(output_dir, exist_ok=True)

    # Run all scenarios
    results = []
    for scenario_name, beta_direct, beta_coop, pi, sigma in scenarios:
        result = run_scenario(scenario_name, beta_direct, beta_coop, pi, sigma, output_dir)
        results.append(result)

    # Save summary
    results_df = pd.DataFrame(results)
    results_df.to_csv(os.path.join(output_dir, 'validation_summary.csv'), index=False)

    # Print summary
    print(f"\n{'='*60}")
    print("VALIDATION SUMMARY")
    print(f"{'='*60}")
    print(results_df[['scenario', 'classification', 'expected_class', 'correct_classification',
                      'beta_coop_error', 'pi_error', 'delta_waic']].to_string(index=False))

    # Overall accuracy
    accuracy = results_df['correct_classification'].mean()
    print(f"\nClassification accuracy: {accuracy:.1%}")

    # Mean absolute errors
    mae_beta_coop = results_df[results_df['expected_class'] == 'cooperative']['beta_coop_error'].mean()
    mae_pi = results_df[results_df['expected_class'] == 'cooperative']['pi_error'].mean()

    print(f"\nParameter recovery (cooperative scenarios only):")
    print(f"  MAE(β_coop): {mae_beta_coop:.2f}")
    print(f"  MAE(π): {mae_pi:.3f}")

    print(f"\nFull results saved to {output_dir}/validation_summary.csv")


if __name__ == '__main__':
    main()
