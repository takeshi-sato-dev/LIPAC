"""
Bayesian Mixture Causal Model for Cooperative Lipid Dynamics

This module implements Bayesian mixture causal models to detect cooperative
lipid effects in MD simulations. When linear models produce sign-inconsistent
effects across replicates, the mixture model decomposes bound-state observations
into latent subpopulations with direct vs. cooperative binding modes.

Author: Takeshi Sato
Kyoto Pharmaceutical University
"""

import argparse
import os
import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import warnings
warnings.filterwarnings('ignore')


def load_stage1_data(input_path):
    """Load Stage 1 output CSV

    Parameters
    ----------
    input_path : str
        Path to Stage 1 CSV file

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: frame, protein_copy, binding_state, lipid_type, contact_count
    """
    df = pd.read_csv(input_path)
    required_cols = ['frame', 'protein_copy', 'binding_state', 'lipid_type', 'contact_count']

    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    return df


def fit_linear_model(Y, T, protein_idx, n_proteins):
    """Fit hierarchical linear model: Y_i = α_p + β_p · T_i + ε_i

    Parameters
    ----------
    Y : np.ndarray
        Contact counts
    T : np.ndarray
        Treatment (binding state: 0=unbound, 1=bound)
    protein_idx : np.ndarray
        Protein copy indices
    n_proteins : int
        Number of protein copies

    Returns
    -------
    pm.Model, arviz.InferenceData
        PyMC model and trace
    """
    with pm.Model() as linear_model:
        # Hyperpriors
        Y_unbound = Y[T == 0]
        mu_alpha = pm.Normal('mu_alpha', mu=np.mean(Y_unbound), sigma=2*np.std(Y))
        sigma_alpha = pm.HalfNormal('sigma_alpha', sigma=np.std(Y))
        mu_beta = pm.Normal('mu_beta', mu=0, sigma=2*np.std(Y))
        sigma_beta = pm.HalfNormal('sigma_beta', sigma=np.std(Y))

        # Per-protein parameters
        alpha_p = pm.Normal('alpha_p', mu=mu_alpha, sigma=sigma_alpha, shape=n_proteins)
        beta_p = pm.Normal('beta_p', mu=mu_beta, sigma=sigma_beta, shape=n_proteins)

        # Observation noise
        sigma = pm.HalfNormal('sigma', sigma=np.std(Y))

        # Linear model
        mu = alpha_p[protein_idx] + beta_p[protein_idx] * T

        # Likelihood
        Y_obs = pm.Normal('Y_obs', mu=mu, sigma=sigma, observed=Y)

        # Sampling
        trace = pm.sample(draws=2000, tune=1000, chains=4, random_seed=42,
                         target_accept=0.9, return_inferencedata=True)

    return linear_model, trace


def fit_mixture_model(Y, T, protein_idx, n_proteins):
    """Fit mixture model for bound state

    Unbound: Y ~ N(α_p, σ)
    Bound: Y ~ π·N(α_p + β_direct,p + β_coop, σ) + (1-π)·N(α_p + β_direct,p, σ)

    Parameters
    ----------
    Y : np.ndarray
        Contact counts
    T : np.ndarray
        Treatment (binding state)
    protein_idx : np.ndarray
        Protein copy indices
    n_proteins : int
        Number of protein copies

    Returns
    -------
    pm.Model, arviz.InferenceData
        PyMC model and trace
    """
    with pm.Model() as mixture_model:
        # Hyperpriors
        Y_unbound = Y[T == 0]
        mu_alpha = pm.Normal('mu_alpha', mu=np.mean(Y_unbound), sigma=2*np.std(Y))
        sigma_alpha = pm.HalfNormal('sigma_alpha', sigma=np.std(Y))
        mu_beta_direct = pm.Normal('mu_beta_direct', mu=0, sigma=2*np.std(Y))
        sigma_beta_direct = pm.HalfNormal('sigma_beta_direct', sigma=np.std(Y))

        # Per-protein parameters
        alpha_p = pm.Normal('alpha_p', mu=mu_alpha, sigma=sigma_alpha, shape=n_proteins)
        beta_direct_p = pm.Normal('beta_direct_p', mu=mu_beta_direct,
                                   sigma=sigma_beta_direct, shape=n_proteins)

        # Shared cooperative effect
        beta_coop = pm.Normal('beta_coop', mu=0, sigma=2*np.std(Y))
        pi = pm.Beta('pi', alpha=2, beta=2)

        # Observation noise
        sigma = pm.HalfNormal('sigma', sigma=np.std(Y))

        # Separate unbound and bound observations
        unbound_mask = T == 0
        bound_mask = T == 1

        Y_unbound_obs = Y[unbound_mask]
        Y_bound_obs = Y[bound_mask]
        protein_idx_unbound = protein_idx[unbound_mask]
        protein_idx_bound = protein_idx[bound_mask]

        # Unbound likelihood
        mu_unbound = alpha_p[protein_idx_unbound]
        Y_unbound_likelihood = pm.Normal('Y_unbound', mu=mu_unbound, sigma=sigma,
                                          observed=Y_unbound_obs)

        # Bound likelihood (mixture)
        # Component 1: direct only
        mu_direct = alpha_p[protein_idx_bound] + beta_direct_p[protein_idx_bound]
        # Component 2: direct + cooperative
        mu_coop = alpha_p[protein_idx_bound] + beta_direct_p[protein_idx_bound] + beta_coop

        # Log-sum-exp for numerical stability
        logp_direct = pm.logp(pm.Normal.dist(mu=mu_direct, sigma=sigma), Y_bound_obs)
        logp_coop = pm.logp(pm.Normal.dist(mu=mu_coop, sigma=sigma), Y_bound_obs)

        # Mixture log-likelihood
        logp_mixture = pm.math.logsumexp(
            pm.math.stack([
                pm.math.log(1 - pi) + logp_direct,
                pm.math.log(pi) + logp_coop
            ], axis=0),
            axis=0
        )

        # Add mixture likelihood as potential
        pm.Potential('Y_bound', logp_mixture)

        # Sampling
        trace = pm.sample(draws=3000, tune=1500, chains=4, random_seed=42,
                         target_accept=0.95, return_inferencedata=True)

    return mixture_model, trace


def compute_waic(model, trace, subsample_size=500):
    """Compute WAIC for model comparison

    Parameters
    ----------
    model : pm.Model
        PyMC model
    trace : arviz.InferenceData
        MCMC trace
    subsample_size : int
        Number of posterior draws to subsample

    Returns
    -------
    float
        WAIC value
    """
    try:
        # Subsample posterior for computational efficiency
        n_draws = len(trace.posterior.draw)
        if n_draws > subsample_size:
            rng = np.random.default_rng(42)
            subsample_idx = rng.choice(n_draws, size=subsample_size, replace=False)
            trace_sub = trace.sel(draw=subsample_idx)
        else:
            trace_sub = trace

        waic = az.waic(trace_sub, pointwise=True)
        return waic.elpd_waic
    except Exception as e:
        print(f"Warning: WAIC computation failed: {e}")
        return np.nan


def check_convergence(trace):
    """Check convergence diagnostics

    Parameters
    ----------
    trace : arviz.InferenceData
        MCMC trace

    Returns
    -------
    dict
        Convergence diagnostics: max_rhat, min_ess
    """
    summary = az.summary(trace)
    max_rhat = summary['r_hat'].max()
    min_ess = summary['ess_bulk'].min()

    return {
        'max_rhat': max_rhat,
        'min_ess': min_ess,
        'converged': max_rhat < 1.01 and min_ess > 400
    }


def classify_lipid(delta_waic, beta_coop_hdi):
    """Classify lipid as linear or cooperative

    Parameters
    ----------
    delta_waic : float
        ΔWAIC = WAIC_linear - WAIC_mixture
    beta_coop_hdi : tuple
        95% HDI of β_coop

    Returns
    -------
    str
        'cooperative' or 'linear'
    """
    if delta_waic > 2 and (beta_coop_hdi[0] > 0 or beta_coop_hdi[1] < 0):
        return 'cooperative'
    else:
        return 'linear'


def analyze_lipid_type(df_lipid, lipid_name, output_dir):
    """Analyze single lipid type with both models

    Parameters
    ----------
    df_lipid : pd.DataFrame
        Data for single lipid type
    lipid_name : str
        Name of lipid type
    output_dir : str
        Output directory

    Returns
    -------
    dict
        Analysis results
    """
    print(f"\n===== Analyzing {lipid_name} =====")

    # Prepare data
    Y = df_lipid['contact_count'].values
    T = df_lipid['binding_state'].values
    protein_labels = df_lipid['protein_copy'].unique()
    protein_map = {p: i for i, p in enumerate(protein_labels)}
    protein_idx = df_lipid['protein_copy'].map(protein_map).values
    n_proteins = len(protein_labels)

    print(f"N observations: {len(Y)}")
    print(f"N proteins: {n_proteins}")
    print(f"Bound fraction: {np.mean(T):.2%}")

    # Fit linear model
    print("\nFitting linear model...")
    linear_model, linear_trace = fit_linear_model(Y, T, protein_idx, n_proteins)
    linear_conv = check_convergence(linear_trace)
    print(f"Linear model: R̂_max={linear_conv['max_rhat']:.3f}, ESS_min={linear_conv['min_ess']:.0f}")

    # Fit mixture model
    print("\nFitting mixture model...")
    mixture_model, mixture_trace = fit_mixture_model(Y, T, protein_idx, n_proteins)
    mixture_conv = check_convergence(mixture_trace)
    print(f"Mixture model: R̂_max={mixture_conv['max_rhat']:.3f}, ESS_min={mixture_conv['min_ess']:.0f}")

    # Compute WAIC
    print("\nComputing WAIC...")
    waic_linear = compute_waic(linear_model, linear_trace)
    waic_mixture = compute_waic(mixture_model, mixture_trace)
    delta_waic = waic_linear - waic_mixture
    print(f"WAIC_linear: {waic_linear:.1f}")
    print(f"WAIC_mixture: {waic_mixture:.1f}")
    print(f"ΔWAIC: {delta_waic:.1f}")

    # Extract parameters
    mixture_summary = az.summary(mixture_trace, hdi_prob=0.95)

    beta_coop_mean = mixture_summary.loc['beta_coop', 'mean']
    beta_coop_hdi = (mixture_summary.loc['beta_coop', 'hdi_2.5%'],
                     mixture_summary.loc['beta_coop', 'hdi_97.5%'])
    pi_mean = mixture_summary.loc['pi', 'mean']
    pi_hdi = (mixture_summary.loc['pi', 'hdi_2.5%'],
              mixture_summary.loc['pi', 'hdi_97.5%'])

    # Beta_direct per protein
    beta_direct_cols = [col for col in mixture_summary.index if col.startswith('beta_direct_p[')]
    beta_direct_means = mixture_summary.loc[beta_direct_cols, 'mean'].values

    # Classification
    classification = classify_lipid(delta_waic, beta_coop_hdi)

    print(f"\nClassification: {classification.upper()}")
    print(f"β_coop: {beta_coop_mean:.2f} [{beta_coop_hdi[0]:.2f}, {beta_coop_hdi[1]:.2f}]")
    print(f"π: {pi_mean:.2f} [{pi_hdi[0]:.2f}, {pi_hdi[1]:.2f}]")
    print(f"β_direct (mean across proteins): {np.mean(beta_direct_means):.2f} ± {np.std(beta_direct_means):.2f}")

    # Save traces
    lipid_dir = os.path.join(output_dir, lipid_name)
    os.makedirs(lipid_dir, exist_ok=True)

    az.to_netcdf(linear_trace, os.path.join(lipid_dir, 'linear_trace.nc'))
    az.to_netcdf(mixture_trace, os.path.join(lipid_dir, 'mixture_trace.nc'))

    linear_summary = az.summary(linear_trace, hdi_prob=0.95)
    linear_summary.to_csv(os.path.join(lipid_dir, 'linear_summary.csv'))
    mixture_summary.to_csv(os.path.join(lipid_dir, 'mixture_summary.csv'))

    # Return results
    return {
        'lipid_name': lipid_name,
        'classification': classification,
        'delta_waic': delta_waic,
        'waic_linear': waic_linear,
        'waic_mixture': waic_mixture,
        'beta_coop_mean': beta_coop_mean,
        'beta_coop_hdi': beta_coop_hdi,
        'pi_mean': pi_mean,
        'pi_hdi': pi_hdi,
        'beta_direct_mean': np.mean(beta_direct_means),
        'beta_direct_std': np.std(beta_direct_means),
        'linear_converged': linear_conv['converged'],
        'mixture_converged': mixture_conv['converged']
    }


def main():
    """Main CLI interface"""
    parser = argparse.ArgumentParser(
        description='Bayesian Mixture Causal Analysis for Lipid Dynamics'
    )
    parser.add_argument('--input', type=str, required=True,
                       help='Path to Stage 1 output CSV')
    parser.add_argument('--output', type=str, default='mixture_output',
                       help='Output directory')
    parser.add_argument('--lipid-types', type=str, nargs='+', default=None,
                       help='Specific lipid types to analyze (default: all)')

    args = parser.parse_args()

    # Load data
    print("Loading data...")
    df = load_stage1_data(args.input)

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Get lipid types to analyze
    if args.lipid_types is None:
        lipid_types = df['lipid_type'].unique()
    else:
        lipid_types = args.lipid_types

    print(f"Analyzing {len(lipid_types)} lipid types: {', '.join(lipid_types)}")

    # Analyze each lipid type
    results = []
    for lipid_name in lipid_types:
        df_lipid = df[df['lipid_type'] == lipid_name].copy()

        if len(df_lipid) < 100:
            print(f"\nSkipping {lipid_name}: insufficient data (N={len(df_lipid)})")
            continue

        try:
            result = analyze_lipid_type(df_lipid, lipid_name, args.output)
            results.append(result)
        except Exception as e:
            print(f"\nError analyzing {lipid_name}: {e}")
            continue

    # Save summary
    results_df = pd.DataFrame(results)
    results_df.to_csv(os.path.join(args.output, 'results_summary.csv'), index=False)

    print("\n===== Summary =====")
    print(results_df[['lipid_name', 'classification', 'delta_waic',
                      'beta_coop_mean', 'pi_mean']].to_string(index=False))

    print(f"\nResults saved to {args.output}")


if __name__ == '__main__':
    main()
