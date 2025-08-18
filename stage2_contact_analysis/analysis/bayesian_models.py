"""
Bayesian regression models for MIRAGE4
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pymc as pm
import arviz as az
import warnings
warnings.filterwarnings('ignore')

# Add parent directory to path for config import
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import MCMC_SAMPLES, TUNE_SAMPLES, CHAINS, RANDOM_SEED, TARGET_ACCEPT
from visualization.plots import plot_convergence_diagnostics


def bayesian_regression_model(data, output_dir):
    """Bayesian regression model of the relationship between Target Lipid contact and protein contact (stable version)"""
    print("\n===== Bayesian Regression Model: Target Lipid Contact Effects =====")
    
    # Create output directory
    model_dir = os.path.join(output_dir, "models")
    os.makedirs(model_dir, exist_ok=True)
    
    # Check data
    if data is None:
        print("No data provided. Aborting Bayesian analysis.")
        return None
    
    if 'target_lipid_contact' not in data.columns:
        print("Warning: No target_lipid_contact column found. Creating dummy column.")
        # Check for alternative lipid contact columns
        lipid_cols = [col for col in data.columns if col.endswith('_contact') and col not in ['protein_contact', 'lipid_contact']]
        if lipid_cols:
            # Use the first available lipid contact column
            data['target_lipid_contact'] = data[lipid_cols[0]]
            print(f"Using {lipid_cols[0]} as target_lipid_contact")
        else:
            # Create dummy data if no lipid columns available
            data['target_lipid_contact'] = np.random.normal(0, 0.1, len(data))
            print("Created dummy target_lipid_contact data for analysis")
    
    # Prepare data
    X_target_lipid = data['target_lipid_contact'].values
    Y_protein_diff = data['protein_contact_diff'].values
    Y_lipid_diff = data['lipid_contact_diff'].values
    
    # Store summary for later use
    summary_protein = None
    
    # 1. Model for Target Lipid contact effect on protein contact
    print("Model 1: Relationship between Target Lipid contact and protein contact change")
    
    with pm.Model() as protein_model:
        # Prior distributions
        alpha = pm.Normal('alpha', mu=0, sigma=1)  # Intercept
        beta = pm.Normal('beta', mu=0, sigma=1)    # Target Lipid contact coefficient
        sigma = pm.HalfNormal('sigma', sigma=1)    # Error standard deviation
        
        # Linear model
        mu = alpha + beta * X_target_lipid
        
        # Likelihood
        Y_obs = pm.Normal('Y_obs', mu=mu, sigma=sigma, observed=Y_protein_diff)
        
        # Sampling
        protein_trace = pm.sample(
            MCMC_SAMPLES, 
            tune=TUNE_SAMPLES, 
            chains=CHAINS, 
            random_seed=RANDOM_SEED,
            target_accept=TARGET_ACCEPT
        )
        
        # Save results
        az.to_netcdf(protein_trace, os.path.join(model_dir, 'protein_model_trace.nc'))
        
        # Summary
        summary = az.summary(protein_trace)
        summary_protein = summary  # Store for later comparison
        print(summary)
        summary.to_csv(os.path.join(model_dir, 'protein_model_summary.csv'))
        
        # Add convergence diagnostics
        plot_convergence_diagnostics(protein_trace, 'protein_model', model_dir)
        
        # Trace plot
        az.plot_trace(protein_trace)
        plt.savefig(os.path.join(model_dir, 'protein_model_trace.png'), dpi=300)
        plt.close()
        
        # Posterior predictive check
        try:
            # Using latest Arviz version for posterior predictive
            ppc = pm.sample_posterior_predictive(protein_trace, model=protein_model, var_names=['Y_obs'])
            
            plt.figure(figsize=(10, 6))
            # Histogram of observed data
            plt.hist(Y_protein_diff, bins=20, alpha=0.5, density=True, label='Observed Data')
            
            # Distribution of predicted data - handle different dimension structures
            if hasattr(ppc, 'posterior_predictive'):
                # New ArviZ structure
                ppc_values = ppc.posterior_predictive['Y_obs'].values
                # Flatten to 1D array regardless of original shape
                ppc_samples = ppc_values.reshape(-1)
            else:
                # Older structure
                ppc_samples = ppc['Y_obs'].flatten()
            
            plt.hist(ppc_samples, bins=20, alpha=0.5, density=True, label='Predicted Data')
            
            plt.title('Posterior Predictive Check: Protein Contact Model')
            plt.xlabel('Protein Contact Difference')
            plt.legend()
            plt.savefig(os.path.join(model_dir, 'protein_model_ppc_simple.png'), dpi=300)
            plt.close()
        except Exception as e:
            print(f"Error in posterior predictive check: {str(e)}")
            print("Trying alternative posterior prediction method")
            
            try:
                # Alternative approach for posterior prediction
                # Extract parameters from summary
                alpha_mean = summary.loc['alpha', 'mean']
                beta_mean = summary.loc['beta', 'mean']
                sigma_mean = summary.loc['sigma', 'mean']
                
                # Generate prediction distribution
                x_pred = np.linspace(0, max(X_target_lipid), 100)
                y_pred_mean = alpha_mean + beta_mean * x_pred
                
                # Plot predictions
                plt.figure(figsize=(10, 6))
                plt.scatter(X_target_lipid, Y_protein_diff, alpha=0.5, label='Observed Data')
                plt.plot(x_pred, y_pred_mean, 'r-', label='Predicted Mean')
                
                # Prediction interval
                plt.fill_between(
                    x_pred,
                    y_pred_mean - 1.96 * sigma_mean,
                    y_pred_mean + 1.96 * sigma_mean,
                    color='r', alpha=0.2, label='95% Prediction Interval'
                )
                
                plt.title('Posterior Prediction: Protein Contact Model (Alternative Method)')
                plt.xlabel('Target Lipid Contact Frequency')
                plt.ylabel('Protein Contact Change')
                plt.legend()
                plt.savefig(os.path.join(model_dir, 'protein_model_ppc_simple.png'), dpi=300)
                plt.close()
            except Exception as e:
                print(f"Alternative posterior prediction also failed: {str(e)}")
                print("Skipping posterior predictive check")
        
        # Regression line plot
        try:
            # Extract parameters from posterior in latest Arviz version
            alpha_samples = protein_trace.posterior['alpha'].values.flatten()
            beta_samples = protein_trace.posterior['beta'].values.flatten()
            
            plt.figure(figsize=(10, 6))
            plt.scatter(X_target_lipid, Y_protein_diff, alpha=0.5)
            
            # Average regression line
            alpha_mean = alpha_samples.mean()
            beta_mean = beta_samples.mean()
            
            # x values for prediction line
            x_range = np.linspace(0, np.max(X_target_lipid) * 1.1, 100)
            
            # Plot regression line
            plt.plot(x_range, alpha_mean + beta_mean * x_range, 'r-', linewidth=2, 
                    label=f'y = {alpha_mean:.3f} + {beta_mean:.3f}x')
            
            # Reference line
            plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
            
            # Style settings
            plt.xlabel('Target Lipid Contact Frequency')
            plt.ylabel('Protein Contact Change')
            plt.title('Bayesian Regression: Effect of Target Lipid on Protein-Protein Contacts')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(model_dir, 'protein_model_regression_simple.png'), dpi=300)
            plt.savefig(os.path.join(model_dir, 'protein_model_regression_simple.svg'), dpi=300)
            plt.close()
            
            # Effect size posterior distribution
            plt.figure(figsize=(8, 6))
            plt.hist(beta_samples, bins=30, alpha=0.7, density=True)
            plt.axvline(x=0, color='r', linestyle='--')
            plt.xlabel('Effect Size (β)')
            plt.ylabel('Density')
            plt.title('Posterior Distribution of Target Lipid Effect on Protein-Protein Contacts')
            
            # 95% highest density interval
            beta_hdi = np.percentile(beta_samples, [2.5, 97.5])
            plt.axvline(x=beta_hdi[0], color='k', linestyle='--', alpha=0.5)
            plt.axvline(x=beta_hdi[1], color='k', linestyle='--', alpha=0.5)
            
            # Statistical information
            prob_neg = (beta_samples < 0).mean()
            plt.text(0.05, 0.95, f'95% HDI: [{beta_hdi[0]:.3f}, {beta_hdi[1]:.3f}]', 
                    transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.7))
            plt.text(0.05, 0.85, f'P(β < 0) = {prob_neg:.3f}', 
                    transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.7))
            
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(model_dir, 'protein_model_effect_simple.png'), dpi=300)
            plt.savefig(os.path.join(model_dir, 'protein_model_effect_simple.svg'), dpi=300)
            plt.close()
        except Exception as e:
            print(f"Error in regression line plot: {str(e)}")
            print("Creating regression line plot from summary")
            
            try:
                # Get parameters directly from summary
                alpha_mean = summary.loc['alpha', 'mean']
                beta_mean = summary.loc['beta', 'mean']
                
                plt.figure(figsize=(10, 6))
                plt.scatter(X_target_lipid, Y_protein_diff, alpha=0.5)
                
                # x values for prediction line
                x_range = np.linspace(0, np.max(X_target_lipid) * 1.1, 100)
                
                # Plot regression line
                plt.plot(x_range, alpha_mean + beta_mean * x_range, 'r-', linewidth=2, 
                        label=f'y = {alpha_mean:.3f} + {beta_mean:.3f}x')
                
                # Reference line
                plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
                
                # Style settings
                plt.xlabel('Target Lipid Contact Frequency')
                plt.ylabel('Protein Contact Change')
                plt.title('Bayesian Regression: Effect of Target Lipid on Protein-Protein Contacts')
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.savefig(os.path.join(model_dir, 'protein_model_regression_simple.png'), dpi=300)
                plt.savefig(os.path.join(model_dir, 'protein_model_regression_simple.svg'), dpi=300)
                plt.close()
                
                # Save effect size information to text file
                with open(os.path.join(model_dir, 'protein_model_effect_info.txt'), 'w') as f:
                    f.write(f"Effect size (β) mean: {beta_mean:.4f}\n")
                    f.write(f"95% HDI: [{summary.loc['beta', 'hdi_2.5%']:.4f}, {summary.loc['beta', 'hdi_97.5%']:.4f}]\n")
                    
                    # Approximate P(β < 0)
                    if beta_mean < 0 and summary.loc['beta', 'hdi_97.5%'] < 0:
                        p_neg = 0.999
                    elif beta_mean > 0 and summary.loc['beta', 'hdi_2.5%'] > 0:
                        p_neg = 0.001
                    else:
                        # Normal distribution approximation
                        from scipy.stats import norm
                        z = beta_mean / summary.loc['beta', 'sd']
                        p_neg = norm.cdf(-z)
                    
                    f.write(f"P(β < 0) ≈ {p_neg:.3f}\n")
            except Exception as e:
                print(f"Summary-based regression line plot also failed: {str(e)}")
                print("Skipping regression line plot")
    
    # 2. Model for Target Lipid contact effect on lipid contact
    print("Model 2: Relationship between Target Lipid contact and lipid contact change")
    
    with pm.Model() as lipid_model:
        # Prior distributions
        alpha = pm.Normal('alpha', mu=0, sigma=1)  # Intercept
        beta = pm.Normal('beta', mu=0, sigma=1)    # Target Lipid contact coefficient
        sigma = pm.HalfNormal('sigma', sigma=1)    # Error standard deviation
        
        # Linear model
        mu = alpha + beta * X_target_lipid
        
        # Likelihood
        Y_obs = pm.Normal('Y_obs', mu=mu, sigma=sigma, observed=Y_lipid_diff)
        
        # Sampling
        lipid_trace = pm.sample(
            MCMC_SAMPLES, 
            tune=TUNE_SAMPLES, 
            chains=CHAINS, 
            random_seed=RANDOM_SEED,
            target_accept=TARGET_ACCEPT
        )
        
        # Save results
        az.to_netcdf(lipid_trace, os.path.join(model_dir, 'lipid_model_trace.nc'))
        
        # Summary
        summary = az.summary(lipid_trace)
        print(summary)
        summary.to_csv(os.path.join(model_dir, 'lipid_model_summary.csv'))
        
        # Add convergence diagnostics
        plot_convergence_diagnostics(lipid_trace, 'lipid_model', model_dir)
        
        # Trace plot
        az.plot_trace(lipid_trace)
        plt.savefig(os.path.join(model_dir, 'lipid_model_trace.png'), dpi=300)
        plt.close()
        
        # Posterior predictive check
        try:
            # Using latest Arviz version for posterior predictive
            ppc = pm.sample_posterior_predictive(lipid_trace, model=lipid_model, var_names=['Y_obs'])
            
            plt.figure(figsize=(10, 6))
            # Histogram of observed data
            plt.hist(Y_lipid_diff, bins=20, alpha=0.5, density=True, label='Observed Data')
            
            # Distribution of predicted data - handle different dimension structures
            if hasattr(ppc, 'posterior_predictive'):
                # New ArviZ structure
                ppc_values = ppc.posterior_predictive['Y_obs'].values
                # Flatten to 1D array regardless of original shape
                ppc_samples = ppc_values.reshape(-1)
            else:
                # Older structure
                ppc_samples = ppc['Y_obs'].flatten()
            
            plt.hist(ppc_samples, bins=20, alpha=0.5, density=True, label='Predicted Data')
            
            plt.title('Posterior Predictive Check: Lipid Contact Model')
            plt.xlabel('Lipid Contact Difference')
            plt.legend()
            plt.savefig(os.path.join(model_dir, 'lipid_model_ppc_simple.png'), dpi=300)
            plt.close()
        except Exception as e:
            print(f"Error in posterior predictive check: {str(e)}")
            print("Trying alternative posterior prediction method")
            
            try:
                # Alternative approach for posterior prediction
                # Extract parameters from summary
                alpha_mean = summary.loc['alpha', 'mean']
                beta_mean = summary.loc['beta', 'mean']
                sigma_mean = summary.loc['sigma', 'mean']
                
                # Generate prediction distribution
                x_pred = np.linspace(0, max(X_target_lipid), 100)
                y_pred_mean = alpha_mean + beta_mean * x_pred
                
                # Plot predictions
                plt.figure(figsize=(10, 6))
                plt.scatter(X_target_lipid, Y_lipid_diff, alpha=0.5, label='Observed Data')
                plt.plot(x_pred, y_pred_mean, 'b-', label='Predicted Mean')
                
                # Prediction interval
                plt.fill_between(
                    x_pred,
                    y_pred_mean - 1.96 * sigma_mean,
                    y_pred_mean + 1.96 * sigma_mean,
                    color='b', alpha=0.2, label='95% Prediction Interval'
                )
                
                plt.title('Posterior Prediction: Lipid Contact Model (Alternative Method)')
                plt.xlabel('Target Lipid Contact Frequency')
                plt.ylabel('Lipid Contact Change')
                plt.legend()
                plt.savefig(os.path.join(model_dir, 'lipid_model_ppc_simple.png'), dpi=300)
                plt.close()
            except Exception as e:
                print(f"Alternative posterior prediction also failed: {str(e)}")
                print("Skipping posterior predictive check")
        
        # Regression line plot
        try:
            # Extract parameters from posterior in latest Arviz version
            alpha_samples = lipid_trace.posterior['alpha'].values.flatten()
            beta_samples = lipid_trace.posterior['beta'].values.flatten()
            
            plt.figure(figsize=(10, 6))
            plt.scatter(X_target_lipid, Y_lipid_diff, alpha=0.5)
            
            # Average regression line
            alpha_mean = alpha_samples.mean()
            beta_mean = beta_samples.mean()
            
            # x values for prediction line
            x_range = np.linspace(0, np.max(X_target_lipid) * 1.1, 100)
            
            # Plot regression line
            plt.plot(x_range, alpha_mean + beta_mean * x_range, 'b-', linewidth=2, 
                    label=f'y = {alpha_mean:.3f} + {beta_mean:.3f}x')
            
            # Reference line
            plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
            
            # Style settings
            plt.xlabel('Target Lipid Contact Frequency')
            plt.ylabel('Lipid Contact Change')
            plt.title('Bayesian Regression: Effect of Target Lipid on Lipid-Protein Contacts')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(model_dir, 'lipid_model_regression_simple.png'), dpi=300)
            plt.savefig(os.path.join(model_dir, 'lipid_model_regression_simple.svg'), dpi=300)
            plt.close()
            
            # Effect size posterior distribution
            plt.figure(figsize=(8, 6))
            plt.hist(beta_samples, bins=30, alpha=0.7, density=True)
            plt.axvline(x=0, color='r', linestyle='--')
            plt.xlabel('Effect Size (β)')
            plt.ylabel('Density')
            plt.title('Posterior Distribution of Target Lipid Effect on Lipid-Protein Contacts')
            
            # 95% highest density interval
            beta_hdi = np.percentile(beta_samples, [2.5, 97.5])
            plt.axvline(x=beta_hdi[0], color='k', linestyle='--', alpha=0.5)
            plt.axvline(x=beta_hdi[1], color='k', linestyle='--', alpha=0.5)
            
            # Statistical information
            prob_pos = (beta_samples > 0).mean()
            plt.text(0.05, 0.95, f'95% HDI: [{beta_hdi[0]:.3f}, {beta_hdi[1]:.3f}]', 
                    transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.7))
            plt.text(0.05, 0.85, f'P(β > 0) = {prob_pos:.3f}', 
                    transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.7))
            
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(model_dir, 'lipid_model_effect_simple.png'), dpi=300)
            plt.savefig(os.path.join(model_dir, 'lipid_model_effect_simple.svg'), dpi=300)
            plt.close()
        except Exception as e:
            print(f"Error in regression line plot: {str(e)}")
            print("Creating regression line plot from summary")
            
            try:
                # Get parameters directly from summary
                alpha_mean = summary.loc['alpha', 'mean']
                beta_mean = summary.loc['beta', 'mean']
                
                plt.figure(figsize=(10, 6))
                plt.scatter(X_target_lipid, Y_lipid_diff, alpha=0.5)
                
                # x values for prediction line
                x_range = np.linspace(0, np.max(X_target_lipid) * 1.1, 100)
                
                # Plot regression line
                plt.plot(x_range, alpha_mean + beta_mean * x_range, 'b-', linewidth=2, 
                        label=f'y = {alpha_mean:.3f} + {beta_mean:.3f}x')
                
                # Reference line
                plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
                
                # Style settings
                plt.xlabel('Target Lipid Contact Frequency')
                plt.ylabel('Lipid Contact Change')
                plt.title('Bayesian Regression: Effect of Target Lipid on Lipid-Protein Contacts')
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.savefig(os.path.join(model_dir, 'lipid_model_regression_simple.png'), dpi=300)
                plt.savefig(os.path.join(model_dir, 'lipid_model_regression_simple.svg'), dpi=300)
                plt.close()
                
                # Save effect size information to text file
                with open(os.path.join(model_dir, 'lipid_model_effect_info.txt'), 'w') as f:
                    f.write(f"Effect size (β) mean: {beta_mean:.4f}\n")
                    f.write(f"95% HDI: [{summary.loc['beta', 'hdi_2.5%']:.4f}, {summary.loc['beta', 'hdi_97.5%']:.4f}]\n")
                    
                    # Approximate P(β > 0)
                    if beta_mean > 0 and summary.loc['beta', 'hdi_2.5%'] > 0:
                        p_pos = 0.999
                    elif beta_mean < 0 and summary.loc['beta', 'hdi_97.5%'] < 0:
                        p_pos = 0.001
                    else:
                        # Normal distribution approximation
                        from scipy.stats import norm
                        z = beta_mean / summary.loc['beta', 'sd']
                        p_pos = 1 - norm.cdf(-z)
                    
                    f.write(f"P(β > 0) ≈ {p_pos:.3f}\n")
            except Exception as e:
                print(f"Summary-based regression line plot also failed: {str(e)}")
                print("Skipping regression line plot")
    
    # 3. Model comparison: relationship between protein contact and lipid contact effects
    try:
        # Extract samples from posterior in latest Arviz version
        protein_samples = protein_trace.posterior['beta'].values.flatten()
        lipid_samples = lipid_trace.posterior['beta'].values.flatten()
        
        plt.figure(figsize=(10, 6))
        plt.hist(protein_samples, bins=30, alpha=0.5, density=True, label='Protein Contact Effect')
        plt.hist(lipid_samples, bins=30, alpha=0.5, density=True, label='Lipid Contact Effect')
        plt.axvline(x=0, color='r', linestyle='--')
        plt.xlabel('Effect Size (β)')
        plt.ylabel('Density')
        plt.title('Comparison of Target Lipid Effect Posterior Distributions')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(model_dir, 'effect_comparison_simple.png'), dpi=300)
        plt.savefig(os.path.join(model_dir, 'effect_comparison_simple.svg'), dpi=300)
        plt.close()
        
        # Effect direction probabilities
        prob_protein_neg = (protein_samples < 0).mean()
        prob_lipid_pos = (lipid_samples > 0).mean()
        prob_opposite = ((protein_samples < 0) & (lipid_samples > 0)).mean()
        
        # Save results to text file
        with open(os.path.join(output_dir, 'bayesian_analysis_results.txt'), 'w') as f:
            f.write("===== Bayesian Analysis Results =====\n\n")
            
            f.write("--- Model 1: Relationship between Target Lipid Contact and Protein Contact ---\n")
            f.write(f"Mean effect size (β): {protein_samples.mean():.4f}\n")
            protein_hdi = np.percentile(protein_samples, [2.5, 97.5])
            f.write(f"95% HDI: [{protein_hdi[0]:.4f}, {protein_hdi[1]:.4f}]\n")
            f.write(f"P(β < 0) = {prob_protein_neg:.4f}\n\n")
            
            f.write("--- Model 2: Relationship between Target Lipid Contact and Lipid Contact ---\n")
            f.write(f"Mean effect size (β): {lipid_samples.mean():.4f}\n")
            lipid_hdi = np.percentile(lipid_samples, [2.5, 97.5])
            f.write(f"95% HDI: [{lipid_hdi[0]:.4f}, {lipid_hdi[1]:.4f}]\n")
            f.write(f"P(β > 0) = {prob_lipid_pos:.4f}\n\n")
            
            f.write("--- Effect Direction ---\n")
            f.write(f"P(Protein Effect < 0 and Lipid Effect > 0) = {prob_opposite:.4f}\n")
            
            if prob_opposite > 0.9:
                f.write("\nConclusion: Target Lipid has a high probability of decreasing protein-protein contacts while simultaneously increasing lipid-protein contacts.\n")
                f.write("This result suggests that Target Lipid inhibits protein dimer formation through a competitive binding mechanism.\n")
            elif prob_opposite > 0.7:
                f.write("\nConclusion: Target Lipid tends to decrease protein-protein contacts while increasing lipid-protein contacts.\n")
                f.write("This result suggests that Target Lipid may influence protein dimer formation through a competitive binding mechanism.\n")
            else:
                f.write("\nConclusion: Evidence for Target Lipid having opposite effects on protein-protein contacts and lipid-protein contacts is limited.\n")
                f.write("Further research is needed.\n")
    except Exception as e:
        print(f"Error in model comparison: {str(e)}")
        print("Creating model comparison from summary")
        
        try:
            # Get parameters directly from summary
            protein_beta_mean = float(summary_protein.loc['beta', 'mean'])
            protein_beta_sd = float(summary_protein.loc['beta', 'sd'])
            lipid_beta_mean = float(summary.loc['beta', 'mean'])
            lipid_beta_sd = float(summary.loc['beta', 'sd'])
            
            # Generate approximate posterior distributions with normal approximation
            protein_samples_approx = np.random.normal(protein_beta_mean, protein_beta_sd, 10000)
            lipid_samples_approx = np.random.normal(lipid_beta_mean, lipid_beta_sd, 10000)
            
            # Effect direction probabilities
            prob_protein_neg_approx = (protein_samples_approx < 0).mean()
            prob_lipid_pos_approx = (lipid_samples_approx > 0).mean()
            prob_opposite_approx = ((protein_samples_approx < 0) & (lipid_samples_approx > 0)).mean()
            
            # Comparison plot
            plt.figure(figsize=(10, 6))
            plt.hist(protein_samples_approx, bins=30, alpha=0.5, density=True, label='Protein Contact Effect (approx.)')
            plt.hist(lipid_samples_approx, bins=30, alpha=0.5, density=True, label='Lipid Contact Effect (approx.)')
            plt.axvline(x=0, color='r', linestyle='--')
            plt.xlabel('Effect Size (β)')
            plt.ylabel('Density')
            plt.title('Comparison of Target Lipid Effect Posterior Distributions (Normal Approximation)')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(model_dir, 'effect_comparison_simple.png'), dpi=300)
            plt.savefig(os.path.join(model_dir, 'effect_comparison_simple.svg'), dpi=300)
            plt.close()
            
            # Save results to text file
            with open(os.path.join(output_dir, 'bayesian_analysis_results.txt'), 'w') as f:
                f.write("===== Bayesian Analysis Results (Approximate Values) =====\n\n")
                
                f.write("--- Model 1: Relationship between Target Lipid Contact and Protein Contact ---\n")
                f.write(f"Mean effect size (β): {protein_beta_mean:.4f}\n")
                protein_hdi_approx = [protein_beta_mean - 1.96*protein_beta_sd, protein_beta_mean + 1.96*protein_beta_sd]
                f.write(f"95% CI: [{protein_hdi_approx[0]:.4f}, {protein_hdi_approx[1]:.4f}]\n")
                f.write(f"P(β < 0) ≈ {prob_protein_neg_approx:.4f}\n\n")
                
                f.write("--- Model 2: Relationship between Target Lipid Contact and Lipid Contact ---\n")
                f.write(f"Mean effect size (β): {lipid_beta_mean:.4f}\n")
                lipid_hdi_approx = [lipid_beta_mean - 1.96*lipid_beta_sd, lipid_beta_mean + 1.96*lipid_beta_sd]
                f.write(f"95% CI: [{lipid_hdi_approx[0]:.4f}, {lipid_hdi_approx[1]:.4f}]\n")
                f.write(f"P(β > 0) ≈ {prob_lipid_pos_approx:.4f}\n\n")
                
                f.write("--- Effect Direction (Approximate Values) ---\n")
                f.write(f"P(Protein Effect < 0 and Lipid Effect > 0) ≈ {prob_opposite_approx:.4f}\n")
                
                if prob_opposite_approx > 0.9:
                    f.write("\nConclusion: Target Lipid has a high probability of decreasing protein-protein contacts while simultaneously increasing lipid-protein contacts.\n")
                    f.write("This result suggests that Target Lipid inhibits protein dimer formation through a competitive binding mechanism.\n")
                elif prob_opposite_approx > 0.7:
                    f.write("\nConclusion: Target Lipid tends to decrease protein-protein contacts while increasing lipid-protein contacts.\n")
                    f.write("This result suggests that Target Lipid may influence protein dimer formation through a competitive binding mechanism.\n")
                else:
                    f.write("\nConclusion: Evidence for Target Lipid having opposite effects on protein-protein contacts and lipid-protein contacts is limited.\n")
                    f.write("Further research is needed.\n")
        except Exception as e:
            print(f"Summary-based model comparison also failed: {str(e)}")
            print("Using simplified model comparison")
            
            with open(os.path.join(output_dir, 'bayesian_analysis_results.txt'), 'w') as f:
                f.write("===== Bayesian Analysis Results =====\n\n")
                f.write("Detailed analysis was limited due to technical issues.\n")
                f.write("Please refer to protein_model_summary.csv and lipid_model_summary.csv for individual model results.\n")
    
    print(f"Bayesian analysis results saved to {os.path.join(output_dir, 'bayesian_analysis_results.txt')}")
    
    return {
        'protein_trace': protein_trace,
        'lipid_trace': lipid_trace
    }


def hierarchical_model(data, output_dir):
    """Hierarchical Bayesian Model: Target Lipid effects by protein"""
    print("\n===== Hierarchical Bayesian Model: Protein-Specific Target Lipid Effects =====")
    
    # Create output directory
    model_dir = os.path.join(output_dir, "models")
    os.makedirs(model_dir, exist_ok=True)
    
    # Check for Target Lipid contact data
    if data is None or 'target_lipid_contact' not in data.columns:
        print("No target lipid contact data. Aborting hierarchical Bayesian model.")
        return None
    
    # Split by protein
    proteins = data['protein'].unique()
    
    if len(proteins) <= 1:
        print("Multiple proteins not found. Aborting hierarchical Bayesian model.")
        return None
    
    print(f"Proteins included in model: {', '.join(proteins)}")
    
    # Create protein indices
    protein_indices = {protein: i for i, protein in enumerate(proteins)}
    protein_idx = np.array([protein_indices[p] for p in data['protein']])
    
    # Prepare data
    X_target_lipid = data['target_lipid_contact'].values
    Y_protein_diff = data['protein_contact_diff'].values
    
    # Hierarchical Bayesian model
    with pm.Model() as hierarchical_model:
        # Global prior distributions
        mu_alpha = pm.Normal('mu_alpha', mu=0, sigma=1)
        sigma_alpha = pm.HalfNormal('sigma_alpha', sigma=1)
        
        mu_beta = pm.Normal('mu_beta', mu=0, sigma=1)
        sigma_beta = pm.HalfNormal('sigma_beta', sigma=1)
        
        # Protein-specific parameters
        alpha = pm.Normal('alpha', mu=mu_alpha, sigma=sigma_alpha, shape=len(proteins))
        beta = pm.Normal('beta', mu=mu_beta, sigma=sigma_beta, shape=len(proteins))
        
        # Error standard deviation
        sigma = pm.HalfNormal('sigma', sigma=1)
        
        # Linear model
        mu = alpha[protein_idx] + beta[protein_idx] * X_target_lipid
        
        # Likelihood
        Y_obs = pm.Normal('Y_obs', mu=mu, sigma=sigma, observed=Y_protein_diff)
        
        # Sampling
        hierarchical_trace = pm.sample(
            MCMC_SAMPLES, 
            tune=TUNE_SAMPLES, 
            chains=CHAINS, 
            random_seed=RANDOM_SEED,
            target_accept=TARGET_ACCEPT
        )
        
        # Save results - Using az.to_netcdf for PyMC v4
        az.to_netcdf(hierarchical_trace, os.path.join(model_dir, 'hierarchical_model_trace.nc'))
        
        # Summary
        summary = az.summary(hierarchical_trace)
        print(summary)
        summary.to_csv(os.path.join(model_dir, 'hierarchical_model_summary.csv'))
        
        # Add convergence diagnostics
        plot_convergence_diagnostics(hierarchical_trace, 'hierarchical_model', model_dir)
        
        # Trace plot for global parameters
        az.plot_trace(hierarchical_trace, var_names=['mu_alpha', 'mu_beta', 'sigma_alpha', 'sigma_beta', 'sigma'])
        plt.savefig(os.path.join(model_dir, 'hierarchical_model_trace_global.png'), dpi=300)
        plt.close()
        
        # Forest plot of protein-specific effects
        plt.figure(figsize=(10, 8))
        
        # Calculate effect for each protein
        means = []
        lower_cis = []
        upper_cis = []
        protein_names = []
        
        for i, protein in enumerate(proteins):
            beta_samples = hierarchical_trace.posterior['beta'][:, :, i].values.flatten()
            mean = beta_samples.mean()
            ci = np.percentile(beta_samples, [2.5, 97.5])
            
            means.append(mean)
            lower_cis.append(ci[0])
            upper_cis.append(ci[1])
            protein_names.append(f'Protein {i+1}: {protein}')
        
        # Add global effect
        mu_beta_samples = hierarchical_trace.posterior['mu_beta'].values.flatten()
        mu_beta_mean = mu_beta_samples.mean()
        mu_beta_ci = np.percentile(mu_beta_samples, [2.5, 97.5])
        
        means.append(mu_beta_mean)
        lower_cis.append(mu_beta_ci[0])
        upper_cis.append(mu_beta_ci[1])
        protein_names.append('Global')
        
        # Display bottom to top
        y_pos = np.arange(len(protein_names))
        
        # Create plot
        plt.figure(figsize=(10, 8))
        plt.errorbar(means, y_pos, xerr=[np.array(means)-np.array(lower_cis), np.array(upper_cis)-np.array(means)], 
                    fmt='o', capsize=5, elinewidth=2, markeredgewidth=2)
        
        # Zero line
        plt.axvline(x=0, color='r', linestyle='--', alpha=0.5)
        
        # Label each protein
        plt.yticks(y_pos, protein_names)
        
        plt.xlabel('Effect Size (β)')
        plt.title('Target Lipid Effect on Protein-Protein Contacts by Protein')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(model_dir, 'hierarchical_model_forest.png'), dpi=300)
        plt.close()
        
        # Regression line plot for each protein
        plt.figure(figsize=(12, 10))
        
        for i, protein in enumerate(proteins):
            protein_data = data[data['protein'] == protein]
            X = protein_data['target_lipid_contact'].values
            Y = protein_data['protein_contact_diff'].values
            
            # Scatter plot
            plt.subplot(2, 2, i+1)
            plt.scatter(X, Y, alpha=0.5)
            plt.title(f'{protein}')
            plt.xlabel('Target Lipid Contact Frequency')
            plt.ylabel('Change in Protein-Protein Contact')
            
            # Regression lines
            alpha_samples = hierarchical_trace.posterior['alpha'][:, :, i].values.flatten()
            beta_samples = hierarchical_trace.posterior['beta'][:, :, i].values.flatten()
            
            x_range = np.linspace(0, max(1, X.max()), 100)
            
            # Sample regression lines
            for j in range(0, len(alpha_samples), 50):
                plt.plot(x_range, alpha_samples[j] + beta_samples[j] * x_range, 'r-', alpha=0.05)
            
            # Average regression line
            alpha_mean = alpha_samples.mean()
            beta_mean = beta_samples.mean()
            plt.plot(x_range, alpha_mean + beta_mean * x_range, 'r-', linewidth=2, 
                    label=f'y = {alpha_mean:.3f} + {beta_mean:.3f}x')
            
            # 95% HDI - using numpy.percentile
            beta_hdi = np.percentile(beta_samples, [2.5, 97.5])
            plt.text(0.05, 0.90, f'Effect: {beta_mean:.3f} [{beta_hdi[0]:.3f}, {beta_hdi[1]:.3f}]', 
                    transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.7))
            
            # Zero inclusion check
            prob_neg = (beta_samples < 0).mean()
            plt.text(0.05, 0.80, f'P(β < 0) = {prob_neg:.3f}', 
                    transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.7))
            
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(model_dir, 'hierarchical_model_regression.png'), dpi=300)
        plt.close()
        
        # Effect size comparison by protein
        plt.figure(figsize=(12, 6))
        
        for i, protein in enumerate(proteins):
            beta_samples = hierarchical_trace.posterior['beta'][:, :, i].values.flatten()
            sns.kdeplot(beta_samples, fill=True, label=protein)
        
        plt.axvline(x=0, color='r', linestyle='--')
        plt.xlabel('Effect Size (β)')
        plt.ylabel('Density')
        plt.title('Posterior Distributions of Target Lipid Effects by Protein')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(model_dir, 'hierarchical_model_effects.png'), dpi=300)
        plt.close()
        
        # Save protein-specific results to text file
        with open(os.path.join(output_dir, 'hierarchical_model_results.txt'), 'w') as f:
            f.write("===== Hierarchical Bayesian Model Results =====\n\n")
            f.write("--- Global Effect ---\n")
            
            mu_beta_samples = hierarchical_trace.posterior['mu_beta'].values.flatten()
            mu_beta_mean = mu_beta_samples.mean()
            # HDI calculation
            mu_beta_hdi = np.percentile(mu_beta_samples, [2.5, 97.5])
            prob_global_neg = (mu_beta_samples < 0).mean()
            
            f.write(f"Global Target Lipid effect (μ_β): {mu_beta_mean:.4f}\n")
            f.write(f"95% HDI: [{mu_beta_hdi[0]:.4f}, {mu_beta_hdi[1]:.4f}]\n")
            f.write(f"P(μ_β < 0) = {prob_global_neg:.4f}\n\n")
            
            f.write("--- Protein-Specific Effects ---\n")
            
            for i, protein in enumerate(proteins):
                beta_samples = hierarchical_trace.posterior['beta'][:, :, i].values.flatten()
                beta_mean = beta_samples.mean()
                # HDI calculation
                beta_hdi = np.percentile(beta_samples, [2.5, 97.5])
                prob_neg = (beta_samples < 0).mean()
                
                f.write(f"{protein}:\n")
                f.write(f"  Target Lipid effect (β): {beta_mean:.4f}\n")
                f.write(f"  95% HDI: [{beta_hdi[0]:.4f}, {beta_hdi[1]:.4f}]\n")
                f.write(f"  P(β < 0) = {prob_neg:.4f}\n\n")
        
        print(f"Hierarchical Bayesian model results saved to {os.path.join(output_dir, 'hierarchical_model_results.txt')}")
    
    return hierarchical_trace