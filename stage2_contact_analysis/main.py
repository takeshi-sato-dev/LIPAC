#!/usr/bin/env python
"""
Bayesian analysis of lipid-protein interactions
Main execution script

Author: Takeshi Sato, PhD
Kyoto Pharmaceutical University
2024
"""

import os
import sys
import argparse
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import from local package structure
import config  # config.pyはメインディレクトリにある
from core.data_loader import load_data, prepare_residue_analysis_data
from analysis.exploratory import exploratory_analysis
from analysis.bayesian_models import bayesian_regression_model, hierarchical_model
from analysis.residue_analysis import analyze_residue_contacts
from analysis.lipid_analysis import analyze_lipid_specific_effects, analyze_per_peptide_effects
from analysis.causal_analysis import load_causal_data, perform_causal_bayesian_analysis
from utils.reports import generate_summary_report


def main():
    """Main execution function - Bayesian statistical analysis of lipid-protein interactions"""
    
    parser = argparse.ArgumentParser(description='Bayesian statistical analysis of lipid-protein interactions')
    parser.add_argument('--with-lipid-data', default=config.DEFAULT_WITH_LIPID_DATA, help='CSV file for system with target lipid')
    parser.add_argument('--without-lipid-data', default=config.DEFAULT_WITHOUT_LIPID_DATA, help='CSV file for system without target lipid')
    parser.add_argument('--output-dir', default=config.DEFAULT_OUTPUT_DIR, help='Output directory for analysis results')
    parser.add_argument('--skip-exploratory', action='store_true', help='Skip exploratory analysis')
    parser.add_argument('--skip-bayesian', action='store_true', help='Skip Bayesian regression analysis')
    parser.add_argument('--skip-hierarchical', action='store_true', help='Skip hierarchical Bayesian model')
    parser.add_argument('--skip-residue-analysis', action='store_true', help='Skip residue-level analysis')
    parser.add_argument('--skip-lipid-specific', action='store_true', help='Skip lipid-specific analysis')
    parser.add_argument('--skip-per-peptide', action='store_true', help='Skip per-peptide analysis')
    parser.add_argument('--skip-causal', action='store_true', help='Skip causal Bayesian analysis')
    parser.add_argument('--causal-data-dir', default=None, help='Directory containing GM3 causal data from stage1')
    parser.add_argument('--mcmc-samples', type=int, default=config.MCMC_SAMPLES, help='Number of MCMC samples')
    parser.add_argument('--tune-samples', type=int, default=config.TUNE_SAMPLES, help='Number of tuning samples')
    parser.add_argument('--chains', type=int, default=config.CHAINS, help='Number of chains')
    parser.add_argument('--seed', type=int, default=config.RANDOM_SEED, help='Random seed')
    parser.add_argument('--verbose', action='store_true', help='Display detailed output')
    
    args = parser.parse_args()
    
    # Update config with command-line parameters
    config.MCMC_SAMPLES = args.mcmc_samples
    config.TUNE_SAMPLES = args.tune_samples
    config.CHAINS = args.chains
    config.RANDOM_SEED = args.seed
    
    # Set random seed (for reproducibility)
    np.random.seed(config.RANDOM_SEED)
    
    print("===== Bayesian Statistical Analysis of Lipid-Protein Interactions =====")
    print(f"MCMC Samples: {config.MCMC_SAMPLES}, Tune: {config.TUNE_SAMPLES}, Chains: {config.CHAINS}")
    print(f"Random Seed: {config.RANDOM_SEED}")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Create subdirectories
    subdirs = ['figures', 'models', 'residue_analysis', 'lipid_specific_analysis', 'causal_analysis']
    for subdir in subdirs:
        os.makedirs(os.path.join(args.output_dir, subdir), exist_ok=True)
    
    # Load data
    print("\nLoading data...")
    combined_df, matched_df = load_data(args.with_lipid_data, args.without_lipid_data)
    
    if combined_df is None or matched_df is None:
        print("Data loading failed. Aborting analysis.")
        return 1
    
    # Save data summaries
    print("\nSaving data summaries...")
    combined_df.to_csv(os.path.join(args.output_dir, 'combined_data.csv'), index=False)
    matched_df.to_csv(os.path.join(args.output_dir, 'matched_data.csv'), index=False)
    
    # Track completed analyses
    completed_analyses = []
    
    # 1. Exploratory data analysis
    if not args.skip_exploratory:
        print("\n" + "="*50)
        print("Starting Exploratory Analysis")
        print("="*50)
        try:
            exploratory_analysis(combined_df, matched_df, args.output_dir)
            completed_analyses.append("exploratory")
        except Exception as e:
            print(f"Error in exploratory analysis: {str(e)}")
            if args.verbose:
                import traceback
                traceback.print_exc()
    else:
        print("\nExploratory analysis was skipped")
    
    # 2. Bayesian regression analysis
    traces = None
    if not args.skip_bayesian:
        print("\n" + "="*50)
        print("Starting Bayesian Regression Analysis")
        print("="*50)
        try:
            traces = bayesian_regression_model(matched_df, args.output_dir)
            if traces:
                completed_analyses.append("bayesian_regression")
        except Exception as e:
            print(f"Error in Bayesian regression analysis: {str(e)}")
            if args.verbose:
                import traceback
                traceback.print_exc()
    else:
        print("\nBayesian regression analysis was skipped")
    
    # 3. Hierarchical Bayesian model
    hierarchical_trace = None
    if not args.skip_hierarchical:
        print("\n" + "="*50)
        print("Starting Hierarchical Bayesian Model")
        print("="*50)
        try:
            hierarchical_trace = hierarchical_model(matched_df, args.output_dir)
            if hierarchical_trace:
                completed_analyses.append("hierarchical")
        except Exception as e:
            print(f"Error in hierarchical Bayesian model: {str(e)}")
            if args.verbose:
                import traceback
                traceback.print_exc()
    else:
        print("\nHierarchical Bayesian model was skipped")
    
    # 4. Residue-level analysis
    if not args.skip_residue_analysis:
        print("\n" + "="*50)
        print("Starting Residue-Level Analysis")
        print("="*50)
        try:
            # Prepare data for residue analysis
            with_lipid_df = combined_df[combined_df['has_target_lipid'] == True].copy()
            without_lipid_df = combined_df[combined_df['has_target_lipid'] == False].copy()
            
            # Check and prepare partner protein information
            if 'partner_protein' not in with_lipid_df.columns:
                print("Preprocessing data for residue analysis...")
                with_lipid_df = prepare_residue_analysis_data(with_lipid_df)
                without_lipid_df = prepare_residue_analysis_data(without_lipid_df)
            
            # Run residue-level analysis
            residue_result = analyze_residue_contacts(with_lipid_df, without_lipid_df, args.output_dir)
            if residue_result:
                completed_analyses.append("residue_analysis")
            
        except Exception as e:
            print(f"Error in residue-level analysis: {str(e)}")
            if args.verbose:
                import traceback
                traceback.print_exc()
    else:
        print("\nResidue-level analysis was skipped")
    
    # 5. Lipid-specific analysis
    if not args.skip_lipid_specific:
        print("\n" + "="*50)
        print("Starting Lipid-Specific Analysis")
        print("="*50)
        try:
            lipid_result = analyze_lipid_specific_effects(combined_df, matched_df, args.output_dir)
            if lipid_result:
                completed_analyses.append("lipid_specific")
        except Exception as e:
            print(f"Error in lipid-specific analysis: {str(e)}")
            if args.verbose:
                import traceback
                traceback.print_exc()
    else:
        print("\nLipid-specific analysis was skipped")
    
    # 6. Per-peptide analysis
    if not args.skip_per_peptide:
        print("\n" + "="*50)
        print("Starting Per-Peptide Analysis")
        print("="*50)
        try:
            per_peptide_result = analyze_per_peptide_effects(combined_df, matched_df, args.output_dir)
            if per_peptide_result:
                completed_analyses.append("per_peptide")
        except Exception as e:
            print(f"Error in per-peptide analysis: {str(e)}")
            if args.verbose:
                import traceback
                traceback.print_exc()
    else:
        print("\nPer-peptide analysis was skipped")
    
    # 7. Causal Bayesian analysis
    if not args.skip_causal:
        print("\n" + "="*50)
        print("Starting Causal Bayesian Analysis")
        print("="*50)
        try:
            # Determine causal data directory
            if args.causal_data_dir:
                causal_data_dir = args.causal_data_dir
            else:
                # Try default locations based on stage1 config
                possible_dirs = [
                    config.DEFAULT_CAUSAL_DATA_DIR,  # Default from config
                    os.path.join(os.path.dirname(args.output_dir), 'stage1_output'),
                ]
                
                # Try to import stage1 config to get the correct path
                try:
                    stage1_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'stage1_contact_analysis')
                    sys.path.insert(0, stage1_path)
                    import config as stage1_config
                    
                    # Add stage1 output directory to possible locations
                    if hasattr(stage1_config, 'DEFAULT_WITH_LIPID_OUTPUT'):
                        possible_dirs.append(stage1_config.DEFAULT_WITH_LIPID_OUTPUT)
                    if hasattr(stage1_config, 'BASE_OUTPUT_DIR'):
                        possible_dirs.append(os.path.join(stage1_config.BASE_OUTPUT_DIR, 'with_target_lipid'))
                        
                    sys.path.remove(stage1_path)
                except Exception as e:
                    print(f"Warning: Could not import stage1 config: {e}")
                
                # Add fallback locations
                possible_dirs.extend([
                    '../stage1_contact_analysis/lipac_results/with_target_lipid',
                    'lipac_results/with_target_lipid'
                ])
                causal_data_dir = None
                for dir_path in possible_dirs:
                    if os.path.exists(dir_path):
                        # Check if target lipid causal data files exist
                        import glob
                        pattern = os.path.join(dir_path, '*_causal_data_*.csv')
                        if glob.glob(pattern):
                            causal_data_dir = dir_path
                            break
                
                if causal_data_dir is None:
                    print("Warning: No target lipid causal data directory found. Please specify --causal-data-dir")
                    print("Skipping causal analysis. This is normal if stage1 hasn't been re-run after the GM3 -> target_lipid modifications.")
                    completed_analyses.append("causal_skipped")
                    causal_data_dir = None
            
            if causal_data_dir:
                print(f"Loading target lipid causal data from: {causal_data_dir}")
                
                # Determine target lipid name from config or files
                target_lipid = getattr(config, 'TARGET_LIPID', 'DPG3')
                
                # Load causal data
                causal_data = load_causal_data(causal_data_dir, target_lipid.lower())
                
                if causal_data:
                    # Perform causal analysis
                    causal_results = perform_causal_bayesian_analysis(causal_data, args.output_dir, target_lipid)
                    
                    if causal_results:
                        completed_analyses.append("causal_analysis")
                        print("Causal Bayesian analysis completed successfully!")
                    else:
                        print("Causal Bayesian analysis failed.")
                else:
                    print("Failed to load causal data. This is expected if stage1 was run before the GM3->target_lipid modifications.")
                    print("To enable causal analysis, please re-run stage1 with the updated code.")
                    completed_analyses.append("causal_data_unavailable")
                    
        except Exception as e:
            print(f"Error in causal Bayesian analysis: {str(e)}")
            if args.verbose:
                import traceback
                traceback.print_exc()
    else:
        print("\nCausal Bayesian analysis was skipped")
    
    # 8. Generate summary report
    print("\n" + "="*50)
    print("Generating Summary Report")
    print("="*50)
    try:
        summary_result = generate_summary_report(args.output_dir)
        if summary_result:
            completed_analyses.append("summary_report")
    except Exception as e:
        print(f"Error generating summary report: {str(e)}")
        if args.verbose:
            import traceback
            traceback.print_exc()
    
    # Final summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nCompleted analyses: {', '.join(completed_analyses)}")
    print(f"\nResults saved in: {os.path.abspath(args.output_dir)}")
    
    # List generated files
    print("\nGenerated files:")
    for root, dirs, files in os.walk(args.output_dir):
        level = root.replace(args.output_dir, '').count(os.sep)
        indent = ' ' * 2 * level
        subdir = os.path.basename(root)
        if subdir:
            print(f"{indent}{subdir}/")
        sub_indent = ' ' * 2 * (level + 1)
        for file in sorted(files):
            if file.endswith(('.png', '.csv', '.txt', '.nc', '.ps', '.svg')):
                print(f"{sub_indent}{file}")
    
    # Return status
    if len(completed_analyses) == 0:
        print("\nWARNING: No analyses were completed successfully!")
        return 1
    elif len(completed_analyses) < 5 and not any(args.__dict__.get(f'skip_{x}', False) for x in ['exploratory', 'bayesian', 'hierarchical', 'residue_analysis', 'lipid_specific']):
        print(f"\nWARNING: Only {len(completed_analyses)} out of 5 analyses completed successfully.")
        return 2
    else:
        print("\nAll requested analyses completed successfully!")
        return 0


if __name__ == "__main__":
    sys.exit(main())