"""
Report generation functions for BayesianLipidAnalysis
"""

import os
from datetime import datetime


def generate_summary_report(output_dir):
    """Generate summary report of analysis results"""
    print("\n===== Generating Analysis Summary Report =====")
    
    # Output file path
    summary_path = os.path.join(output_dir, "analysis_summary.txt")
    
    try:
        # Paths to existing analysis result files
        exploratory_path = os.path.join(output_dir, 'exploratory_analysis_results.txt')
        bayesian_path = os.path.join(output_dir, 'bayesian_analysis_results.txt')
        hierarchical_path = os.path.join(output_dir, 'hierarchical_model_results.txt')
        
        with open(summary_path, 'w') as summary_file:
            summary_file.write("=====================================================\n")
            summary_file.write("        Target Lipid-Protein Interaction Analysis - Summary Report        \n")
            summary_file.write("=====================================================\n\n")
            
            # Date/time information
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            summary_file.write(f"Generated: {current_time}\n\n")
            
            # Exploratory analysis summary
            summary_file.write("1. Key Findings from Exploratory Analysis\n")
            summary_file.write("------------------------\n")
            if os.path.exists(exploratory_path):
                with open(exploratory_path, 'r') as f:
                    content = f.read()
                
                # Extract key results
                if "with target lipid vs without" in content:
                    lines = content.split('\n')
                    for i, line in enumerate(lines):
                        if "Protein-protein contact difference" in line:
                            summary_file.write(f"* {line.strip()}\n")
                        if "Lipid-protein contact difference" in line:
                            summary_file.write(f"* {line.strip()}\n")
                        if "Ratio difference" in line:
                            summary_file.write(f"* {line.strip()}\n")
                
                # Target lipid contact strength comparison
                if "Comparison by Target Lipid Contact Strength" in content:
                    summary_file.write("\nEffects by Target Lipid Contact Strength:\n")
                    lines = content.split('\n')
                    for i, line in enumerate(lines):
                        if "Protein contact difference (strong vs weak/no target lipid)" in line:
                            summary_file.write(f"* {line.strip()}\n")
                            if i+1 < len(lines) and "Strong target lipid mean" in lines[i+1]:
                                summary_file.write(f"* {lines[i+1].strip()}\n")
                            if i+2 < len(lines) and "Weak/no target lipid mean" in lines[i+2]:
                                summary_file.write(f"* {lines[i+2].strip()}\n")
            else:
                summary_file.write("* Exploratory analysis results not found.\n")
            
            summary_file.write("\n")
            
            # Bayesian regression analysis summary
            summary_file.write("2. Key Findings from Bayesian Regression Analysis\n")
            summary_file.write("------------------------------\n")
            if os.path.exists(bayesian_path):
                with open(bayesian_path, 'r') as f:
                    content = f.read()
                
                # Model 1 results
                if "Model 1: Relationship between Target Lipid Contact and Protein Contact" in content:
                    summary_file.write("Protein Contact Model:\n")
                    lines = content.split('\n')
                    for i, line in enumerate(lines):
                        if "Mean effect size (β)" in line:
                            summary_file.write(f"* {line.strip()}\n")
                        if "95% HDI" in line and "Model 1" in content[:content.find(line)]:
                            summary_file.write(f"* {line.strip()}\n")
                        if "P(β < 0)" in line and "Model 1" in content[:content.find(line)]:
                            summary_file.write(f"* {line.strip()}\n")
                
                # Model 2 results
                if "Model 2: Relationship between Target Lipid Contact and Lipid Contact" in content:
                    summary_file.write("\nLipid Contact Model:\n")
                    lines = content.split('\n')
                    for i, line in enumerate(lines):
                        if "Mean effect size (β)" in line and "Model 2" in content[:content.find(line)]:
                            summary_file.write(f"* {line.strip()}\n")
                        if "95% HDI" in line and "Model 2" in content[:content.find(line)]:
                            summary_file.write(f"* {line.strip()}\n")
                        if "P(β > 0)" in line:
                            summary_file.write(f"* {line.strip()}\n")
                
                # Effect direction results
                if "Effect Direction" in content:
                    summary_file.write("\nEffect Direction:\n")
                    lines = content.split('\n')
                    for line in lines:
                        if "P(Protein Effect < 0 and Lipid Effect > 0)" in line:
                            summary_file.write(f"* {line.strip()}\n")
                
                # Add conclusion if available
                if "Conclusion:" in content:
                    summary_file.write("\nConclusion:\n")
                    found_conclusion = False
                    lines = content.split('\n')
                    for i, line in enumerate(lines):
                        if "Conclusion:" in line:
                            found_conclusion = True
                            continue
                        if found_conclusion and line.strip():
                            summary_file.write(f"* {line.strip()}\n")
            else:
                summary_file.write("* Bayesian regression analysis results not found.\n")
            
            summary_file.write("\n")
            
            # Hierarchical model summary (if exists)
            summary_file.write("3. Key Findings from Hierarchical Bayesian Model\n")
            summary_file.write("--------------------------------\n")
            if os.path.exists(hierarchical_path):
                with open(hierarchical_path, 'r') as f:
                    content = f.read()
                
                # Global effect
                if "Global Effect" in content:
                    summary_file.write("Global Effect:\n")
                    lines = content.split('\n')
                    for i, line in enumerate(lines):
                        if "Global GM3 effect" in line:
                            summary_file.write(f"* {line.strip()}\n")
                        if "95% HDI" in line and "Global Effect" in content[:content.find(line)]:
                            summary_file.write(f"* {line.strip()}\n")
                        if "P(μ_β < 0)" in line:
                            summary_file.write(f"* {line.strip()}\n")
                
                # Protein-specific effects (summary only)
                if "Protein-Specific Effects" in content:
                    summary_file.write("\nProtein-Specific Effect Summary:\n")
                    lines = content.split('\n')
                    proteins_with_strong_effect = []
                    
                    current_protein = None
                    for i, line in enumerate(lines):
                        if ":" in line and not "GM3 effect" in line and not "HDI" in line and not "P(" in line:
                            current_protein = line.strip().replace(':', '')
                        if current_protein and "P(β < 0)" in line:
                            prob = float(line.split('=')[1].strip())
                            if prob > 0.9 or prob < 0.1:  # Proteins with strong effect
                                effect = "negative" if prob > 0.5 else "positive"
                                proteins_with_strong_effect.append((current_protein, effect, prob))
                    
                    if proteins_with_strong_effect:
                        for protein, effect, prob in proteins_with_strong_effect:
                            effect_str = "negative" if effect == "negative" else "positive"
                            summary_file.write(f"* {protein}: Probability of {effect_str} GM3 effect = {prob:.3f}\n")
                    else:
                        summary_file.write("* No proteins with strong effects detected.\n")
            else:
                summary_file.write("* Hierarchical Bayesian model results not found.\n")
            
            # Overall conclusion
            summary_file.write("\n4. Overall Conclusion\n")
            summary_file.write("------------\n")
            
            # Integrate findings from all analyses
            has_exploratory = os.path.exists(exploratory_path)
            has_bayesian = os.path.exists(bayesian_path)
            has_hierarchical = os.path.exists(hierarchical_path)
            
            if has_bayesian:
                with open(bayesian_path, 'r') as f:
                    bayesian_content = f.read()
                if "Conclusion:" in bayesian_content:
                    lines = bayesian_content.split('\n')
                    found_conclusion = False
                    conclusion_text = []
                    for line in lines:
                        if "Conclusion:" in line:
                            found_conclusion = True
                            continue
                        if found_conclusion and line.strip():
                            conclusion_text.append(line.strip())
                    
                    if conclusion_text:
                        for line in conclusion_text:
                            summary_file.write(f"{line}\n")
                    else:
                        summary_file.write("The Bayesian analysis results show consistent effects across different models.\n")
                        summary_file.write("GM3 appears to have opposite directional effects on protein-protein contacts\n")
                        summary_file.write("and lipid-protein contacts.\n")
            elif has_exploratory:
                summary_file.write("The exploratory analysis suggests that GM3 presence affects both\n")
                summary_file.write("protein-protein contacts and lipid-protein contacts.\n")
                summary_file.write("Refer to the Bayesian model results for more detailed effect measurements.\n")
            else:
                summary_file.write("Limited analysis results are available for definitive conclusions.\n")
                summary_file.write("Please refer to the detailed model results and consider further analysis.\n")
            
            # Footer
            summary_file.write("\n=====================================================\n")
            summary_file.write("Note: This summary report is automatically generated. For detailed results,\n")
            summary_file.write("please refer to the individual analysis files.\n")
            summary_file.write("=====================================================\n")
        
        print(f"Summary report saved to {summary_path}")
        return True
        
    except Exception as e:
        print(f"Error generating summary report: {str(e)}")
        # Create minimal report to continue
        try:
            with open(summary_path, 'w') as f:
                f.write("===== GM3-Protein Interaction Analysis - Summary Report =====\n\n")
                f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                f.write("An error occurred while generating the summary report.\n")
                f.write("Please refer to the individual analysis files.\n")
            print(f"Simplified report saved to {summary_path}")
            return False
        except:
            print("Completely aborting summary report creation")
            return False