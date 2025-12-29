#!/usr/bin/env python3
"""
Test script to verify that the complementarity.py duplicate bug is fixed
"""

import sys
import os
import pandas as pd
import pickle

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from analysis.complementarity import analyze_contact_complementarity

# Load test data
with_lipid_dir = 'lipac_results/with_target_lipid'
without_lipid_dir = 'lipac_results/without_target_lipid'

# Choose one to test
test_dir = with_lipid_dir if os.path.exists(with_lipid_dir) else without_lipid_dir

if not os.path.exists(test_dir):
    print("Error: Cannot find test data directory")
    sys.exit(1)

print(f"Testing with data from: {test_dir}")

# Load protein contact results
protein_contact_file = os.path.join(test_dir, 'protein_contacts.pkl')
lipid_contact_file = os.path.join(test_dir, 'lipid_contacts.pkl')

if not os.path.exists(protein_contact_file) or not os.path.exists(lipid_contact_file):
    print("Error: Cannot find contact pickle files")
    sys.exit(1)

with open(protein_contact_file, 'rb') as f:
    protein_contacts = pickle.load(f)

with open(lipid_contact_file, 'rb') as f:
    lipid_contacts = pickle.load(f)

print("\n" + "="*60)
print("Running complementarity analysis with fixed code...")
print("="*60)

# Run the analysis
result_df = analyze_contact_complementarity(protein_contacts, lipid_contacts)

if result_df is not None:
    print("\n" + "="*60)
    print("Analysis Results:")
    print("="*60)
    print(f"Total entries: {len(result_df)}")
    print(f"Unique proteins: {result_df['protein'].nunique()}")

    # Check for duplicates
    duplicates = result_df.duplicated(subset=['protein', 'residue'], keep=False)
    n_duplicates = duplicates.sum()

    print(f"\n🔍 Duplicate Check:")
    if n_duplicates > 0:
        print(f"❌ FAILURE: Found {n_duplicates} duplicate entries")
        duplicate_df = result_df[duplicates].sort_values(['protein', 'residue'])
        print("\nDuplicate examples:")
        print(duplicate_df[['protein', 'residue', 'protein_pair']].head(10))
    else:
        print(f"✅ SUCCESS: No duplicates found!")

    # Count entries per protein
    print(f"\n📊 Entries per protein:")
    protein_counts = result_df['protein'].value_counts()
    for protein, count in protein_counts.items():
        n_unique = result_df[result_df['protein'] == protein]['residue'].nunique()
        print(f"  {protein}: {count} entries, {n_unique} unique residues")
        if count != n_unique:
            print(f"    ⚠️ WARNING: Entry count doesn't match unique residue count!")

    # Save cleaned result
    output_file = os.path.join(test_dir, 'complementarity_fixed.csv')
    result_df.to_csv(output_file, index=False)
    print(f"\n💾 Saved fixed results to: {output_file}")

    # Compare with original if it exists
    original_file = os.path.join(test_dir, 'complementarity.csv')
    if os.path.exists(original_file):
        original_df = pd.read_csv(original_file)
        print(f"\n📈 Comparison with original:")
        print(f"  Original: {len(original_df)} entries")
        print(f"  Fixed: {len(result_df)} entries")
        print(f"  Difference: {len(original_df) - len(result_df)} entries removed")
else:
    print("❌ Analysis failed - no results returned")

print("\n" + "="*60)
print("Test complete!")
print("="*60)