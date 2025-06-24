#!/usr/bin/env python3
"""
Generate synthetic test data for FracFixR package with known ground truth.

This script creates RNA-seq count data that simulates a fractionation experiment
with known fraction proportions, allowing validation of FracFixR's estimates.

Usage:
    python generate_data.py --total_reads 1000000 --n_transcripts 10000 
                          --fraction_weights "0.3,0.5" --lost_fraction 0.2
                          --output_dir test_data/

Author: Alice Cleynen, Agin Ravindran, Nikolay Shirokikh
License: CC BY-NC-ND 4.0
"""

import numpy as np
import pandas as pd
import argparse
import os
import json
from typing import List, Tuple, Dict


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate synthetic RNA-seq fractionation data with known ground truth'
    )
    
    parser.add_argument('--total_reads', type=int, default=1000000,
                        help='Total number of reads in the dataset')
    parser.add_argument('--n_transcripts', type=int, default=1000,
                        help='Number of transcripts to simulate')
    parser.add_argument('--fraction_weights', type=str, default="0.3,0.5",
                        help='Comma-separated fraction proportions (must sum to < 1)')
    parser.add_argument('--lost_fraction', type=float, default=0.2,
                        help='Proportion of lost/unrecoverable material')
    parser.add_argument('--n_conditions', type=int, default=2,
                        help='Number of experimental conditions')
    parser.add_argument('--n_replicates', type=int, default=2,
                        help='Number of replicates per condition')
    parser.add_argument('--recovery_rates', type=str, default=None,
                        help='Comma-separated recovery rates for each fraction')
    parser.add_argument('--seed', type=int, default=123,
                        help='Random seed for reproducibility')
    parser.add_argument('--output_dir', type=str, default='test_data',
                        help='Output directory for generated files')
    
    return parser.parse_args()


def generate_transcript_expression(n_transcripts: int, seed: int = 123) -> np.ndarray:
    """
    Generate realistic transcript expression levels using log-normal distribution.
    
    Returns array of expression values for each transcript.
    """
    np.random.seed(seed)
    
    # Use log-normal distribution to simulate realistic expression patterns
    # Most transcripts have low expression, few have very high expression
    log_mean = 4.0  # Mean of log(expression)
    log_std = 2.0   # Standard deviation of log(expression)
    
    expression = np.random.lognormal(log_mean, log_std, n_transcripts)
    
    # Add some zeros (unexpressed genes)
    zero_prob = 0.1
    zero_mask = np.random.random(n_transcripts) < zero_prob
    expression[zero_mask] = 0
    
    return expression


def distribute_to_fractions(expression: np.ndarray, 
                          fraction_weights: List[float],
                          lost_fraction: float,
                          seed: int = 123) -> Dict[str, np.ndarray]:
    """
    Distribute transcript expression across fractions.
    
    Each transcript can have different propensities for each fraction.
    """
    np.random.seed(seed + 1)
    n_transcripts = len(expression)
    n_fractions = len(fraction_weights)
    
    # Ensure weights sum to <= 1
    total_weight = sum(fraction_weights) + lost_fraction
    if total_weight > 1.001:  # Allow small numerical error
        raise ValueError(f"Fraction weights + lost fraction must sum to <= 1, got {total_weight}")
    
    # Generate transcript-specific fraction preferences using Dirichlet distribution
    # This creates realistic variation in how transcripts distribute across fractions
    alpha = np.ones(n_fractions + 1) * 2  # +1 for lost fraction
    transcript_preferences = np.random.dirichlet(alpha, n_transcripts)
    
    # Scale preferences by global fraction weights
    global_weights = fraction_weights + [lost_fraction]
    for i, weight in enumerate(global_weights):
        transcript_preferences[:, i] *= weight
    
    # Normalize so each transcript's preferences sum to 1
    transcript_preferences /= transcript_preferences.sum(axis=1, keepdims=True)
    
    # Distribute expression
    fractions = {}
    for i, fraction_name in enumerate([f"Fraction{i+1}" for i in range(n_fractions)] + ["Lost"]):
        fractions[fraction_name] = expression * transcript_preferences[:, i]
    
    return fractions


def apply_library_effects(counts: np.ndarray, 
                        recovery_rate: float,
                        library_size: int,
                        seed: int = 123) -> np.ndarray:
    """
    Apply library preparation effects: recovery rate and sequencing depth.
    """
    np.random.seed(seed)
    
    # Apply recovery rate (some material is lost during library prep)
    recovered_counts = counts * recovery_rate
    
    # Apply sequencing depth effect (sampling from the library)
    if recovered_counts.sum() > 0:
        # Calculate sampling probabilities
        probs = recovered_counts / recovered_counts.sum()
        
        # Sample reads according to multinomial distribution
        sampled_counts = np.random.multinomial(library_size, probs)
    else:
        sampled_counts = np.zeros_like(counts)
    
    return sampled_counts


def generate_dataset(args) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Generate complete dataset with counts, annotation, and ground truth.
    """
    # Parse fraction weights
    fraction_weights = [float(x) for x in args.fraction_weights.split(',')]
    n_fractions = len(fraction_weights)
    
    # Set default recovery rates if not provided
    if args.recovery_rates is None:
        recovery_rates = [0.8, 0.9] + [0.7] * n_fractions  # Total has higher recovery
    else:
        recovery_rates = [float(x) for x in args.recovery_rates.split(',')]
    
    # Generate base expression levels
    expression = generate_transcript_expression(args.n_transcripts, args.seed)
    
    # Initialize storage
    all_counts = {}
    all_ground_truth = {}
    sample_names = []
    
    # Generate data for each condition and replicate
    sample_id = 0
    for cond_idx in range(args.n_conditions):
        condition = f"Condition{cond_idx + 1}"
        
        # Add condition-specific effects for some transcripts
        cond_expression = expression.copy()
        if cond_idx > 0:
            # Make 5% of transcripts differentially distributed
            n_diff = int(0.05 * args.n_transcripts)
            diff_indices = np.random.choice(args.n_transcripts, n_diff, replace=False)
            
            # Modify expression or fraction preference
            cond_expression[diff_indices] *= np.random.uniform(0.5, 2.0, n_diff)
        
        for rep_idx in range(args.n_replicates):
            replicate = f"Rep{rep_idx + 1}"
            
            # Distribute to fractions with some replicate-specific variation
            rep_seed = args.seed + cond_idx * 100 + rep_idx
            fractions = distribute_to_fractions(
                cond_expression, fraction_weights, args.lost_fraction, rep_seed
            )
            
            # Generate Total sample (sum of all fractions including lost)
            total_counts = sum(fractions.values())
            
            # Apply library effects to Total
            sample_name = f"Sample{sample_id + 1}"
            library_size = int(args.total_reads / (args.n_conditions * args.n_replicates * (n_fractions + 1)))
            all_counts[sample_name] = apply_library_effects(
                total_counts, recovery_rates[0], library_size, rep_seed
            )
            sample_names.append(sample_name)
            sample_id += 1
            
            # Store ground truth for Total
            all_ground_truth[f"{condition}_{replicate}_Total"] = {
                'true_counts': total_counts,
                'fraction_weights': {f"Fraction{i+1}": w for i, w in enumerate(fraction_weights)},
                'lost_fraction': args.lost_fraction
            }
            
            # Apply library effects to each fraction
            for frac_idx, (frac_name, frac_counts) in enumerate(
                [(f"Fraction{i+1}", fractions[f"Fraction{i+1}"]) for i in range(n_fractions)]
            ):
                sample_name = f"Sample{sample_id + 1}"
                all_counts[sample_name] = apply_library_effects(
                    frac_counts, recovery_rates[frac_idx + 1], library_size, rep_seed + frac_idx + 1
                )
                sample_names.append(sample_name)
                sample_id += 1
                
                # Store ground truth for fraction
                all_ground_truth[f"{condition}_{replicate}_{frac_name}"] = {
                    'true_counts': frac_counts,
                    'true_proportion': frac_counts / total_counts
                }
    
    # Create count matrix
    count_matrix = pd.DataFrame(all_counts)
    count_matrix.index = [f"Gene{i+1}" for i in range(args.n_transcripts)]
    
    # Create annotation
    annotation_data = []
    sample_id = 1
    for cond_idx in range(args.n_conditions):
        for rep_idx in range(args.n_replicates):
            # Total sample
            annotation_data.append({
                'Sample': f"Sample{sample_id}",
                'Condition': f"Condition{cond_idx + 1}",
                'Type': 'Total',
                'Replicate': f"Rep{rep_idx + 1}"
            })
            sample_id += 1
            
            # Fraction samples
            for frac_idx in range(n_fractions):
                annotation_data.append({
                    'Sample': f"Sample{sample_id}",
                    'Condition': f"Condition{cond_idx + 1}",
                    'Type': f"Fraction{frac_idx + 1}",
                    'Replicate': f"Rep{rep_idx + 1}"
                })
                sample_id += 1
    
    annotation = pd.DataFrame(annotation_data)
    
    # Create ground truth summary
    ground_truth_summary = pd.DataFrame([
        {
            'parameter': 'fraction_weights',
            'value': ','.join(map(str, fraction_weights))
        },
        {
            'parameter': 'lost_fraction',
            'value': str(args.lost_fraction)
        },
        {
            'parameter': 'recovery_rates',
            'value': ','.join(map(str, recovery_rates))
        },
        {
            'parameter': 'total_reads',
            'value': str(args.total_reads)
        },
        {
            'parameter': 'n_transcripts',
            'value': str(args.n_transcripts)
        }
    ])
    
    return count_matrix, annotation, ground_truth_summary, all_ground_truth


def main():
    """Main function to generate and save synthetic data."""
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"Generating synthetic data with {args.n_transcripts} transcripts...")
    print(f"Fraction weights: {args.fraction_weights}")
    print(f"Lost fraction: {args.lost_fraction}")
    
    # Generate dataset
    counts, annotation, ground_truth_summary, ground_truth_detailed = generate_dataset(args)
    
    # Save files
    counts_file = os.path.join(args.output_dir, 'counts_matrix.csv')
    annotation_file = os.path.join(args.output_dir, 'annotation.csv')
    ground_truth_file = os.path.join(args.output_dir, 'ground_truth.csv')
    ground_truth_json = os.path.join(args.output_dir, 'ground_truth_detailed.json')
    
    counts.to_csv(counts_file)
    annotation.to_csv(annotation_file, index=False)
    ground_truth_summary.to_csv(ground_truth_file, index=False)
    
    # Save detailed ground truth as JSON
    with open(ground_truth_json, 'w') as f:
        json.dump(ground_truth_detailed, f, indent=2, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else x)
    
    print(f"\nGenerated files:")
    print(f"  - {counts_file} ({counts.shape[0]} genes x {counts.shape[1]} samples)")
    print(f"  - {annotation_file}")
    print(f"  - {ground_truth_file}")
    print(f"  - {ground_truth_json}")
    
    # Print summary statistics
    print(f"\nSummary statistics:")
    print(f"  Total counts: {counts.sum().sum():,}")
    print(f"  Mean counts per sample: {counts.sum(axis=0).mean():.0f}")
    print(f"  Genes with zero counts: {(counts.sum(axis=1) == 0).sum()}")
    
    # Validation check
    total_weight = sum([float(x) for x in args.fraction_weights.split(',')]) + args.lost_fraction
    if abs(total_weight - 1.0) > 0.001:
        print(f"\nWarning: Fraction weights + lost fraction = {total_weight:.3f} (should equal 1.0)")


if __name__ == "__main__":
    main()