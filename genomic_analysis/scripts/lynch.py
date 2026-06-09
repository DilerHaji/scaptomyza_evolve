import polars as pl
import argparse
from typing import List, Tuple
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process allele frequency and effective sample size data')
    parser.add_argument('--af_file', required=True, help='Path to allele frequency CSV file')
    parser.add_argument('--neff_file', required=True, help='Path to effective sample size CSV file')
    parser.add_argument('--samples', required=True, help='Ordered sample names separated by |')
    parser.add_argument('--output', required=True, help='Output CSV file path (trajectory)')
    parser.add_argument('--mean_output', help='Output CSV file path for mean_s (optional)')
    parser.add_argument('--s_output', help='Optional CSV path to save per-interval s-values')
    return parser.parse_args()

def load_and_validate_data(af_file: str, neff_file: str, samples: List[str]) -> Tuple[pl.DataFrame, pl.DataFrame]:
    af_df = pl.read_csv(af_file)
    neff_df = pl.read_csv(neff_file)
    
    assert af_df.columns == neff_df.columns, "Column mismatch between AF and Neff files"
    assert (af_df['CHROM'] == neff_df['CHROM']).all() and (af_df['POS'] == neff_df['POS']).all(), \
        "CHROM and POS must match between files"
    
    for sample in samples:
        assert sample in af_df.columns and sample in neff_df.columns, f"Sample {sample} not found in data"
    
    return af_df, neff_df

def calculate_s_and_trajectory(af_df: pl.DataFrame, neff_df: pl.DataFrame, samples: List[str]) -> Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    s_values_df = af_df.select(['CHROM', 'POS'])
    traj_df = af_df.select(['CHROM', 'POS'])
    af0 = af_df[samples[0]].clip(1e-10, 1 - 1e-10).round(2)
    prev_af = af0

    for i in range(len(samples) - 1):
        current_sample = samples[i]
        next_sample = samples[i+1]

        af1 = af_df[current_sample].clip(1e-10, 1 - 1e-10)
        af2 = af_df[next_sample].clip(1e-10, 1 - 1e-10)
        neff1 = (neff_df[current_sample] / 2.0)
        neff2 = (neff_df[next_sample] / 2.0)

        term1 = (af2 - af1) / (af1 * (1 - af2))
        term2 = (af2 * (1 - af1)) / (2 * af1 * (1 - af2)) * (
            1 / (neff1 * af1) + 1 / (neff2 * (1 - af2))
        )
        exp1 = -(10 ** (0.51 - (0.000986 * neff1) + (0.0000184 * (neff1 ** 2))))
        exp2 = (-4.78 + (0.134 * neff1) - (0.000927 * (neff1 ** 2))) * af1.log10()
        exp3 = (-1.633 + (0.084 * neff1) - (0.000594 * (neff1 ** 2))) * (af1.log10() ** 2)
        term3 = (exp1 + exp2 + exp3).pow(10)

        s_value = (term1 - term2 + term3).round(2)
        s_values_df = s_values_df.with_columns(s_value.alias(f"s_{current_sample}_{next_sample}"))

        new_af = ((prev_af * (1 + s_value)) / (1 + s_value * prev_af)).round(2)
        traj_df = traj_df.with_columns(new_af.alias(next_sample))
        prev_af = new_af

    s_cols = [c for c in s_values_df.columns if c.startswith("s_")]
    mean_s_df = s_values_df.with_columns(
        pl.mean_horizontal(s_cols).round(2).alias("mean_s")
    ).select(['CHROM', 'POS', 'mean_s'])

    return traj_df, mean_s_df, s_values_df

def main():
    args = parse_arguments()
    samples = args.samples.split('|')

    af_df, neff_df = load_and_validate_data(args.af_file, args.neff_file, samples)

    traj_df, mean_s_df, s_values_df = calculate_s_and_trajectory(af_df, neff_df, samples)

    traj_df.write_csv(args.output)

    if args.mean_output:
        mean_s_df.write_csv(args.mean_output)
    else:
        base_name = os.path.splitext(args.output)[0]
        mean_output_path = f"{base_name}_mean_s.csv"
        mean_s_df.write_csv(mean_output_path)

    if args.s_output:
        s_values_df.write_csv(args.s_output)

if __name__ == "__main__":
    main()