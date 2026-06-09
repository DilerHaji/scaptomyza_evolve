import polars as pl
import argparse
from pathlib import Path
import re

def parse_sample_map(sample_map_path):
    sample_map = {}
    with open(sample_map_path, 'r') as f:
        for line in f:
            number_id, sample_name = line.strip().split(',')
            number = number_id.split('.')[-1]
            sample_map[number] = sample_name
    return sample_map

def process_genomics_data(input_csv, sample_map_path, statistic, output_path):
    sample_map = parse_sample_map(sample_map_path)
    base_cols = ['chrom', 'start', 'end']
    select_exprs = [
        pl.col(col) for col in base_cols
    ]
    df = pl.scan_csv(input_csv)
    stat_pattern = re.compile(f"(\\d+)\\.{statistic}$")
    stat_cols = [col for col in df.collect_schema().names() if stat_pattern.match(col)]
    value_vars = []
    for col in stat_cols:
        sample_num = stat_pattern.match(col).group(1)
        if sample_num in sample_map:
            value_vars.append(col)
    df = df.select(base_cols + value_vars)
    df = df.collect().melt(
        id_vars=base_cols,
        value_vars=value_vars,
        variable_name="variable",
        value_name=statistic
    )
    df = df.with_columns([
        pl.col('variable')
        .map_elements(
            lambda x: sample_map[stat_pattern.match(x).group(1)],
            return_dtype=pl.Utf8
        )
        .alias('sample')
    ]).drop('variable')
    summary_df = df.group_by('sample').agg([
        pl.col(statistic).mean().alias('mean'),
        pl.col(statistic).std().alias('std_dev'),
        pl.col(statistic).count().alias('n_windows')
    ])
    df.select(base_cols + ['sample', statistic]).write_csv(output_path)
    summary_output = str(Path(output_path).with_suffix('')) + "_means.csv"
    summary_df.write_csv(summary_output)

def main():
    parser = argparse.ArgumentParser(description='Process genomics data file')
    parser.add_argument('input_csv', help='Input CSV file path')
    parser.add_argument('sample_map', help='Sample mapping file path')
    parser.add_argument('statistic', help='Statistic to extract (e.g., theta_watterson_abs)')
    parser.add_argument('output_csv', help='Output CSV file path')
    
    args = parser.parse_args()
    
    process_genomics_data(
        args.input_csv,
        args.sample_map,
        args.statistic,
        args.output_csv
    )

if __name__ == "__main__":
    main()