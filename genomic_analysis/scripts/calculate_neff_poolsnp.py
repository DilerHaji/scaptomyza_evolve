import pandas as pd
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True)
parser.add_argument('--header', type=str, required=True)
parser.add_argument('--poolsizes', type=str, required=True)
parser.add_argument('--pool', type=str, required=True)
parser.add_argument('--output', type=str, required=True)

args = parser.parse_args()

def load_poolsizes(poolsizes):
    pool_dict = {}
    try:
        with open(poolsizes, 'r') as file:
            for line in file:
                parts = line.strip().split(',')
                if len(parts) == 2:
                    pool_dict[parts[0]] = int(parts[1])
    except Exception as e:
        pass
    return pool_dict

def vectorized_calculate_values(df, n, col_name):
    if col_name not in df.columns:
        sys.exit(1)

    df_copy = df.copy()
    df_copy['x_values'] = pd.to_numeric(df_copy[col_name].astype(str).str.split(':', expand=True)[3], errors='coerce')
    valid_x = df_copy['x_values'] > 0
    df_copy.loc[valid_x, 'Czech2023'] = np.round((n * (1 - (((n - 1) / n)**df_copy.loc[valid_x, 'x_values']))), 0)
    df_copy.loc[valid_x, 'Bergland2014'] = np.round((((n*df_copy.loc[valid_x, 'x_values']) -1) / (n + df_copy.loc[valid_x, 'x_values'])), 0)
    df_copy.drop(columns=['x_values'], inplace=True)
    return df_copy

def get_header_from_file(header_file):
    with open(header_file, 'r') as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.strip().split('\t')
    sys.exit(1)

def process_file(filename, header_file, pool, output_filename, poolsizes, chunksize=100000):
    pool_dict = load_poolsizes(poolsizes)
    n = pool_dict.get(pool)

    if n is None:
        sys.exit(1)

    header = get_header_from_file(header_file)

    if pool not in header:
        sys.exit(1)

    first_chunk = True
    
    for chunk in pd.read_csv(filename, chunksize=chunksize, delimiter='\t', header=None, names=header):
        calculated_chunk = vectorized_calculate_values(chunk, n, pool)
        output_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', pool, 'Czech2023', 'Bergland2014']
        final_chunk = calculated_chunk[output_columns].copy()
        final_chunk_cleaned = final_chunk.dropna(subset=['Czech2023', 'Bergland2014'])
        mode = 'w' if first_chunk else 'a'
        header = first_chunk
        final_chunk_cleaned.to_csv(output_filename, mode=mode, index=False, header=header, sep='\t')
        first_chunk = False

if __name__ == "__main__":
    process_file(args.input, args.header, args.pool, args.output, args.poolsizes)
