import polars as pl
import argparse
import sys

def create_genomic_matrix(input_file, output_file, value_column):
    # Read the CSV file
    df = pl.read_csv(input_file)

    # Create a unique identifier for each genomic site
    df = df.with_columns(pl.concat_str(["CHROM", "POS"], separator="_").alias("CHROM_POS"))

    # Pivot the dataframe
    pivot_df = df.pivot(
        values=value_column,
        index="CHROM_POS",
        on="SourceFile",
        aggregate_function="first"
    )

    # Write the result to a CSV file
    pivot_df.write_csv(output_file)

def main():
    parser = argparse.ArgumentParser(description="Create a genomic matrix from CSV data.")
    parser.add_argument("input_file", help="Path to the input CSV file")
    parser.add_argument("output_file", help="Path to the output CSV file")
    parser.add_argument("value_column", help="Column to use for matrix values")
    args = parser.parse_args()

    try:
        create_genomic_matrix(args.input_file, args.output_file, args.value_column)
    except Exception as e:
        sys.exit(1)

if __name__ == "__main__":
    main()