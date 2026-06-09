import polars as pl
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Subset a large CSV file based on a list of chromosomes.")
    parser.add_argument("input_file", help="Path to the input CSV file")
    parser.add_argument("output_file", help="Path to the output CSV file")
    # Changed argument name from 'region' to 'contigs' to reflect new usage
    parser.add_argument("contigs", help="Comma-separated list of chromosomes (e.g. chr1,chr2)")
    parser.add_argument("correction", help="Correction column to include in the output")
    args = parser.parse_args()

    target_chroms = args.contigs.split(',')

    try:
        df = pl.scan_csv(args.input_file).select(
            ["CHROM", "POS", "1.FREQ", args.correction, "gen", "trt", "pop", "SourceFile"]
        )
    except Exception as e:
        sys.exit(1)

    filtered_df = df.filter(
        pl.col("CHROM").is_in(target_chroms)
    )

    filtered_df.collect().write_csv(args.output_file)

if __name__ == "__main__":
    main()
    
    
    