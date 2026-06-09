import polars as pl
import sys

def create_afmat_samps_sites(df, cols_for_samps, freq_col, neff_column, output_prefix):
    sites = df.select(['CHROM', 'POS']).unique()
    sites = sites.sort(['CHROM', 'POS'])
    
    sites.write_csv(output_prefix + '_sites.csv', include_header=True)
    afmat = df.pivot(
        index=['CHROM', 'POS'], 
        on='SourceFile', 
        values=freq_col,
        aggregate_function=None # Fail if duplicates exist
    ).sort(['CHROM', 'POS']).fill_null(0)

    if not afmat.select(['CHROM', 'POS']).equals(sites):
        raise ValueError("Critical Error: AF Matrix row order does not match Sites order.")

    sample_order = afmat.columns[2:]
    
    afmat.select(sample_order).write_csv(output_prefix + '_afmat.csv', include_header=True)

    neff_mat = df.pivot(
        index=['CHROM', 'POS'], 
        on='SourceFile', 
        values=neff_column,
        aggregate_function=None
    ).sort(['CHROM', 'POS']).fill_null(0)

    if not neff_mat.select(['CHROM', 'POS']).equals(sites):
        raise ValueError("Critical Error: Neff Matrix row order does not match Sites order.")

    if neff_mat.columns[2:] != sample_order:
        raise ValueError("Critical Error: Sample order differs between Neff mat and AF mat")

    neff_mat.select(sample_order).write_csv(output_prefix + '_neff_mat.csv', include_header=True)

    unique_meta = df.select(cols_for_samps).unique()
    order_df = pl.DataFrame({"SourceFile": sample_order})
    
    samps = order_df.join(unique_meta, on="SourceFile", how="left")

    if samps.height != len(sample_order):
        raise ValueError(f"Metadata rows ({samps.height}) do not match matrix columns ({len(sample_order)})")

    samps.write_csv(output_prefix + '_samps.csv', include_header=True)
    
    with open(output_prefix + '_sample_order.txt', 'w') as f:
        f.write('\n'.join(sample_order))

if __name__ == "__main__":
    try:
        file_path = sys.argv[1]
        freq_col = sys.argv[2]
        neff_column = sys.argv[3]
        output_prefix = sys.argv[4]
        cols_for_samps = [x.strip() for x in sys.argv[5].split(',')]

        df = pl.read_csv(file_path)
        
        create_afmat_samps_sites(df, cols_for_samps, freq_col, neff_column, output_prefix)
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        sys.exit(1)
