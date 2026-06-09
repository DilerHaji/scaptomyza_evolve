import polars as pl

chunk_files = snakemake.input.chunks
candidate_file = snakemake.input.candidates
output_file = snakemake.output.final_csv

if not chunk_files:
    with open(output_file, 'w') as f:
        f.write("chrom,pos,LRT_chisq,PB_p_val,error,Region_ID\n")
else:
    q = [pl.scan_csv(f) for f in chunk_files]
    results_df = pl.concat(q).collect()
    meta_df = pl.read_csv(candidate_file, separator="\t")
    rename_map = {}
    for col in meta_df.columns:
        if col.lower() == "chrom":
            rename_map[col] = "chrom"
        elif col.lower() == "pos":
            rename_map[col] = "pos"
    
    if rename_map:
        meta_df = meta_df.rename(rename_map)
    
    results_df = results_df.with_columns([
        pl.col("chrom").cast(pl.Utf8).str.strip_chars(),
        pl.col("pos").cast(pl.Int64)
    ])
    
    meta_df = meta_df.with_columns([
        pl.col("chrom").cast(pl.Utf8).str.strip_chars(),
        pl.col("pos").cast(pl.Int64)
    ])

    final_df = results_df.join(meta_df, on=["chrom", "pos"], how="left")
    
    final_df.write_csv(output_file)
