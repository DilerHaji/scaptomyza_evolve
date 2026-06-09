import polars as pl
import sys

inputs = snakemake.input
output = snakemake.output[0]

dfs = []
for f in inputs:
    try:
        d = pl.read_csv(f)
        if d.height > 0:
            dfs.append(d)
    except:
        pass

if dfs:
    pl.concat(dfs).write_csv(output)
else:
    with open(output, 'w') as f:
        f.write("chrom,start,end,slope_lmm,se_lmm,z_score,p_value\n")