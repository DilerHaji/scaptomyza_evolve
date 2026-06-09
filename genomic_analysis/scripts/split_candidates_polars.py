import polars as pl
import math
import os

if not os.path.exists(snakemake.output.dir):
    os.makedirs(snakemake.output.dir)

df = pl.read_csv(snakemake.input.candidates, separator="\t")
total_rows = df.height
chunk_size = snakemake.params.chunk_size

num_chunks = math.ceil(total_rows / chunk_size)

for i in range(num_chunks):
    offset = i * chunk_size
    chunk = df.slice(offset, chunk_size)
    fname = os.path.join(snakemake.output.dir, f"chunk_{i:06d}.tsv")
    chunk.write_csv(fname, separator="\t")