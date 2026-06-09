# scripts/melt_permuted.py
import argparse
import sys
import polars as pl
import itertools

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--perm-id", type=int, required=True)
    parser.add_argument("--target-reps", required=True, type=lambda x: x.split(','))
    parser.add_argument("--ref-reps", required=True, type=lambda x: x.split(','))
    parser.add_argument("--generations", required=True, type=lambda x: x.split(','))
    return parser.parse_args()

def get_col_name(rep1, gen1, rep2, gen2, cols):
    p1s = [f"{rep1}G{gen1}", f"{rep1}_{gen1}", f"{rep1}{gen1}"]
    p2s = [f"{rep2}G{gen2}", f"{rep2}_{gen2}", f"{rep2}{gen2}"]
    for p1 in p1s:
        for p2 in p2s:
            cands = [f"{p1}:{p2}.fst", f"{p2}:{p1}.fst", f"{p1}.{p2}.fst", f"{p2}.{p1}.fst", 
                     f"fst_{p1}_{p2}", f"fst_{p2}_{p1}"]
            for c in cands:
                if c in cols: return c
    return None

def main():
    args = parse_args()
    
    if len(args.target_reps) != len(args.ref_reps):
         sys.exit("Error: Target and Ref lists must be equal length.")

    refs_sorted = sorted(args.ref_reps)
    all_perms = list(itertools.permutations(refs_sorted))

    if args.perm_id >= len(all_perms):
        sys.exit(f"Error: Perm ID {args.perm_id} out of range (Max {len(all_perms)-1})")

    current_ref_order = all_perms[args.perm_id]
    pairs = list(zip(args.target_reps, current_ref_order))

    try:
        lf = pl.scan_csv(args.input, null_values=["NA", "nan", "."], infer_schema_length=0)
        available_cols = set(lf.collect_schema().names())
    except Exception as e:
        sys.exit(f"Error reading input schema: {e}")

    coord_cols = [c for c in ["chrom", "start", "end"] if c in available_cols]

    expressions = []
    
    found_any = False

    for gen in args.generations:
        for (t, r) in pairs:
            col = get_col_name(t, gen, r, gen, available_cols)
            if col:
                found_any = True
                chunk_lf = lf.select(
                    [pl.col(c) for c in coord_cols] + 
                    [pl.col(col).cast(pl.Float64, strict=False).alias("fst")]
                ).with_columns([
                    pl.lit(int(gen)).alias("gen"),
                    pl.lit(f"{t}_vs_{r}").alias("rep_pair")
                ])
                
                expressions.append(chunk_lf)

    if not found_any:
        sys.exit("No matching columns found for this permutation.")

    try:
        final_df = pl.concat(expressions).filter(pl.col("fst").is_not_null()).collect(streaming=True)
        final_df.write_csv(args.output)
    except Exception as e:
        sys.exit(f"Error during materialization: {e}")

if __name__ == "__main__":
    main()