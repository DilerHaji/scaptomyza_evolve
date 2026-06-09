#!/usr/bin/env python3
import argparse
import sys
import polars as pl
import itertools
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Raw FST CSV")
    parser.add_argument("--output", required=True)
    parser.add_argument("--target-reps", required=True)
    parser.add_argument("--ref-reps", required=True)
    parser.add_argument("--outgroup-reps", required=True, help="Replicates for the outgroup (e.g. Mixture)")
    parser.add_argument("--generations", required=True)
    parser.add_argument("--outgroup-gen", default="00", help="Generation of outgroup (use 'dynamic' if it evolves)")
    parser.add_argument("--regions", help="Semicolon-separated regions")
    parser.add_argument("--founder-reps", default="", help="Founder replicates (e.g. F1,F2,F3,F4). If provided and founder-gen appears in generations, use first 3 distinct founder reps as shared PBS anchor at gen 0.")
    parser.add_argument("--founder-gen", default="00", help="Generation code for the founder timepoint (default: 00)")
    parser.add_argument("--mode", default="combinatorial", choices=["combinatorial", "per_rep"],
                        help="combinatorial: all T×B×M trios (legacy). per_rep: matched Ti↔Bi pairs with pooled B/M averages (rigorous).")
    return parser.parse_args()

def safe_log_transform(fst_col):
    return -1 * (1 - fst_col.clip(0.0, 0.999)).log()

def get_fst_col_name(schema_cols, rep1, gen1, rep2, gen2):
    p1 = f"{rep1}G{gen1}"; p2 = f"{rep2}G{gen2}"
    c1 = f"{p1}:{p2}.fst"; c2 = f"{p2}:{p1}.fst"
    if c1 in schema_cols: return c1
    elif c2 in schema_cols: return c2

    p1_u = f"{rep1}_{gen1}"; p2_u = f"{rep2}_{gen2}"
    c3 = f"{p1_u}:{p2_u}.fst"; c4 = f"{p2_u}:{p1_u}.fst"
    if c3 in schema_cols: return c3
    elif c4 in schema_cols: return c4
    return None

def avg_fst_expr(col_names):
    exprs = [pl.col(c).cast(pl.Float64) for c in col_names]
    return pl.sum_horizontal(exprs) / len(exprs)

def _build_combinatorial_frames(lf, valid_combinations, output_schema_cols):
    frames = []
    for combo in valid_combinations:
        q = lf.select([
            pl.col("chrom"), pl.col("start"), pl.col("end"),
            pl.col(combo["col_TB"]).cast(pl.Float64).alias("fst_tb"),
            pl.col(combo["col_TM"]).cast(pl.Float64).alias("fst_tm"),
            pl.col(combo["col_BM"]).cast(pl.Float64).alias("fst_bm")
        ]).with_columns([
            safe_log_transform(pl.col("fst_tb")).alias("T_tb"),
            safe_log_transform(pl.col("fst_tm")).alias("T_tm"),
            safe_log_transform(pl.col("fst_bm")).alias("T_bm"),
        ]).with_columns([
            ((pl.col("T_tb") + pl.col("T_tm") - pl.col("T_bm")) / 2).alias("pbs_target"),
            ((pl.col("T_tb") + pl.col("T_bm") - pl.col("T_tm")) / 2).alias("pbs_ref"),
            ((pl.col("T_tm") + pl.col("T_bm") - pl.col("T_tb")) / 2).alias("pbs_outgroup")
        ]).with_columns([
            (pl.col("pbs_target") / (1 + pl.col("pbs_target") + pl.col("pbs_ref") + pl.col("pbs_outgroup"))).alias("pbsn1_target"),
            (pl.col("pbs_ref") / (1 + pl.col("pbs_target") + pl.col("pbs_ref") + pl.col("pbs_outgroup"))).alias("pbsn1_ref"),
            pl.lit(combo["gen"]).cast(pl.Int32).alias("gen"),
            pl.lit(combo["lineage"]).alias("lineage")
        ])
        q = q.select([
            "chrom", "start", "end", "gen", "lineage",
            pl.col("T_tb").alias("total_divergence"),
            "T_tm", "T_bm",
            "pbs_target", "pbs_ref",
            "pbsn1_target", "pbsn1_ref"
        ])
        frames.append(q)
    return frames

def _build_per_rep_frames(lf, per_rep_combos):
    frames = []
    for combo in per_rep_combos:
        q = lf.select([
            pl.col("chrom"), pl.col("start"), pl.col("end"),
            avg_fst_expr(combo["ti_b_cols"]).alias("fst_ti_b"),
            avg_fst_expr(combo["ti_m_cols"]).alias("fst_ti_m"),
            avg_fst_expr(combo["bm_cols"]).alias("fst_bm"),
            avg_fst_expr(combo["bi_t_cols"]).alias("fst_bi_t"),
            avg_fst_expr(combo["bi_m_cols"]).alias("fst_bi_m"),
            avg_fst_expr(combo["tm_cols"]).alias("fst_tm_ref"),
        ]).with_columns([
            safe_log_transform(pl.col("fst_ti_b")).alias("T_ti_b"),
            safe_log_transform(pl.col("fst_ti_m")).alias("T_ti_m"),
            safe_log_transform(pl.col("fst_bm")).alias("T_bm"),
            safe_log_transform(pl.col("fst_bi_t")).alias("T_bi_t"),
            safe_log_transform(pl.col("fst_bi_m")).alias("T_bi_m"),
            safe_log_transform(pl.col("fst_tm_ref")).alias("T_tm_ref"),
        ]).with_columns([
            ((pl.col("T_ti_b") + pl.col("T_ti_m") - pl.col("T_bm")) / 2).alias("pbs_target"),
            ((pl.col("T_bi_t") + pl.col("T_bi_m") - pl.col("T_tm_ref")) / 2).alias("pbs_ref"),
            ((pl.col("T_ti_m") + pl.col("T_bm") - pl.col("T_ti_b")) / 2).alias("pbs_outgroup"),
        ]).with_columns([
            (pl.col("pbs_target") / (1 + pl.col("pbs_target") + pl.col("pbs_ref") + pl.col("pbs_outgroup"))).alias("pbsn1_target"),
            (pl.col("pbs_ref") / (1 + pl.col("pbs_target") + pl.col("pbs_ref") + pl.col("pbs_outgroup"))).alias("pbsn1_ref"),
            pl.lit(combo["gen"]).cast(pl.Int32).alias("gen"),
            pl.lit(combo["lineage"]).alias("lineage"),
        ])
        q = q.select([
            "chrom", "start", "end", "gen", "lineage",
            pl.col("T_ti_b").alias("total_divergence"),  # Ti vs pooled B
            pl.col("T_ti_m").alias("T_tm"),              # Ti vs pooled M
            pl.col("T_bm").alias("T_bm"),                # pooled B vs pooled M
            "pbs_target", "pbs_ref",
            "pbsn1_target", "pbsn1_ref",
        ])
        frames.append(q)
    return frames

def main():
    args = parse_args()

    targets = [x.strip() for x in args.target_reps.split(",") if x.strip()]
    refs = [x.strip() for x in args.ref_reps.split(",") if x.strip()]
    outgroups = [x.strip() for x in args.outgroup_reps.split(",") if x.strip()]
    gens = [x.strip() for x in args.generations.split(",") if x.strip()]
    o_gen_arg = args.outgroup_gen
    founders = [x.strip() for x in args.founder_reps.split(",") if x.strip()]
    founder_gen = args.founder_gen.strip()

    try:
        lf = pl.scan_csv(args.input, null_values=["NA", "nan", "."])
        schema_cols = set(lf.collect_schema().names())
    except Exception as e:
        pl.DataFrame({c: pl.Series([], dtype=pl.Utf8) for c in [
            "chrom","start","end","gen","lineage","total_divergence",
            "T_tm","T_bm","pbs_target","pbs_ref","pbsn1_target","pbsn1_ref"
        ]}).write_csv(args.output)
        sys.exit(0)

    output_schema_cols = [
        "chrom","start","end","gen","lineage","total_divergence",
        "T_tm","T_bm", "pbs_target","pbs_ref","pbsn1_target","pbsn1_ref"
    ]

    if args.regions:
        region_filters = []
        for region in args.regions.split(';'):
            if not region.strip(): continue
            try:
                r_chrom, r_start, r_end = region.split(':')
                region_filters.append((pl.col("chrom") == r_chrom) & (pl.col("start") >= int(r_start)) & (pl.col("start") < int(r_end)))
            except: pass
        if region_filters:
            combined = region_filters[0]
            for f in region_filters[1:]: combined = combined | f
            lf = lf.filter(combined)

    if args.mode == "per_rep":
        if len(targets) != len(refs):
            sys.exit(f"per_rep mode requires equal-length --target-reps and --ref-reps (got {len(targets)} vs {len(refs)})")

        per_rep_combos = []

        for t_rep, b_rep in zip(targets, refs):
            lineage_id = t_rep  # e.g. "T1"

            for gen in gens:
                if founders and gen == founder_gen:
                    continue

                cur_o_gen = gen if o_gen_arg == "dynamic" else o_gen_arg

                # Ti's triangle columns
                ti_b_cols = [c for c in (get_fst_col_name(schema_cols, t_rep, gen, b, gen) for b in refs) if c]
                ti_m_cols = [c for c in (get_fst_col_name(schema_cols, t_rep, gen, m, cur_o_gen) for m in outgroups) if c]
                bm_cols   = [c for c in (get_fst_col_name(schema_cols, b, gen, m, cur_o_gen) for b in refs for m in outgroups) if c]

                # Bi's triangle columns
                bi_t_cols  = [c for c in (get_fst_col_name(schema_cols, b_rep, gen, t, gen) for t in targets) if c]
                bi_m_cols  = [c for c in (get_fst_col_name(schema_cols, b_rep, gen, m, cur_o_gen) for m in outgroups) if c]
                tm_cols    = [c for c in (get_fst_col_name(schema_cols, t, gen, m, cur_o_gen) for t in targets for m in outgroups) if c]

                if ti_b_cols and ti_m_cols and bm_cols and bi_t_cols and bi_m_cols and tm_cols:
                    per_rep_combos.append({
                        "lineage": lineage_id,
                        "gen": gen,
                        "ti_b_cols": ti_b_cols,
                        "ti_m_cols": ti_m_cols,
                        "bm_cols": bm_cols,
                        "bi_t_cols": bi_t_cols,
                        "bi_m_cols": bi_m_cols,
                        "tm_cols": tm_cols,
                    })

        if founders and founder_gen in gens:
            distinct_f = list(dict.fromkeys(founders))
            if len(distinct_f) >= 3:
                f_t, f_r, f_o = distinct_f[0], distinct_f[1], distinct_f[2]
            elif len(distinct_f) == 2:
                f_t, f_r, f_o = distinct_f[0], distinct_f[1], distinct_f[0]
            else:
                f_t = f_r = f_o = distinct_f[0]

            col_f1 = get_fst_col_name(schema_cols, f_t, founder_gen, f_r, founder_gen)
            col_f2 = get_fst_col_name(schema_cols, f_t, founder_gen, f_o, founder_gen)
            col_f3 = get_fst_col_name(schema_cols, f_r, founder_gen, f_o, founder_gen)
            if col_f1 and col_f2 and col_f3:
                existing_lineages = list(dict.fromkeys(c["lineage"] for c in per_rep_combos))
                for lin in existing_lineages:
                    per_rep_combos.append({
                        "lineage": lin,
                        "gen": founder_gen,
                        "ti_b_cols": [col_f1],
                        "ti_m_cols": [col_f2],
                        "bm_cols":   [col_f3],
                        "bi_t_cols": [col_f1],  # symmetric: same near-zero values
                        "bi_m_cols": [col_f2],
                        "tm_cols":   [col_f3],
                    })

        if not per_rep_combos:
            with open(args.output, "w") as f:
                f.write(",".join(output_schema_cols) + "\n")
            sys.exit(0)

        frames = _build_per_rep_frames(lf, per_rep_combos)

    else:
        valid_combinations = []

        for t_rep, r_rep, o_rep in itertools.product(targets, refs, outgroups):
            lineage_id = f"{t_rep}_{r_rep}_{o_rep}"

            for gen in gens:
                if founders and gen == founder_gen:
                    continue

                cur_o_gen = gen if o_gen_arg == "dynamic" else o_gen_arg

                col_TB = get_fst_col_name(schema_cols, t_rep, gen, r_rep, gen)
                col_TM = get_fst_col_name(schema_cols, t_rep, gen, o_rep, cur_o_gen)
                col_BM = get_fst_col_name(schema_cols, r_rep, gen, o_rep, cur_o_gen)

                if col_TB and col_TM and col_BM:
                    valid_combinations.append({
                        "lineage": lineage_id, "gen": gen,
                        "col_TB": col_TB, "col_TM": col_TM, "col_BM": col_BM
                    })

        if founders and founder_gen in gens:
            distinct_f = list(dict.fromkeys(founders))
            if len(distinct_f) >= 3:
                f_t, f_r, f_o = distinct_f[0], distinct_f[1], distinct_f[2]
            elif len(distinct_f) == 2:
                f_t, f_r, f_o = distinct_f[0], distinct_f[1], distinct_f[0]
            else:
                f_t = f_r = f_o = distinct_f[0]
            col_TB_f = get_fst_col_name(schema_cols, f_t, founder_gen, f_r, founder_gen)
            col_TM_f = get_fst_col_name(schema_cols, f_t, founder_gen, f_o, founder_gen)
            col_BM_f = get_fst_col_name(schema_cols, f_r, founder_gen, f_o, founder_gen)
            if col_TB_f and col_TM_f and col_BM_f:
                existing_lineages = set(c["lineage"] for c in valid_combinations)
                for lin in existing_lineages:
                    valid_combinations.append({
                        "lineage": lin, "gen": founder_gen,
                        "col_TB": col_TB_f, "col_TM": col_TM_f, "col_BM": col_BM_f
                    })

        if not valid_combinations:
            with open(args.output, "w") as f:
                f.write(",".join(output_schema_cols) + "\n")
            sys.exit(0)

        frames = _build_combinatorial_frames(lf, valid_combinations, output_schema_cols)

    if frames:
        pl.concat(frames).collect().write_csv(args.output)
    else:
        with open(args.output, "w") as f:
            f.write(",".join(output_schema_cols) + "\n")

if __name__ == "__main__":
    main()
