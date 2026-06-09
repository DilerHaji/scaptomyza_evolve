#!/usr/bin/env python3
import argparse
import polars as pl
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
import matplotlib.font_manager as font_manager
from scipy.ndimage import gaussian_filter1d
import sys
import os
import copy
import json

font_dirs = ['/usr/share/fonts/truetype/msttcorefonts', '/usr/share/fonts/dejavu']
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
for font_file in font_files:
    try:
        font_manager.fontManager.addfont(font_file)
    except:
        pass

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rcParams['xtick.major.width'] = 0.5
matplotlib.rcParams['ytick.major.width'] = 0.5
matplotlib.rcParams['axes.spines.top'] = False
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.fonttype'] = 42

def parse_window_size(filepath):
    folder_name = os.path.basename(os.path.dirname(filepath))
    parts = folder_name.split('_')
    ints = [int(p) for p in parts if p.isdigit()]
    if ints: return max(ints)
    return 0

def parse_fasta_lengths(fasta_path):
    if not os.path.exists(fasta_path): sys.exit(f"Error: FASTA not found: {fasta_path}")
    sizes = {}
    total_len = 0
    current_chrom = None
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    current_chrom = line[1:].split()[0]
                    sizes[current_chrom] = 0
                else:
                    if current_chrom:
                        sizes[current_chrom] += len(line)
                        total_len += len(line)
    except Exception as e: sys.exit(f"Error parsing FASTA: {e}")
    return sizes, total_len

def load_chrom_map(map_file):
    if not map_file: return None
    mapping = {}
    try:
        with open(map_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2: mapping[parts[0]] = parts[1]
    except: pass
    return mapping

def mask_gaps(data_array, layout_items, bin_size):
    if data_array is None: return None
    n_bins = len(data_array)
    valid_mask = np.zeros(n_bins, dtype=bool)
    for item in layout_items:
        s_bin = max(0, min(int(item['start'] / bin_size), n_bins-1))
        e_bin = max(0, min(int((item['start'] + item['width']) / bin_size), n_bins-1))
        valid_mask[s_bin:e_bin+1] = True
    masked_data = data_array.copy()
    masked_data[~valid_mask] = np.nan
    return masked_data

def load_bed_annotations(bed_file, chrom_order, offsets):
    if not bed_file: return [], []
    annot_data, annot_labels = [], []
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith("#") or not line.strip(): continue
                parts = line.strip().split()
                if len(parts) < 3: continue
                chrom = parts[0]
                if chrom not in chrom_order: continue
                try:
                    start, end = int(parts[1]), int(parts[2])
                    label = parts[3] if len(parts) > 3 else ""
                    annot_data.append((start + offsets[chrom], end + offsets[chrom], label))
                    annot_labels.append(label)
                except: continue
    except: pass
    return annot_data, annot_labels


def load_pbs_manhattan(gen_agg_file, stat_col, gen1, gen2, chrom_order, offsets):
    try:
        df = pl.read_csv(gen_agg_file, null_values=['NA', 'nan', 'NaN', ''],
                         schema_overrides={'start': pl.Int64, 'end': pl.Int64,
                                           'gen': pl.Int64})
    except Exception as e:
        return np.array([]), np.array([]), np.array([]), np.array([])

    is_fst = (stat_col == 'fst_tb')
    if is_fst:
        if 'pbs_target' not in df.columns or 'pbs_ref' not in df.columns:
            return np.array([]), np.array([]), np.array([]), np.array([])
    elif stat_col not in df.columns:
        return np.array([]), np.array([]), np.array([]), np.array([])

    x1, y1, x2, y2 = [], [], [], []
    for gen_val, xl, yl in [(gen1, x1, y1), (gen2, x2, y2)]:
        df_g = df.filter(pl.col('gen') == gen_val).filter(
            pl.col('chrom').is_in(chrom_order))
        if df_g.height == 0:
            continue
        for chrom in chrom_order:
            c = df_g.filter(pl.col('chrom') == chrom)
            if c.height == 0: continue
            offset = offsets[chrom]
            mids = ((c['start'] + c['end']) / 2 + offset).to_numpy()
            if is_fst:
                t_tb = c['pbs_target'].to_numpy() + c['pbs_ref'].to_numpy()
                vals = 1.0 - np.exp(-t_tb)
            else:
                vals = c[stat_col].to_numpy()
            xl.append(mids)
            yl.append(vals)

    x1 = np.concatenate(x1) if x1 else np.array([])
    y1 = np.concatenate(y1) if y1 else np.array([])
    x2 = np.concatenate(x2) if x2 else np.array([])
    y2 = np.concatenate(y2) if y2 else np.array([])
    return x1, y1, x2, y2

def _remap_scatter(xs, ys, keep, cum_map, bin_size):
    if len(xs) == 0:
        return xs, ys
    n = len(keep)
    bin_idx = np.clip((xs / bin_size).astype(int), 0, n - 1)
    valid = keep[bin_idx]
    new_x = cum_map[bin_idx[valid]] * bin_size
    return new_x, ys[valid]

def build_genome_structure(heatmap_files, top_n, gap, genome_fasta=None):
    fasta_sizes, _ = (None, 0.0)
    if genome_fasta: fasta_sizes, _ = parse_fasta_lengths(genome_fasta)
    df = pl.read_csv(heatmap_files[0])
    chrom_lens = df.group_by("chrom").agg(pl.max("end").alias("length"))
    top_chroms = chrom_lens.sort("length", descending=True).head(top_n)
    ordered_chroms = top_chroms["chrom"].to_list()
    offsets, current_x, layout_items = {}, 0, []
    for chrom in ordered_chroms:
        L = fasta_sizes[chrom] if fasta_sizes and chrom in fasta_sizes else top_chroms.filter(pl.col("chrom")==chrom)["length"][0]
        offsets[chrom] = current_x
        parts = chrom.split("_")
        name = parts[2] if len(parts) > 2 and parts[1].startswith("Sc") else chrom
        layout_items.append({"start": current_x, "width": L, "mid": current_x + L/2, "label": f"{name}", "chrom": chrom})
        current_x += L + gap
    return ordered_chroms, offsets, layout_items, current_x, fasta_sizes

def process_heatmap_binned(files, col_primary, col_secondary, chrom_order, offsets, total_width, bin_size, use_percentile, min_win, max_win, layout_items):
    file_map = sorted([(parse_window_size(f), f) for f in files if min_win <= parse_window_size(f) <= max_win], key=lambda x: x[0])
    sizes, sorted_files = [x[0] for x in file_map], [x[1] for x in file_map]
    n_bins = int(total_width / bin_size) + 1

    matrix_1 = np.full((len(sizes), n_bins), np.nan)
    matrix_2 = np.full((len(sizes), n_bins), np.nan) if col_secondary else None

    for r_idx, fpath in enumerate(sorted_files):
        try:
            q = pl.scan_csv(fpath).filter(pl.col("chrom").is_in(chrom_order))
            df = q.collect()
            if df.height == 0: continue

            def fill_matrix_row(matrix, column_name):
                if column_name not in df.columns: return
                vals = df[column_name].to_numpy()
                valid_vals = vals[~np.isnan(vals)]
                if len(valid_vals) == 0: return
                if use_percentile:
                    df_c = df.with_columns((pl.col(column_name).rank() / pl.col(column_name).count()).alias("trans_val"))
                else:
                    mean, std = np.nanmean(valid_vals), np.nanstd(valid_vals)
                    df_c = df.with_columns(((pl.col(column_name) - mean) / (std + 1e-9)).alias("trans_val"))
                for chrom in chrom_order:
                    c_data = df_c.filter(pl.col("chrom") == chrom)
                    if c_data.height == 0: continue
                    offset = offsets[chrom]
                    c_data = c_data.with_columns([(((pl.col("start") + pl.col("end"))/2 + offset)/bin_size).floor().cast(pl.Int64).alias("bin_idx")])
                    agg = c_data.group_by("bin_idx").agg(pl.col("trans_val").mean())
                    idx, val = agg["bin_idx"].to_numpy(), agg["trans_val"].to_numpy()
                    valid = (idx >= 0) & (idx < n_bins)
                    matrix[r_idx, idx[valid]] = val[valid]

            fill_matrix_row(matrix_1, col_primary)
            if col_secondary and matrix_2 is not None:
                fill_matrix_row(matrix_2, col_secondary)
        except Exception as e: pass

    return matrix_1, matrix_2, sizes, sorted_files

def process_repeats_binned(repeat_file, chrom_map, chrom_order, offsets, total_width, bin_size, fasta_sizes):
    if not repeat_file: return None, None
    n_bins = int(total_width / bin_size) + 1
    coverage = np.zeros(n_bins, dtype=np.float32)
    try:
        with open(repeat_file, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) < 7 or not parts[0][0].isdigit(): continue
                chrom = chrom_map.get(parts[4], parts[4]) if chrom_map else parts[4]
                if chrom not in chrom_order: continue
                try:
                    s, e = int(parts[5]), int(parts[6])
                    g_s, g_e = min(s, e) + offsets[chrom], max(s, e) + offsets[chrom]
                    b_s, b_e = int(g_s // bin_size), int(g_e // bin_size)
                    for b in range(b_s, min(b_e + 1, n_bins)):
                        coverage[b] += max(0, min(g_e, (b+1)*bin_size) - max(g_s, b*bin_size))
                except: continue
    except: pass
    return np.arange(n_bins) * bin_size, np.clip(coverage / bin_size, 0, 1)

def process_glm_smoothed_mean(glm_file, pval_col, chrom_order, offsets, total_width, bin_size, min_maf, repeat_densities, layout_items, smooth_sigma=1.0):
    if not glm_file: return None, None, None, None, None, 0, (0,1), None

    q = pl.scan_csv(glm_file)
    schema = q.collect_schema().names()

    if "converged" in schema: q = q.filter(pl.col("converged") == True)
    if "average_freq" in schema: q = q.filter((pl.col("average_freq") >= min_maf) & (pl.col("average_freq") <= (1-min_maf)))

    eff_col = f"coef_{pval_col[2:]}" if pval_col.startswith("p_") else "beta"
    beta_expr = pl.col(eff_col).fill_null(0.0) if eff_col in schema else pl.lit(0.0)
    af_expr = pl.col("average_freq") if "average_freq" in schema else pl.lit(0.5)

    df = q.select([
        pl.col("chrom"), pl.col("pos"),
        pl.when(pl.col(pval_col) <= 0).then(1e-300).otherwise(pl.col(pval_col)).alias("p"),
        beta_expr.alias("beta"), af_expr.alias("af")
    ]).filter(pl.col("chrom").is_in(chrom_order)).with_columns([
        (-pl.col("p").log10()).alias("log_p"),
        (pl.col("beta") * pl.col("af") * (1 - pl.col("af")) * 100).abs().alias("effect_size")
    ]).collect()

    n_bins = int(total_width / bin_size) + 1
    if repeat_densities is not None: rep_df = pl.DataFrame({"bin_idx": np.arange(len(repeat_densities)), "rep_d": repeat_densities})
    else: rep_df = pl.DataFrame({"bin_idx": np.arange(n_bins), "rep_d": np.zeros(n_bins)})

    binned_data = []
    for chrom in chrom_order:
        c_data = df.filter(pl.col("chrom") == chrom)
        if c_data.height == 0: continue
        offset = offsets[chrom]
        c_data = c_data.with_columns([(((pl.col("pos") + offset) / bin_size).floor().cast(pl.Int64)).alias("bin_idx")])
        agg = c_data.group_by("bin_idx").agg([
            pl.col("log_p").mean().alias("mean_logp"),
            pl.col("log_p").std().fill_null(0).alias("std_logp"),
            pl.col("log_p").count().alias("count_logp"),
            pl.col("effect_size").mean().alias("mean_effect"),
            pl.col("log_p").median().alias("median_logp"),
            pl.col("log_p").quantile(0.25).alias("lower_q_logp"),
            pl.col("log_p").quantile(0.75).alias("upper_q_logp")
        ])
        binned_data.append(agg)

    if not binned_data: return np.arange(n_bins)*bin_size, np.zeros(n_bins), np.zeros(n_bins), np.zeros(n_bins), df, 0, (0,1), None

    full_agg = pl.concat(binned_data)
    full_agg = full_agg.join(rep_df, on="bin_idx", how="left")

    clean_bins = full_agg.filter(pl.col("rep_d") < 0.2)
    clean_bg_val = clean_bins["median_logp"].mean() if clean_bins.height > 0 else full_agg["median_logp"].mean()

    x_arr = np.arange(n_bins) * bin_size
    p_raw = np.full(n_bins, np.nan)
    sem_raw = np.full(n_bins, np.nan)
    e_raw = np.full(n_bins, np.nan)
    p_med = np.full(n_bins, np.nan)

    idx = full_agg["bin_idx"].to_numpy()
    valid = (idx >= 0) & (idx < n_bins)
    v_idx = idx[valid]

    p_raw[v_idx] = full_agg["mean_logp"].to_numpy()[valid]
    sem_raw[v_idx] = full_agg["std_logp"].to_numpy()[valid] / np.sqrt(full_agg["count_logp"].to_numpy()[valid])
    e_raw[v_idx] = full_agg["mean_effect"].to_numpy()[valid]
    p_med[v_idx] = full_agg["median_logp"].to_numpy()[valid]

    p_raw[np.isnan(p_raw)] = clean_bg_val
    p_smooth = gaussian_filter1d(p_raw, sigma=smooth_sigma)

    dist_stats = (mask_gaps(p_med, layout_items, bin_size), None, None)

    return x_arr, mask_gaps(p_smooth, layout_items, bin_size), mask_gaps(sem_raw, layout_items, bin_size), mask_gaps(e_raw, layout_items, bin_size), df, clean_bg_val, (0,1), dist_stats

def process_overlaps_and_candidates(df_glm, heatmap_files, heatmap_col, chrom_order, offsets, total_width, bin_size, fst_pct, glm_pct, overlap_size, layout_items):
    n_bins = int(total_width / bin_size) + 1
    if not heatmap_files or df_glm is None: return np.zeros(n_bins), np.zeros(n_bins), []

    file_map = [(parse_window_size(f), f) for f in heatmap_files]
    target_fst_file = min(file_map, key=lambda x: abs(x[0] - overlap_size))[1]

    df_fst = pl.read_csv(target_fst_file).filter(pl.col("chrom").is_in(chrom_order))

    if heatmap_col not in df_fst.columns:
        return np.zeros(n_bins), np.zeros(n_bins), []

    valid_fst = df_fst[heatmap_col].drop_nans()
    if len(valid_fst) == 0: return np.zeros(n_bins), np.zeros(n_bins), []

    fst_thresh = valid_fst.quantile((100 - fst_pct) / 100.0)
    glm_thresh = df_glm["log_p"].quantile((100 - glm_pct) / 100.0)

    hits_global = []
    for chrom in chrom_order:
        c_fst = df_fst.filter((pl.col("chrom") == chrom) & (pl.col(heatmap_col) >= fst_thresh))
        c_glm = df_glm.filter((pl.col("chrom") == chrom) & (pl.col("log_p") >= glm_thresh))

        if c_fst.height == 0 or c_glm.height == 0: continue

        high_glm_pos = c_glm["pos"].to_numpy()
        starts, ends = c_fst["start"].to_numpy(), c_fst["end"].to_numpy()

        idx = np.searchsorted(ends, high_glm_pos, side='right')
        valid_idx_mask = idx < len(ends)
        potential = valid_idx_mask & (high_glm_pos >= starts[np.clip(idx, 0, len(starts)-1)])

        if np.any(potential):
            hits_global.append(high_glm_pos[potential] + offsets[chrom])

    if not hits_global: return np.arange(n_bins) * bin_size, np.zeros(n_bins), []

    all_hits = np.concatenate(hits_global)
    counts, _ = np.histogram(all_hits, bins=n_bins, range=(0, n_bins*bin_size))
    counts = np.nan_to_num(mask_gaps(counts.astype(float), layout_items, bin_size))

    valid_counts = counts[counts > 0]
    merged_candidates = []

    if len(valid_counts) > 0:
        count_thresh = np.percentile(valid_counts, 90)
        candidate_indices = np.where(counts >= count_thresh)[0]
        if len(candidate_indices) > 0:
            curr_cluster = [candidate_indices[0]]
            for i in candidate_indices[1:]:
                if i <= curr_cluster[-1] + 2: curr_cluster.append(i)
                else:
                    merged_candidates.append((curr_cluster[0] * bin_size, (curr_cluster[-1] + 1) * bin_size))
                    curr_cluster = [i]
            merged_candidates.append((curr_cluster[0] * bin_size, (curr_cluster[-1] + 1) * bin_size))

    return np.arange(n_bins) * bin_size, counts, merged_candidates

def squeeze_dataset(hm_matrix, hm_matrix2, glm_p, glm_sem, glm_e, over_y, rep_y,
                    annot_data, highlights, layout_items, bin_size, threshold,
                    smooth_sigma, dist_stats,
                    pbs_scatter=None):
    if rep_y is None:
        return (hm_matrix, hm_matrix2, glm_p, glm_sem, glm_e, over_y, rep_y,
                annot_data, highlights, layout_items, len(glm_p)*bin_size, dist_stats,
                pbs_scatter)

    n_bins = len(rep_y)
    smooth_rep = np.zeros_like(rep_y)
    for item in layout_items:
        s = max(0, min(int(item['start']/bin_size), n_bins))
        e = max(0, min(int((item['start']+item['width'])/bin_size), n_bins))
        if e > s: smooth_rep[s:e] = gaussian_filter1d(np.nan_to_num(rep_y[s:e]), sigma=smooth_sigma)

    keep = (smooth_rep <= threshold)
    cum_map = np.cumsum(keep) - 1

    new_layout = []
    cur_x = 0
    chrom_map_arr = np.full(n_bins, -1)
    for i, item in enumerate(layout_items):
        s = max(0, min(int(item['start']/bin_size), n_bins))
        e = max(0, min(int((item['start']+item['width'])/bin_size), n_bins))
        chrom_map_arr[s:e] = i

    kept_ids = chrom_map_arr[keep]
    for i in range(len(layout_items)):
        w = np.sum(kept_ids == i) * bin_size
        new_layout.append({"start": cur_x, "width": w, "mid": cur_x+w/2,
                            "label": layout_items[i]['label'], "chrom": layout_items[i]['chrom']})
        cur_x += w

    def remap(intervals):
        res = []
        for item in intervals:
            is_bed = len(item) == 3
            g_s, g_e = item[0], item[1]
            b_s = max(0, min(int(g_s/bin_size), n_bins-1))
            b_e = max(0, min(int(g_e/bin_size), n_bins-1))
            if keep[b_s] or keep[b_e] or (cum_map[b_e] > cum_map[b_s]):
                n_s = max(0, cum_map[b_s] * bin_size)
                n_e = (cum_map[b_e]+1) * bin_size
                res.append((n_s, n_e, item[2]) if is_bed else (n_s, n_e))
        return res

    new_hm  = hm_matrix[:, keep]  if hm_matrix  is not None else None
    new_hm2 = hm_matrix2[:, keep] if hm_matrix2 is not None else None
    new_dist = (dist_stats[0][keep], None, None) if dist_stats else None

    # Remap PBS scatter points
    new_pbs_scatter = None
    if pbs_scatter is not None:
        new_pbs_scatter = [_remap_scatter(sx, sy, keep, cum_map, bin_size)
                           for sx, sy in pbs_scatter]

    return (new_hm, new_hm2,
            glm_p[keep] if glm_p is not None else None,
            glm_sem[keep] if glm_sem is not None else None,
            glm_e[keep]   if glm_e  is not None else None,
            over_y[keep]  if over_y is not None else None,
            rep_y[keep],
            remap(annot_data), remap(highlights), new_layout, np.sum(keep)*bin_size,
            new_dist, new_pbs_scatter)

def export_candidates_tsv(df_glm, chrom_order, output_base, filtered_candidates, offsets):
    summary_file  = f"{output_base}_candidate_sites.tsv"
    detailed_file = f"{output_base}_candidate_sites_detailed.tsv"

    summary_schema = ["Region_ID", "Chrom", "Start", "End", "Max_LogP", "Mean_Pct_AF_Change", "SNP_Count"]
    detail_schema  = ["Region_ID", "Chrom", "Pos", "LogP", "Pct_AF_Change"]

    if not filtered_candidates or df_glm is None:
        pl.DataFrame(schema=summary_schema).write_csv(summary_file, separator="\t")
        pl.DataFrame(schema=detail_schema).write_csv(detailed_file, separator="\t")
        return

    candidates, site_details = [], []
    for i, (g_start, g_end) in enumerate(filtered_candidates):
        target_chrom, local_start, local_end = None, 0, 0
        for chrom in chrom_order:
            if g_start >= offsets[chrom]:
                target_chrom = chrom
                local_start  = g_start - offsets[chrom]
                local_end    = g_end   - offsets[chrom]
            else: break
        if not target_chrom: continue
        c_glm = df_glm.filter((pl.col("chrom") == target_chrom) &
                               (pl.col("pos") >= local_start) &
                               (pl.col("pos") <= local_end))
        if c_glm.height == 0: continue
        snps_p   = c_glm["log_p"].to_numpy()
        snps_e   = c_glm["effect_size"].to_numpy()
        snps_pos = c_glm["pos"].to_numpy()
        reg_id = f"Region_{i+1}"
        candidates.append({
            "Region_ID": reg_id, "Chrom": target_chrom,
            "Start": int(local_start), "End": int(local_end),
            "Max_LogP": np.max(snps_p), "Mean_Pct_AF_Change": np.mean(snps_e),
            "SNP_Count": len(snps_pos)
        })
        site_details.append(pl.DataFrame({
            "Region_ID": [reg_id]*len(snps_pos),
            "Chrom": [target_chrom]*len(snps_pos),
            "Pos": snps_pos, "LogP": snps_p, "Pct_AF_Change": snps_e
        }))

    if candidates:    pl.DataFrame(candidates).write_csv(summary_file, separator="\t")
    if site_details:  pl.concat(site_details).write_csv(detailed_file, separator="\t")

def export_legends(args, output_base, c_line, c_ribbon, c_repeat, c_overlap, c_gene, cmap):
    fig_leg_hm, ax_leg_hm = plt.subplots(figsize=(2.2, 1.6))
    C_T_rgb        = np.array([0xED, 0xB7, 0x2D], dtype=float) / 255
    C_B_rgb        = np.array([0x49, 0x9F, 0xFF], dtype=float) / 255
    C_parallel_rgb = np.array([0xCC, 0x33, 0x33], dtype=float) / 255
    n = 128
    dir_thr   = getattr(args, 'heatmap_dir_threshold', 0.25)
    mag_gamma = getattr(args, 'heatmap_mag_gamma', 2.0)
    dirs = np.linspace(-1, 1, n)
    mags = np.linspace(0, 1, n)
    legend_img = np.ones((n, n, 3), dtype=float)
    for yi, mag in enumerate(mags):
        mag_d = mag ** mag_gamma
        for xi, dirn in enumerate(dirs):
            if dirn > dir_thr:
                base = C_T_rgb
            elif dirn < -dir_thr:
                base = C_B_rgb
            else:
                base = C_parallel_rgb
            legend_img[yi, xi] = np.clip(mag_d * base + (1 - mag_d), 0, 1)
    ax_leg_hm.imshow(legend_img, origin='lower', aspect='auto',
                     extent=[-1, 1, 0, 1], interpolation='bilinear')
    ax_leg_hm.axvline( dir_thr, color='white', lw=0.8, linestyle='--', alpha=0.7)
    ax_leg_hm.axvline(-dir_thr, color='white', lw=0.8, linestyle='--', alpha=0.7)
    ax_leg_hm.set_xlabel("Direction  (B ← → T)", fontsize=6)
    ax_leg_hm.set_ylabel("Magnitude\n(total divergence)", fontsize=6)
    ax_leg_hm.tick_params(labelsize=5)
    for x, label, col in [(-0.85, 'B-specific', '#499FFF'),
                           (0.0,  'parallel',   '#CC3333'),
                           (0.85, 'T-specific', '#EDB72D')]:
        ax_leg_hm.text(x, 0.92, label, ha='center', va='top', fontsize=5,
                       color=col, fontweight='bold', transform=ax_leg_hm.transData)
    ax_leg_hm.text(0.0, 0.05, 'nothing happening', ha='center', fontsize=5,
                   color='gray', transform=ax_leg_hm.transData)
    plt.tight_layout()
    plt.savefig(f"{output_base}_legend_heatmap.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_base}_legend_heatmap.svg", bbox_inches='tight')
    plt.close(fig_leg_hm)

    fig_leg_tr = plt.figure(figsize=(3, 2))
    ax_leg_tr = fig_leg_tr.add_subplot(111)
    ax_leg_tr.axis('off')
    legend_elements = [
        Line2D([0], [0], color=c_line,  lw=1,   label='GLM Mean'),
        Patch(facecolor=c_ribbon, alpha=args.ribbon_alpha, label=f'GLM ± {args.sem_scale:g} SEM'),
        Patch(facecolor=c_repeat,  label='Repeat Density'),
        Patch(facecolor="#0072B2", label='Effect Size (AME)'),
        Patch(facecolor=c_overlap, label='Overlap Count'),
        Patch(facecolor=c_gene,    label='Genes'),
        Line2D([0], [0], color=args.pbs_gen1_color, lw=0, marker='o', markersize=4,
               label=f'PBS G{args.pbs_gen1:02d}'),
        Line2D([0], [0], color=args.pbs_gen2_color, lw=0, marker='o', markersize=4,
               label=f'PBS G{args.pbs_gen2:02d}'),
    ]
    ax_leg_tr.legend(handles=legend_elements, loc='center', frameon=False, fontsize=7)
    plt.savefig(f"{output_base}_legend_tracks.png", bbox_inches='tight', dpi=300)
    plt.close(fig_leg_tr)


def render_final(hm_matrix1, hm_matrix2, win_sizes,
                 glm_x, glm_p, glm_sem, glm_e, gl_bg_mean, gl_bg_range, dist_stats,
                 over_x, over_y, highlight_intervals,
                 rep_x, rep_y,
                 annot_data, annot_labels,
                 layout_items, total_width, bin_size, args, filename,
                 pbs_scatter=None):

    FIG_W, FIG_H = args.fig_width, args.fig_height
    FS_STD = 7 * args.text_scale

    C_GLM_LINE = "#212121"
    C_GLM_RIBBON = "#BDBDBD"
    C_REPEAT   = "#555555"
    C_OVERLAP  = "#009E73"
    C_GLM_LAWN = "#0072B2"
    C_GENES    = "#424242"

    cmap_T = mcolors.LinearSegmentedColormap.from_list(
        'T_seq', ['#F8ECC0', '#EDB72D', '#5C3800'], N=256)
    cmap_B = mcolors.LinearSegmentedColormap.from_list(
        'B_seq', ['#C0D8FF', '#499FFF', '#002060'], N=256)
    cmap_T.set_bad(color='#e8e8e8')
    cmap_B.set_bad(color='#e8e8e8')

    if getattr(args, 'heatmap_diff', False):
        cmap_heat = mcolors.LinearSegmentedColormap.from_list(
            'BT_div', ['#499FFF', '#BB4C78', '#EDB72D'], N=256)
        cmap_heat.set_bad(color='white')
    elif args.percentile:
        cmap_heat = copy.copy(plt.get_cmap("viridis"))
        cmap_heat.set_bad(color='white')
    else:
        cmap_heat = copy.copy(plt.get_cmap("coolwarm"))
        cmap_heat.set_bad(color='white')

    active_tracks = args.tracks
    ratios = []
    for t in active_tracks:
        if t == 'heatmap':
            ratios.append(1.2)
        elif t == 'glm':          ratios.append(0.8)
        elif t == 'pbs_manhattan': ratios.append(0.5)
        elif t == 'footer':        ratios.append(0.1)
        else:                      ratios.append(0.3)

    if args.fig_ratios: ratios = args.fig_ratios

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)
    gs = GridSpec(len(active_tracks), 1, height_ratios=ratios,
                  hspace=0.1, left=0.1, right=0.98, top=0.95, bottom=0.05)

    axes_dict = {}
    share_ax = None
    for i, track_name in enumerate(active_tracks):
        if track_name == 'footer':
            ax = fig.add_subplot(gs[i])
        else:
            if share_ax is None: ax = fig.add_subplot(gs[i]); share_ax = ax
            else:                ax = fig.add_subplot(gs[i], sharex=share_ax)
        axes_dict[track_name] = ax

    alphas = np.ones(int(total_width/bin_size)+1, dtype=float)
    if rep_y is not None:
        n_points = min(len(alphas), len(rep_y))
        slope = 1.0 - args.min_alpha
        alphas[:n_points] = 1.0 - (np.clip(np.nan_to_num(rep_y[:n_points]), 0, 1) * slope)

    if 'heatmap' in axes_dict and hm_matrix1 is not None:
        ax = axes_dict['heatmap']
        clip_hi = getattr(args, 'heatmap_clip_hi', 98)
        clip_lo = 100 - clip_hi

        def _norm_rows(mat, lo_pct, hi_pct):
            out = np.full_like(mat, np.nan, dtype=float)
            for i in range(mat.shape[0]):
                row   = mat[i]
                valid = row[np.isfinite(row)]
                if len(valid) == 0:
                    continue
                lo = np.nanpercentile(valid, lo_pct)
                hi = np.nanpercentile(valid, hi_pct)
                if hi > lo:
                    out[i] = np.clip((row - lo) / (hi - lo), 0, 1)
                else:
                    out[i] = np.where(np.isfinite(row), 0.5, np.nan)
            return out

        t_norm = _norm_rows(hm_matrix1, clip_lo, clip_hi)
        if hm_matrix2 is not None:
            b_norm = _norm_rows(hm_matrix2, clip_lo, clip_hi)
        else:
            b_norm = np.full_like(t_norm, 0.5)

        C_T_rgb        = np.array([0xED, 0xB7, 0x2D], dtype=float) / 255  # gold
        C_B_rgb        = np.array([0x49, 0x9F, 0xFF], dtype=float) / 255  # blue
        C_parallel_rgb = np.array([0xCC, 0x33, 0x33], dtype=float) / 255  # red

        dir_thr = getattr(args, 'heatmap_dir_threshold', 0.25)
        mag_gamma = getattr(args, 'heatmap_mag_gamma', 2.0)

        nan_mask = ~(np.isfinite(t_norm) & np.isfinite(b_norm))
        t_safe = np.where(np.isfinite(t_norm), t_norm, 0.0)
        b_safe = np.where(np.isfinite(b_norm), b_norm, 0.0)

        direction = t_safe - b_safe                          # [-1, +1]
        magnitude = (t_safe + b_safe) / 2.0                 # [0, 1]

        is_T = direction >  dir_thr
        is_B = direction < -dir_thr
        is_P = ~(is_T | is_B)
        base = (is_T[..., None] * C_T_rgb +
                is_B[..., None] * C_B_rgb +
                is_P[..., None] * C_parallel_rgb)

        mag_display = np.power(np.clip(magnitude, 0.0, 1.0), mag_gamma)

        final = mag_display[..., None] * base + (1.0 - mag_display[..., None])

        rgb_image        = np.ones((*t_norm.shape, 4), dtype=float)
        rgb_image[..., :3] = np.clip(final, 0, 1)
        rgb_image[..., 3]  = 1.0
        rgb_image[nan_mask] = [0.91, 0.91, 0.91, 1.0]       # NaN → grey

        ax.imshow(rgb_image, aspect='auto', origin='lower',
                  extent=[0, total_width, 0, len(win_sizes)],
                  interpolation='nearest')

        yticks = np.linspace(0, len(win_sizes) - 1, 3)
        ax.set_yticks(yticks)
        ax.set_yticklabels([f"{int(win_sizes[int(i)] / 1000)}" for i in yticks])
        ax.set_ylabel("Window (kbp)\nT − B", rotation=0, ha='right',
                      va='center', fontsize=FS_STD)

        for item in layout_items:
            ax.add_patch(Rectangle((item["start"], 0), item["width"], len(win_sizes),
                                   linewidth=0.5, edgecolor='black',
                                   facecolor='none', zorder=5))

        ax.tick_params(bottom=False, labelbottom=False)

    if 'glm' in axes_dict and glm_p is not None:
        ax = axes_dict['glm']
        if args.plot_bin_stats and dist_stats:
            d_med, _, _ = dist_stats
            mask = ~np.isnan(d_med)
            if np.any(mask):
                ax.scatter(glm_x[mask], d_med[mask], s=4, c='#222222', alpha=1.0,
                           zorder=5, linewidths=0, rasterized=True)
        if glm_sem is not None and not args.no_sem:
            scaled_sem = glm_sem * args.sem_scale
            rib_h = np.nan_to_num(2 * scaled_sem, nan=0.0)
            rib_b = np.nan_to_num(glm_p - scaled_sem, nan=0.0)
            rib_c = np.zeros((len(glm_x), 4))
            rib_c[:] = mcolors.to_rgba(C_GLM_RIBBON)
            rib_c[:, 3] = args.ribbon_alpha
            ax.bar(glm_x, rib_h, bottom=rib_b, width=bin_size, align='edge',
                   color=rib_c, linewidth=0, rasterized=True, zorder=2)
        if not args.no_trend_line:
            points   = np.array([glm_x, glm_p]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc_colors = np.zeros((len(segments), 4))
            lc_colors[:] = mcolors.to_rgba(C_GLM_LINE)
            lc_colors[:, 3] = alphas[:len(segments)]
            ax.add_collection(LineCollection(segments, colors=lc_colors, linewidth=0.8, zorder=3))
        ax.axhline(y=gl_bg_mean, color='#999999', linestyle=':', linewidth=0.5, zorder=1)
        ax.set_ylabel(r"Avg $-\log_{10}P$", rotation=0, ha='right', va='center', fontsize=FS_STD)
        ax.tick_params(bottom=False, labelbottom=False)
        ax.set_ylim(gl_bg_range[0] * 0.9, gl_bg_range[1] * 1.2)


    if 'pbs_manhattan' in axes_dict and pbs_scatter is not None:
        ax = axes_dict['pbs_manhattan']
        (px1, py1), (px2, py2) = pbs_scatter

        if len(px1) > 0:
            mask1 = ~np.isnan(py1)
            ax.scatter(px1[mask1], py1[mask1],
                       s=args.pbs_point_size, c=args.pbs_gen1_color,
                       alpha=args.pbs_alpha, linewidths=0,
                       rasterized=True, zorder=2, label=f'G{args.pbs_gen1:02d}')

        if len(px2) > 0:
            mask2 = ~np.isnan(py2)
            ax.scatter(px2[mask2], py2[mask2],
                       s=args.pbs_point_size, c=args.pbs_gen2_color,
                       alpha=args.pbs_alpha, linewidths=0,
                       rasterized=True, zorder=3, label=f'G{args.pbs_gen2:02d}')

        ax.axhline(0, color='#cccccc', linewidth=0.4, zorder=1)

        if args.pbs_stat == 'fst_tb':
            stat_label = "FST (B–T)"
        else:
            stat_label = args.pbs_stat.replace('_target', '').replace('_ref', '').upper()
        ax.set_ylabel(f"{stat_label}\n(G{args.pbs_gen1:02d}/G{args.pbs_gen2:02d})",
                      rotation=0, ha='right', va='center', fontsize=FS_STD)
        ax.legend(fontsize=5, frameon=False, loc='upper right',
                  markerscale=1.5, handletextpad=0.3)
        ax.tick_params(bottom=False, labelbottom=False)


    if 'overlap' in axes_dict and over_y is not None:
        ax = axes_dict['overlap']
        if np.max(np.nan_to_num(over_y)) > 0:
            over_colors = np.zeros((len(over_x), 4))
            over_colors[:] = mcolors.to_rgba(C_OVERLAP)
            over_colors[:, 3] = alphas[:len(over_x)]
            ax.bar(over_x, np.nan_to_num(over_y), width=bin_size, align='edge',
                   color=over_colors, linewidth=0, zorder=1)
        ax.set_ylabel("Overlap\nCount", rotation=0, ha='right', va='center', fontsize=FS_STD)
        ax.tick_params(bottom=False, labelbottom=False)

    if 'repeats' in axes_dict and rep_y is not None:
        ax = axes_dict['repeats']
        rep_y_masked = mask_gaps(rep_y, layout_items, bin_size)
        rep_colors = np.zeros((len(rep_x), 4))
        rep_colors[:] = mcolors.to_rgba(C_REPEAT)
        rep_colors[:, 3] = alphas[:len(rep_x)]
        ax.bar(rep_x, rep_y_masked, width=bin_size, align='edge',
               color=rep_colors, linewidth=0, rasterized=True, zorder=1)
        ax.set_ylabel("Repeat\nDensity", rotation=0, ha='right', va='center', fontsize=FS_STD)
        ax.tick_params(bottom=False, labelbottom=False)
        ax.set_ylim(0, 1)

    if 'genes' in axes_dict:
        ax = axes_dict['genes']
        ax.axis('off')
        if annot_data:
            for (g_start, g_end, label) in annot_data:
                ax.add_patch(Rectangle(
                    (g_start, 0.15), max(g_end - g_start, total_width * 0.002), 0.7,
                    facecolor=C_GENES, edgecolor='none', zorder=1))
                if label and (g_end - g_start) / total_width > 0.003:
                    ax.text((g_start + g_end) / 2, 0.85, label,
                            ha='center', va='bottom', fontsize=4.5,
                            color=C_GENES, clip_on=True)
        ax.set_ylabel("Genes", rotation=0, ha='right', va='center', fontsize=FS_STD)

    if 'footer' in axes_dict:
        ax = axes_dict['footer']
        ax.axis('off')
        for i, item in enumerate(layout_items):
            if (item["width"] / total_width) > 0.005:
                ax.plot([item["start"], item["start"] + item["width"]],
                        [0.8, 0.8], color="black", linewidth=1.5)
                ax.text(item["mid"], 0.45 if i % 2 == 0 else 0.0,
                        item["label"], ha='center', va='top', fontsize=FS_STD)

    for ax in axes_dict.values(): ax.set_xlim(0, total_width)

    _ax_bounds = {}
    try:
        _ref_ax = axes_dict.get('heatmap', list(axes_dict.values())[0])
        _pos    = _ref_ax.get_position()
        _ax_bounds = {
            "axes_x0_frac": _pos.x0,
            "axes_x1_frac": _pos.x1,
        }
    except Exception:
        pass

    plt.savefig(filename, bbox_inches='tight', dpi=300)
    plt.close(fig)

    export_legends(args, filename.rsplit('.', 1)[0],
                   C_GLM_LINE, C_GLM_RIBBON, C_REPEAT, C_OVERLAP, C_GENES, cmap_heat)
    return _ax_bounds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--heatmap-files",  nargs='+', required=True)
    parser.add_argument("--output",         required=True)
    parser.add_argument("--heatmap-col",         default="slope_pbs")
    parser.add_argument("--heatmap-col-ref",     default=None)
    parser.add_argument("--heatmap-diff",         action="store_true",
                        help="Plot single diverging difference heatmap (T − B) "
                             "instead of two-row target/ref layout")
    parser.add_argument("--heatmap-diff-vlim",    type=float, default=0.3,
                        help="Symmetric color limit for diff heatmap (default: 0.3)")
    parser.add_argument("--glm",            required=False)
    parser.add_argument("--target-pval",    required=False)
    parser.add_argument("--genome-fasta",   type=str)
    parser.add_argument("--repeat-masker",  type=str)
    parser.add_argument("--chrom-map",      type=str)
    parser.add_argument("--bed-file",       type=str)
    parser.add_argument("--bin-size",       type=int,   default=100000)
    parser.add_argument("--top-n",          type=int,   default=10)
    parser.add_argument("--gap",            type=int,   default=1_000_000)
    parser.add_argument("--min-maf",        type=float, default=0.01)
    parser.add_argument("--min-window",     type=int,   default=20000)
    parser.add_argument("--max-window",     type=int,   default=500000)
    parser.add_argument("--percentile",     action="store_true")
    parser.add_argument("--fst-top-pct",    type=float, default=5.0)
    parser.add_argument("--glm-top-pct",    type=float, default=5.0)
    parser.add_argument("--overlap-target-size", type=int, default=250000)
    parser.add_argument("--y-max",          type=float, default=None)
    parser.add_argument("--y-min",          type=float, default=None)
    parser.add_argument("--heatmap-clip-hi", type=float, default=98,
                        help="Per-row upper percentile clip for heatmap normalisation "
                             "(default 98). Lower end = 100 - this value.")
    parser.add_argument("--heatmap-dir-threshold", type=float, default=0.25,
                        help="Normalised direction threshold for discrete hue assignment. "
                             "direction > thr → T (gold); < -thr → B (blue); else → parallel (red). "
                             "(default 0.25)")
    parser.add_argument("--heatmap-mag-gamma", type=float, default=2.0,
                        help="Power applied to magnitude before colour mixing. "
                             "> 1 compresses low values toward white so high-signal regions "
                             "stand out more strongly. (default 2.0)")
    parser.add_argument("--ribbon-alpha",   type=float, default=0.4)
    parser.add_argument("--smooth-sigma",   type=float, default=1.0)
    parser.add_argument("--sem-scale",      type=float, default=3.0)
    parser.add_argument("--min-alpha",      type=float, default=0.4)
    parser.add_argument("--text-scale",     type=float, default=1.0)
    parser.add_argument("--fig-width",      type=float, default=7.5)
    parser.add_argument("--fig-height",     type=float, default=7.0)
    parser.add_argument("--fig-ratios",     nargs='+',  type=float, default=[])
    parser.add_argument("--tracks",         nargs='+',
                        default=["heatmap", "glm", "pbs_manhattan", "overlap", "genes", "footer"],
                        choices=["heatmap", "glm", "pbs_manhattan", "effect",
                                 "overlap", "repeats", "genes", "footer"])
    parser.add_argument("--remove-repeats",       action="store_true")
    parser.add_argument("--repeat-threshold",     type=float, default=0.4)
    parser.add_argument("--repeat-filter-sigma",  type=float, default=5.0)
    parser.add_argument("--plot-bin-stats",        action="store_true")
    parser.add_argument("--bin-stats-alpha",       type=float, default=1.0)
    parser.add_argument("--stats-lower-q",         type=float, default=0.25)
    parser.add_argument("--stats-upper-q",         type=float, default=0.75)
    parser.add_argument("--no-sem",                action="store_true")
    parser.add_argument("--no-trend-line",         action="store_true")
    parser.add_argument("--connect-points",        action="store_true")

    parser.add_argument("--gen-agg-file",   type=str,   default=None,
                        help="Compact gen_agg CSV for PBS Manhattan track")
    parser.add_argument("--pbs-stat",       type=str,   default="pbs_target",
                        help="Column to plot in PBS Manhattan (e.g. pbs_target, pbsn1_target)")
    parser.add_argument("--pbs-gen1",       type=int,   default=1,
                        help="First generation (default: 1 = G01)")
    parser.add_argument("--pbs-gen2",       type=int,   default=10,
                        help="Second generation (default: 10 = G10)")
    parser.add_argument("--pbs-gen1-color", type=str,   default="#4575B4",
                        help="Color for gen1 scatter points")
    parser.add_argument("--pbs-gen2-color", type=str,   default="#D73027",
                        help="Color for gen2 scatter points")
    parser.add_argument("--pbs-alpha",      type=float, default=0.6,
                        help="Alpha for PBS scatter points")
    parser.add_argument("--pbs-point-size", type=float, default=3.0,
                        help="Marker size for PBS scatter points")

    args = parser.parse_args()

    chroms, offsets, layout, width, fasta_sizes = build_genome_structure(
        args.heatmap_files, args.top_n, args.gap, args.genome_fasta)

    hm_matrix1, hm_matrix2, sizes, sorted_files = process_heatmap_binned(
        args.heatmap_files, args.heatmap_col, args.heatmap_col_ref,
        chroms, offsets, width, args.bin_size, args.percentile,
        args.min_window, args.max_window, layout)

    chrom_map = load_chrom_map(args.chrom_map)
    rep_x, rep_y = None, None
    if args.repeat_masker and fasta_sizes:
        rep_x, rep_y = process_repeats_binned(
            args.repeat_masker, chrom_map, chroms, offsets, width, args.bin_size, fasta_sizes)

    glm_x, glm_p, glm_sem, glm_e, raw_glm, bg_mean, bg_range, dist_stats = (None,)*8
    over_x, over_y, highlights = (None, None, [])
    if args.glm:
        glm_x, glm_p, glm_sem, glm_e, raw_glm, bg_mean, bg_range, dist_stats = \
            process_glm_smoothed_mean(
                args.glm, args.target_pval, chroms, offsets, width,
                args.bin_size, args.min_maf, rep_y, layout, args.smooth_sigma)
        over_x, over_y, highlights = process_overlaps_and_candidates(
            raw_glm, sorted_files, args.heatmap_col, chroms, offsets, width,
            args.bin_size, args.fst_top_pct, args.glm_top_pct,
            args.overlap_target_size, layout)

    annot_data, annot_labels = load_bed_annotations(args.bed_file, chroms, offsets)

    pbs_scatter = None
    if args.gen_agg_file and 'pbs_manhattan' in args.tracks:
        px1, py1, px2, py2 = load_pbs_manhattan(
            args.gen_agg_file, args.pbs_stat,
            args.pbs_gen1, args.pbs_gen2, chroms, offsets)
        pbs_scatter = [(px1, py1), (px2, py2)]

    if args.remove_repeats and rep_y is not None:
        result = squeeze_dataset(
            hm_matrix1, hm_matrix2, glm_p, glm_sem, glm_e, over_y, rep_y,
            annot_data, highlights, layout, args.bin_size,
            args.repeat_threshold, args.repeat_filter_sigma, dist_stats,
            pbs_scatter=pbs_scatter)
        (hm_matrix1, hm_matrix2, glm_p, glm_sem, glm_e, over_y, rep_y,
         annot_data, highlights, layout, width, dist_stats, pbs_scatter) = result
        n_bins = len(rep_y)
        glm_x = rep_x = over_x = np.arange(n_bins) * args.bin_size

    _ax_bounds = render_final(hm_matrix1, hm_matrix2, sizes,
                 glm_x, glm_p, glm_sem, glm_e, bg_mean, bg_range, dist_stats,
                 over_x, over_y, highlights,
                 rep_x, rep_y, annot_data, annot_labels,
                 layout, width, args.bin_size, args, args.output,
                 pbs_scatter=pbs_scatter)

    import json as _json
    from PIL import Image as _PILImage
    _coords_path = args.output.rsplit('.', 1)[0] + '_coords.json'
    _img = _PILImage.open(args.output)
    _coords = {
        "total_width_bp": int(width),
        "image_width_px":  _img.size[0],
        "image_height_px": _img.size[1],
        "layout": [
            {"chrom": it["chrom"], "label": it["label"],
             "start_bp": int(it["start"]), "width_bp": int(it["width"])}
            for it in layout
        ],
    }
    _coords.update(_ax_bounds)
    with open(_coords_path, 'w') as _f:
        _json.dump(_coords, _f, indent=2)

    if args.glm:
        export_candidates_tsv(raw_glm, chroms,
                              os.path.splitext(args.output)[0], highlights, offsets)


if __name__ == "__main__":
    main()
