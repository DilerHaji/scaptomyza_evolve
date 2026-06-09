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
    if len(ints) >= 2: return ints[-2]
    elif len(ints) == 1: return ints[0]
    return 0

def parse_fasta_lengths(fasta_path):
    if not os.path.exists(fasta_path):
        sys.exit(f"Error: Genome FASTA file not found: {fasta_path}")
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
    except Exception as e:
        sys.exit(f"Error parsing FASTA file: {e}")
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
        s_bin = int(item['start'] / bin_size)
        e_bin = int((item['start'] + item['width']) / bin_size)
        s_bin = max(0, min(s_bin, n_bins-1))
        e_bin = max(0, min(e_bin, n_bins-1))
        valid_mask[s_bin:e_bin+1] = True
    masked_data = data_array.copy()
    masked_data[~valid_mask] = np.nan
    return masked_data

def load_bed_annotations(bed_file, chrom_order, offsets):
    if not bed_file: return [], []
    annot_data = []
    annot_labels = []
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
                    global_start = start + offsets[chrom]
                    global_end = end + offsets[chrom]
                    annot_data.append((global_start, global_end, label))
                    annot_labels.append(label)
                except ValueError: continue
    except Exception as e:
        return [], []
    return annot_data, annot_labels


def build_genome_structure(heatmap_files, top_n, gap, genome_fasta=None):
    fasta_sizes, _ = (None, 0.0)
    if genome_fasta:
        fasta_sizes, _ = parse_fasta_lengths(genome_fasta)
    
    df = pl.read_csv(heatmap_files[0])
    chrom_lens_csv = df.group_by("chrom").agg(pl.max("end").alias("length"))
    top_chroms = chrom_lens_csv.sort("length", descending=True).head(top_n)
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

def process_heatmap_binned(files, col, chrom_order, offsets, total_width, bin_size, use_percentile, min_win, max_win, layout_items):
    file_map = [(parse_window_size(f), f) for f in files if min_win <= parse_window_size(f) <= max_win]
    file_map.sort(key=lambda x: x[0])
    sorted_sizes, sorted_files = [x[0] for x in file_map], [x[1] for x in file_map]
    n_bins = int(total_width / bin_size) + 1
    matrix = np.full((len(sorted_sizes), n_bins), np.nan)

    for r_idx, fpath in enumerate(sorted_files):
        try:
            q = pl.scan_csv(fpath).filter(pl.col("chrom").is_in(chrom_order))
            df = q.collect()
            if df.height == 0: continue
            vals = df[col].to_numpy()
            if use_percentile:
                df = df.with_columns((pl.col(col).rank() / pl.col(col).count()).alias("trans_val"))
            else:
                df = df.with_columns(((pl.col(col) - np.nanmean(vals)) / (np.nanstd(vals) + 1e-9)).alias("trans_val"))
            for chrom in chrom_order:
                c_data = df.filter(pl.col("chrom") == chrom)
                if c_data.height == 0: continue
                chrom_offset = offsets[chrom]
                c_data = c_data.with_columns([(((pl.col("start") + pl.col("end")) / 2.0 + chrom_offset) / bin_size).floor().cast(pl.Int64).alias("bin_idx")])
                agg_data = c_data.group_by("bin_idx").agg(pl.col("trans_val").mean())
                b_indices = agg_data["bin_idx"].to_numpy()
                b_vals = agg_data["trans_val"].to_numpy()
                valid_mask = (b_indices >= 0) & (b_indices < n_bins)
                if np.any(valid_mask): matrix[r_idx, b_indices[valid_mask]] = b_vals[valid_mask]
        except Exception as e: pass
    return matrix, sorted_sizes, sorted_files

def process_repeats_binned(repeat_file, chrom_map, chrom_order, offsets, total_width, bin_size, fasta_sizes):
    if not repeat_file: return None, None
    n_bins = int(total_width / bin_size) + 1
    coverage = np.zeros(n_bins, dtype=np.float32)
    try:
        with open(repeat_file, 'r') as f:
            for line in f:
                line = line.strip() 
                if not line or not line[0].isdigit(): continue 
                parts = line.split()
                if len(parts) < 7: continue
                chrom = parts[4]
                if chrom_map and chrom in chrom_map: chrom = chrom_map[chrom]
                if chrom not in chrom_order: continue
                try:
                    s, e = int(parts[5]), int(parts[6])
                    g_start = min(s, e) + offsets[chrom]
                    g_end = max(s, e) + offsets[chrom]
                    bin_s, bin_e = int(g_start // bin_size), int(g_end // bin_size)
                    for b in range(bin_s, bin_e + 1):
                        if b >= n_bins: break
                        overlap = max(0, min(g_end, (b+1)*bin_size) - max(g_start, b*bin_size))
                        coverage[b] += overlap
                except: continue
    except: pass
    densities = np.clip(coverage / bin_size, 0, 1)
    return np.arange(n_bins) * bin_size, densities

def process_glm_smoothed_mean(glm_file, pval_col, chrom_order, offsets, total_width, bin_size, min_maf, repeat_densities, layout_items, smooth_sigma=1.0, lower_q=0.25, upper_q=0.75):
    q = pl.scan_csv(glm_file)
    schema = q.collect_schema().names()

    if "error" in schema:
        q = q.filter(pl.col("error") == "OK")

    if "converged" in schema:
        q = q.filter(pl.col("converged") == True)

    if "singular" in schema:
        q = q.filter(pl.col("singular") == False)

    if "average_freq" in schema:
        q = q.filter((pl.col("average_freq") >= min_maf) & (pl.col("average_freq") <= (1.0 - min_maf)))
        af_expr = pl.col("average_freq")
    else: 
        af_expr = pl.lit(0.5)

    eff_col = "coef_" + pval_col[2:] if pval_col.startswith("p_") else None
    if eff_col and eff_col not in schema: 
        eff_col = None
    
    beta_expr = pl.col(eff_col).fill_null(0.0) if eff_col else pl.lit(0.0)

    df = q.select([
        pl.col("chrom"), pl.col("pos"),
        pl.when(pl.col(pval_col) <= 0).then(1e-300).otherwise(pl.col(pval_col)).alias("p"),
        beta_expr.alias("beta"), af_expr.alias("af")
    ]).filter(pl.col("chrom").is_in(chrom_order)).with_columns([
        (-pl.col("p").log10()).alias("log_p"),
        (pl.col("beta") * pl.col("af") * (1 - pl.col("af")) * 100.0).abs().alias("effect_size")
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
            pl.col("log_p").mean().alias("mean_logp"), # Kept for SEM calc
            pl.col("log_p").std().fill_null(0).alias("std_logp"), # Kept for SEM calc
            pl.col("log_p").count().alias("count_logp"),
            pl.col("effect_size").mean().alias("mean_effect"),
            pl.col("log_p").median().alias("median_logp"),
            pl.col("log_p").quantile(lower_q).alias("lower_q_logp"),
            pl.col("log_p").quantile(upper_q).alias("upper_q_logp")
        ])
        binned_data.append(agg)

    if not binned_data:
        zeros = np.zeros(n_bins)
        return np.arange(n_bins)*bin_size, zeros, zeros, zeros, df, 0.0, (0.0, 1.0), (zeros, zeros, zeros)

    full_agg = pl.concat(binned_data)
    full_agg = full_agg.join(rep_df, on="bin_idx", how="left")
    
    clean_bins = full_agg.filter(pl.col("rep_d") < 0.2)
    clean_bg_val = clean_bins["median_logp"].mean() if clean_bins.height > 0 else full_agg["median_logp"].mean()
    clean_min = clean_bins["median_logp"].min() if clean_bins.height > 0 else full_agg["median_logp"].min()
    clean_max = clean_bins["median_logp"].max() if clean_bins.height > 0 else full_agg["median_logp"].max()

    x_arr = np.arange(n_bins) * bin_size
    
    p_main = np.full(n_bins, np.nan)
    sem_raw = np.full(n_bins, np.nan)
    e_raw = np.full(n_bins, np.nan)
    p_med = np.full(n_bins, np.nan)
    p_low = np.full(n_bins, np.nan)
    p_high = np.full(n_bins, np.nan)

    bins = full_agg["bin_idx"].to_numpy()
    valid = (bins >= 0) & (bins < n_bins)
    valid_bins = bins[valid]

    p_main[valid_bins] = full_agg["median_logp"].to_numpy()[valid]
    
    stds, cts = full_agg["std_logp"].to_numpy()[valid], full_agg["count_logp"].to_numpy()[valid]
    sem_raw[valid_bins] = stds / np.sqrt(cts)
    
    e_raw[valid_bins] = full_agg["mean_effect"].to_numpy()[valid]

    p_med[valid_bins] = full_agg["median_logp"].to_numpy()[valid]
    p_low[valid_bins] = full_agg["lower_q_logp"].to_numpy()[valid]
    p_high[valid_bins] = full_agg["upper_q_logp"].to_numpy()[valid]

    if np.any(np.isnan(p_main)): p_main[np.isnan(p_main)] = np.nanmean(p_main)
    
    p_smooth = gaussian_filter1d(p_main, sigma=smooth_sigma)
    sem_smooth = gaussian_filter1d(sem_raw, sigma=smooth_sigma)
    
    dist_stats = (
        mask_gaps(p_med, layout_items, bin_size),
        mask_gaps(p_low, layout_items, bin_size),
        mask_gaps(p_high, layout_items, bin_size)
    )

    return x_arr, mask_gaps(p_smooth, layout_items, bin_size), mask_gaps(sem_smooth, layout_items, bin_size), mask_gaps(e_raw, layout_items, bin_size), df, clean_bg_val, (clean_min, clean_max), dist_stats

def process_overlaps_and_candidates(df_glm, heatmap_files, heatmap_col, chrom_order, offsets, total_width, bin_size, fst_pct, glm_pct, overlap_size, layout_items):
    n_bins = int(total_width / bin_size) + 1
    file_map = [(parse_window_size(f), f) for f in heatmap_files]
    target_fst_file = min(file_map, key=lambda x: abs(x[0] - overlap_size))[1]
    df_fst = pl.read_csv(target_fst_file).filter(pl.col("chrom").is_in(chrom_order))
    valid_fst = df_fst[heatmap_col].drop_nans()
    if len(valid_fst) == 0: return np.zeros(n_bins), [], []
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
        if np.any(potential): hits_global.append(high_glm_pos[potential] + offsets[chrom])
    if not hits_global: return np.arange(n_bins) * bin_size, np.zeros(n_bins), []
    all_hits = np.concatenate(hits_global)
    counts, _ = np.histogram(all_hits, bins=n_bins, range=(0, n_bins*bin_size))
    counts = np.nan_to_num(mask_gaps(counts.astype(float), layout_items, bin_size))
    valid_counts = counts[counts > 0]
    if len(valid_counts) == 0: return np.arange(n_bins) * bin_size, counts, []
    count_thresh = np.percentile(valid_counts, 90)
    candidate_indices = np.where(counts >= count_thresh)[0]
    merged_candidates = []
    if len(candidate_indices) > 0:
        curr_cluster = [candidate_indices[0]]
        for i in candidate_indices[1:]:
            if i <= curr_cluster[-1] + 2: curr_cluster.append(i)
            else:
                merged_candidates.append((curr_cluster[0] * bin_size, (curr_cluster[-1] + 1) * bin_size))
                curr_cluster = [i]
        merged_candidates.append((curr_cluster[0] * bin_size, (curr_cluster[-1] + 1) * bin_size))
    return np.arange(n_bins) * bin_size, counts, merged_candidates

def export_candidates_tsv(df_glm, chrom_order, output_base, filtered_candidates, offsets):
    if not filtered_candidates: return
    candidates, site_details = [], []
    for i, (g_start, g_end) in enumerate(filtered_candidates):
        target_chrom, local_start, local_end = None, 0, 0
        for chrom in chrom_order:
            if g_start >= offsets[chrom]:
                target_chrom = chrom
                local_start = g_start - offsets[chrom]
                local_end = g_end - offsets[chrom]
            else: break
        if not target_chrom: continue
        c_glm = df_glm.filter((pl.col("chrom") == target_chrom) & (pl.col("pos") >= local_start) & (pl.col("pos") <= local_end))
        if c_glm.height == 0: continue
        snps_p, snps_e, snps_pos = c_glm["log_p"].to_numpy(), c_glm["effect_size"].to_numpy(), c_glm["pos"].to_numpy()
        reg_id = f"Region_{i+1}"
        candidates.append({
            "Region_ID": reg_id, "Chrom": target_chrom, "Start": int(local_start), "End": int(local_end),
            "Max_LogP": np.max(snps_p), "Mean_Pct_AF_Change": np.mean(snps_e), "SNP_Count": len(snps_pos)
        })
        site_details.append(pl.DataFrame({
            "Region_ID": [reg_id]*len(snps_pos), "Chrom": [target_chrom]*len(snps_pos),
            "Pos": snps_pos, "LogP": snps_p, "Pct_AF_Change": snps_e
        }))
    if candidates: pl.DataFrame(candidates).write_csv(f"{output_base}_candidates.tsv", separator="\t")
    if site_details: pl.concat(site_details).write_csv(f"{output_base}_candidate_sites_detailed.tsv", separator="\t")


def export_interactive_payload(output_path, chroms, offsets, layout_items, total_width, hm_matrix, hm_sizes,
                               glm_x, glm_p, glm_sem, glm_e, over_x, over_y, highlights,
                               rep_x, rep_y, annot_data, bin_size, dist_stats=None):
    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            return super(NumpyEncoder, self).default(obj)

    med, low, high = (None, None, None)
    if dist_stats:
        med, low, high = dist_stats

    payload = {
        "bin_size": bin_size,
        "chrom_order": chroms,
        "offsets": offsets,
        "layout": layout_items,
        "genome_width": total_width,
        "heatmap": {
            "window_sizes": hm_sizes,
            "matrix": hm_matrix.tolist() if hm_matrix is not None else None,
        },
        "glm": {
            "x": glm_x.tolist() if glm_x is not None else None,
            "mean": glm_p.tolist() if glm_p is not None else None,
            "sem": glm_sem.tolist() if glm_sem is not None else None,
            "effect": glm_e.tolist() if glm_e is not None else None,
            "median": med.tolist() if med is not None else None,
            "lower_q": low.tolist() if low is not None else None,
            "upper_q": high.tolist() if high is not None else None
        },
        "overlap": {
            "x": over_x.tolist() if over_x is not None else None,
            "count": over_y.tolist() if over_y is not None else None,
            "highlights": highlights,
        },
        "repeats": {
            "x": rep_x.tolist() if rep_x is not None else None,
            "density": rep_y.tolist() if rep_y is not None else None,
        },
        "annotations": [{"start": a[0], "end": a[1], "label": a[2]} for a in annot_data] if annot_data else [],
    }

    with open(output_path, "w") as f:
        json.dump(payload, f, cls=NumpyEncoder)

def squeeze_dataset(hm_matrix, glm_p, glm_sem, glm_e, over_y, rep_y, annot_data, highlights, layout_items, bin_size, threshold, smooth_sigma, dist_stats=None):
    n_bins = len(glm_p)
    if rep_y is None:
        return hm_matrix, glm_p, glm_sem, glm_e, over_y, rep_y, annot_data, highlights, layout_items, n_bins * bin_size, dist_stats

    smoothed_rep_y = np.zeros_like(rep_y)
    
    for item in layout_items:
        s_bin = int(item['start'] / bin_size)
        e_bin = int((item['start'] + item['width']) / bin_size)
        s_bin = max(0, min(s_bin, n_bins))
        e_bin = max(0, min(e_bin, n_bins))
        
        if e_bin > s_bin:
            chrom_reps = rep_y[s_bin:e_bin]
            chrom_reps = np.nan_to_num(chrom_reps, nan=0.0)
            smoothed_slice = gaussian_filter1d(chrom_reps, sigma=smooth_sigma)
            smoothed_rep_y[s_bin:e_bin] = smoothed_slice

    keep_mask = (smoothed_rep_y <= threshold)

    new_rep_y = rep_y[keep_mask]
    new_glm_p = glm_p[keep_mask]
    new_glm_sem = glm_sem[keep_mask]
    new_glm_e = glm_e[keep_mask]
    new_over_y = over_y[keep_mask]
    new_hm = hm_matrix[:, keep_mask]

    new_dist_stats = None
    if dist_stats:
        med, low, high = dist_stats
        new_med = med[keep_mask] if med is not None else None
        new_low = low[keep_mask] if low is not None else None
        new_high = high[keep_mask] if high is not None else None
        new_dist_stats = (new_med, new_low, new_high)

    new_layout_items = []
    current_x = 0
    
    bin_to_chrom_idx = np.full(n_bins, -1)
    for idx, item in enumerate(layout_items):
        s_bin = int(item['start'] / bin_size)
        e_bin = int((item['start'] + item['width']) / bin_size)
        s_bin = max(0, min(s_bin, n_bins))
        e_bin = max(0, min(e_bin, n_bins))
        bin_to_chrom_idx[s_bin:e_bin] = idx

    kept_chrom_ids = bin_to_chrom_idx[keep_mask]
    
    for i in range(len(layout_items)):
        count = np.sum(kept_chrom_ids == i)
        width = count * bin_size
        new_layout_items.append({
            "start": current_x,
            "width": width,
            "mid": current_x + width/2,
            "label": layout_items[i]['label'],
            "chrom": layout_items[i]['chrom']
        })
        current_x += width

    total_new_width = len(new_glm_p) * bin_size

    cumulative_map = np.cumsum(keep_mask) - 1
    
    def remap_intervals(intervals, is_bed=False):
        new_intervals = []
        for item in intervals:
            if is_bed:
                g_start, g_end, label = item
            else:
                g_start, g_end = item
                label = None
            
            b_start = int(g_start / bin_size)
            b_end = int(g_end / bin_size)
            
            b_start = max(0, min(b_start, n_bins-1))
            b_end = max(0, min(b_end, n_bins-1))
            
            new_b_start = cumulative_map[b_start]
            new_b_end = cumulative_map[b_end]
            
            if keep_mask[b_start] or keep_mask[b_end] or (new_b_end > new_b_start):
                 if new_b_start < 0: new_b_start = 0
                 if new_b_end < new_b_start: new_b_end = new_b_start
                 
                 new_g_start = new_b_start * bin_size
                 new_g_end = (new_b_end + 1) * bin_size # Approximation
                 
                 if is_bed:
                     new_intervals.append((new_g_start, new_g_end, label))
                 else:
                     new_intervals.append((new_g_start, new_g_end))
        return new_intervals

    new_annot_data = remap_intervals(annot_data, is_bed=True)
    new_highlights = remap_intervals(highlights, is_bed=False)

    return new_hm, new_glm_p, new_glm_sem, new_glm_e, new_over_y, new_rep_y, new_annot_data, new_highlights, new_layout_items, total_new_width, new_dist_stats


def export_legends(args, output_base, c_line, c_ribbon, c_repeat, c_overlap, c_gene, cmap):
    fig_leg_hm = plt.figure(figsize=(2, 0.5))
    ax_leg_hm = fig_leg_hm.add_axes([0.1, 0.4, 0.8, 0.4])
    norm = matplotlib.colors.Normalize(vmin=0.9 if args.percentile else 0, vmax=1.0 if args.percentile else 3)
    cb = matplotlib.colorbar.ColorbarBase(ax_leg_hm, cmap=cmap, norm=norm, orientation='horizontal')
    cb.set_label("Mean Z-score (Percentile)" if args.percentile else "Mean Z-score", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    plt.savefig(f"{output_base}_legend_heatmap.png", bbox_inches='tight', dpi=300)
    plt.close(fig_leg_hm)

    fig_leg_tr = plt.figure(figsize=(3, 2))
    ax_leg_tr = fig_leg_tr.add_subplot(111)
    ax_leg_tr.axis('off')
    
    legend_elements = []
    
    if not args.no_trend_line:
        legend_elements.append(Line2D([0], [0], color=c_line, lw=1, label='Smoothed Median'))
    
    if not args.no_sem:
        legend_elements.append(Patch(facecolor=c_ribbon, alpha=args.ribbon_alpha, label=f'Smoothed ± {args.sem_scale:g} SEM'))
        
    if args.plot_bin_stats:
        if args.connect_points:
             legend_elements.append(Line2D([0], [0], marker='o', color='#444444', markerfacecolor='#222222', markersize=4, lw=0.5, label='Bin Median (Connected)'))
        else:
             legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='#222222', markersize=4, label='Bin Median'))
        legend_elements.append(Line2D([0], [0], color='#444444', lw=1.0, label='Bin Spread'))

    legend_elements.extend([
        Patch(facecolor=c_repeat, label='Repeat Density'),
        Patch(facecolor="#0072B2", label='Effect Size (AME)'),
        Patch(facecolor=c_overlap, label='Overlap Count'),
        Patch(facecolor=c_gene, label='Genes')
    ])
    
    ax_leg_tr.legend(handles=legend_elements, loc='center', frameon=False, fontsize=7)
    plt.savefig(f"{output_base}_legend_tracks.png", bbox_inches='tight', dpi=300)
    plt.close(fig_leg_tr)

def render_final(hm_matrix, win_sizes, 
                 glm_x, glm_p, glm_sem, glm_e, gl_bg_mean, gl_bg_range, dist_stats,
                 over_x, over_y, highlight_intervals,
                 rep_x, rep_y,
                 annot_data, annot_labels,
                 layout_items, total_width, bin_size, args, filename, clean_mode=False):
    
    FIG_W, FIG_H = args.fig_width, args.fig_height
    FS_STD = 7 * args.text_scale
    FS_SMALL = 6 * args.text_scale
    
    C_OI_BLUE = "#0072B2"
    C_OI_TEAL = "#009E73"
    C_OI_YELLOW = "#F0E442"
    C_GLM_LINE = "#212121"
    C_GLM_RIBBON = "#BDBDBD" 
    C_REPEAT   = "#555555"
    C_OVERLAP  = C_OI_TEAL
    C_GLM_LAWN = C_OI_BLUE 
    C_GENES    = "#424242"
    cmap_heat = copy.copy(plt.get_cmap("viridis"))
    cmap_heat.set_bad(color='white')

    TRACK_DEFAULTS = {
        'heatmap': 1.2,
        'glm': 0.8,
        'effect': 0.25,
        'overlap': 0.4,
        'repeats': 0.25,
        'genes': 0.15,
        'footer': 0.1
    }

    active_tracks = args.tracks
    if len(args.fig_ratios) == len(active_tracks): ratios = args.fig_ratios
    else: ratios = [TRACK_DEFAULTS.get(t, 0.5) for t in active_tracks]

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)
    gs = GridSpec(len(active_tracks), 1, height_ratios=ratios, 
                  hspace=0.08, left=0.1, right=0.98, top=0.95, bottom=0.05)
    
    axes_dict = {}
    share_ax = None
    
    for i, track_name in enumerate(active_tracks):
        if track_name == 'footer':
            ax = fig.add_subplot(gs[i])
        else:
            if share_ax is None:
                ax = fig.add_subplot(gs[i])
                share_ax = ax
            else:
                ax = fig.add_subplot(gs[i], sharex=share_ax)
        axes_dict[track_name] = ax

    def is_active(name): return name in axes_dict

    if not clean_mode and highlight_intervals:
        for name, ax in axes_dict.items():
            if name == 'footer': continue
            for (start_x, end_x) in highlight_intervals:
                ax.axvspan(start_x, end_x, color=C_OI_YELLOW, alpha=0.3, linewidth=0, zorder=0)
        
        first_track = next((t for t in active_tracks if t != 'footer'), None)
        if first_track:
            ax_top = axes_dict[first_track]
            for i, (start_x, end_x) in enumerate(highlight_intervals):
                mid_x = (start_x + end_x) / 2
                if end_x - start_x > total_width * 0.005:
                    ax_top.text(mid_x, 1.05, f"R{i+1}", transform=ax_top.get_xaxis_transform(),
                                 ha='center', va='bottom', fontsize=FS_SMALL, color='black', rotation=90)

    alphas = np.ones(len(glm_x), dtype=float)
    if rep_y is not None:
        n_points = min(len(alphas), len(rep_y))
        current_reps = np.clip(np.nan_to_num(rep_y[:n_points], nan=0.0), 0.0, 1.0)
        slope = 1.0 - args.min_alpha
        alphas[:n_points] = 1.0 - (current_reps * slope)

    if is_active('heatmap'):
        ax_heat = axes_dict['heatmap']
        vmin, vmax = (0.90, 1.00) if args.percentile else (0, 3.0)
        ax_heat.imshow(hm_matrix, aspect='auto', origin='lower', cmap=cmap_heat,
                       extent=[0, total_width, 0, len(win_sizes)], 
                       vmin=vmin, vmax=vmax, interpolation='nearest', zorder=1)
        
        for item in layout_items:
            ax_heat.add_patch(Rectangle((item["start"], 0), item["width"], len(win_sizes),
                             linewidth=0.5, edgecolor='white', facecolor='none', zorder=2))

        n_ticks = 5
        tick_indices = np.linspace(0, len(win_sizes)-1, n_ticks, dtype=int)
        ax_heat.set_yticks([i + 0.5 for i in tick_indices])
        ax_heat.set_yticklabels([f"{int(win_sizes[i]/1000)}" for i in tick_indices])
        ax_heat.set_ylabel("Window\n(kbp)", rotation=0, ha='right', va='center', fontsize=FS_STD)
        ax_heat.tick_params(bottom=False, labelbottom=False)

    if is_active('glm'):
        ax_glm_p = axes_dict['glm']
        
        if args.plot_bin_stats and dist_stats:
            d_med, d_low, d_high = dist_stats
            mask = ~np.isnan(d_med)
            
            if np.any(mask):
                ax_glm_p.vlines(glm_x[mask], d_low[mask], d_high[mask], 
                                colors='#444444', alpha=args.bin_stats_alpha, linewidths=0.8, zorder=4)
                
                if args.connect_points:
                    ax_glm_p.plot(glm_x, d_med, color='#444444', linewidth=0.5, alpha=0.8, zorder=4.5)

                ax_glm_p.scatter(glm_x[mask], d_med[mask], 
                                 s=4, c='#222222', alpha=1.0, zorder=5, linewidths=0, rasterized=True)

        if glm_sem is not None and not args.no_sem:
            scaled_sem = glm_sem * args.sem_scale
            rib_h = 2 * scaled_sem
            rib_b = glm_p - scaled_sem
            rib_h = np.nan_to_num(rib_h, nan=0.0)
            rib_b = np.nan_to_num(rib_b, nan=0.0)
            rib_c = np.zeros((len(glm_x), 4))
            rib_c[:] = mcolors.to_rgba(C_GLM_RIBBON)
            rib_c[:, 3] = args.ribbon_alpha 
            ax_glm_p.bar(glm_x, rib_h, bottom=rib_b, width=bin_size, align='edge',
                         color=rib_c, linewidth=0, rasterized=True, zorder=2)
            
            if not clean_mode:
                ax_glm_p.text(0.01, 0.9, f"SEM x{args.sem_scale:g}", transform=ax_glm_p.transAxes,
                              ha='left', va='top', fontsize=FS_SMALL, fontweight='bold', color='black')

        if not args.no_trend_line:
            points = np.array([glm_x, glm_p]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc_rgba = mcolors.to_rgba(C_GLM_LINE)
            lc_colors = np.zeros((len(segments), 4))
            lc_colors[:] = lc_rgba
            lc_colors[:, 3] = alphas[:len(segments)] 
            lc = LineCollection(segments, colors=lc_colors, linewidth=0.8, zorder=3)
            ax_glm_p.add_collection(lc)

        ax_glm_p.axhline(y=gl_bg_mean, color='#999999', linestyle=':', linewidth=0.5, zorder=1)
        ax_glm_p.set_ylabel(r"Average" + "\n" + r"$-\log_{10}(P)$", rotation=0, ha='right', va='center', fontsize=FS_STD)
        ax_glm_p.tick_params(bottom=False, labelbottom=False)
        
        if args.y_max: ax_glm_p.set_ylim(args.y_min, args.y_max)
        else:
            clean_min, clean_max = gl_bg_range
            if clean_max == 0: clean_max = 1.0
            ax_glm_p.set_ylim(clean_min * 0.9 if clean_min > 0 else 0, clean_max * 1.2)

    if is_active('effect'):
        ax_glm_e = axes_dict['effect']
        valid_e = ~np.isnan(glm_e)
        if np.any(valid_e):
            eff_colors = np.zeros((len(glm_x), 4))
            eff_colors[:] = mcolors.to_rgba(C_GLM_LAWN)
            eff_colors[:, 3] = alphas 
            safe_e = np.nan_to_num(glm_e, nan=0.0)
            ax_glm_e.bar(glm_x, safe_e, width=bin_size, align='edge', 
                         color=eff_colors, linewidth=0, rasterized=True, zorder=1)
        
        if not clean_mode:
            ax_glm_e.text(0.01, 0.9, "Average Marginal Effect", transform=ax_glm_e.transAxes,
                         ha='left', va='top', fontsize=FS_SMALL, fontweight='normal', color='black')
        
        ax_glm_e.set_ylabel(r"% $\Delta$AF/gen", rotation=0, ha='right', va='center', fontsize=FS_SMALL)
        ax_glm_e.tick_params(bottom=False, labelbottom=False)
        ax_glm_e.yaxis.set_major_locator(mticker.MaxNLocator(nbins=2))

    if is_active('overlap'):
        ax_over = axes_dict['overlap']
        safe_over = np.nan_to_num(over_y, nan=0.0)
        if np.max(safe_over) > 0:
            over_colors = np.zeros((len(glm_x), 4))
            over_colors[:] = mcolors.to_rgba(C_OVERLAP)
            over_colors[:, 3] = alphas 
            ax_over.bar(over_x, safe_over, width=bin_size, align='edge', 
                        color=over_colors, linewidth=0, zorder=1)
        
        ax_over.set_ylabel("Overlap\nCount", rotation=0, ha='right', va='center', fontsize=FS_SMALL)
        
        if not clean_mode:
            fst_str = f"{args.fst_top_pct:g}"
            glm_str = f"{args.glm_top_pct:g}"
            annotation_text = rf"Top {fst_str}% $\mathrm{{F}}_{{ST}}$ & {glm_str}% GLM"
            ax_over.text(0.99, 0.9, annotation_text, transform=ax_over.transAxes, 
                         ha='right', va='top', fontsize=FS_SMALL, fontweight='normal')
        
        ax_over.tick_params(bottom=False, labelbottom=False)
        ax_over.yaxis.set_major_locator(mticker.MaxNLocator(nbins=3, integer=True))

    if is_active('repeats'):
        ax_rep = axes_dict['repeats']
        if rep_x is not None:
            rep_y_masked = mask_gaps(rep_y, layout_items, bin_size)
            rep_colors = np.zeros((len(rep_x), 4))
            rep_colors[:] = mcolors.to_rgba(C_REPEAT)
            rep_colors[:, 3] = alphas[:len(rep_x)] 
            ax_rep.bar(rep_x, rep_y_masked, width=bin_size, align='edge',
                       color=rep_colors, linewidth=0, rasterized=True, zorder=1)
            
        ax_rep.set_ylim(0, 1)
        ax_rep.set_yticks([0, 0.5, 1])
        ax_rep.set_yticklabels(["0", "", "1"], fontsize=FS_SMALL)
        ax_rep.set_ylabel("Repeat\nDensity", rotation=0, ha='right', va='center', fontsize=FS_STD)
        ax_rep.tick_params(bottom=False, labelbottom=False)

    if is_active('genes'):
        ax_annot = axes_dict['genes']
        ax_annot.set_ylim(0, 1)
        ax_annot.axis('off')
        if annot_data:
            bar_y, bar_h = 0.3, 0.4
            min_w = total_width * 0.002
            for (g_start, g_end, label) in annot_data:
                w = max(g_end - g_start, min_w)
                ax_annot.add_patch(Rectangle((g_start, bar_y), w, bar_h, facecolor=C_GENES, edgecolor='none', zorder=1))
        ax_annot.text(-0.01, 0.5, "Genes", transform=ax_annot.transAxes, ha='right', va='center', fontsize=FS_STD, fontweight='normal')

    if is_active('footer'):
        ax_foot = axes_dict['footer']
        ax_foot.set_xlim(0, total_width)
        ax_foot.set_ylim(0, 1)
        ax_foot.axis('off')
        for i, item in enumerate(layout_items):
            if (item["width"]/total_width) > 0.005: 
                ax_foot.plot([item["start"], item["start"] + item["width"]], [0.8, 0.8], color="black", linewidth=1.5, solid_capstyle='butt')
                y_pos = 0.45 if i % 2 == 0 else 0.0
                ax_foot.text(item["mid"], y_pos, item["label"], ha='center', va='top', fontsize=FS_STD)

    for name, ax in axes_dict.items():
        ax.set_xlim(0, total_width)

    plt.savefig(filename, bbox_inches='tight', dpi=300)
    plt.close(fig)

    export_legends(args, filename.rsplit('.', 1)[0], C_GLM_LINE, C_GLM_RIBBON, C_REPEAT, C_OVERLAP, C_GENES, cmap_heat)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--heatmap-files", nargs='+', required=True)
    parser.add_argument("--glm", required=True)
    parser.add_argument("--target-pval", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--svg-output", type=str, help="Optional path to save SVG version of the plot.")
    parser.add_argument("--json-output", type=str, help="Optional path to save plot-ready data for interactive editing.")
    parser.add_argument("--genome-fasta", type=str)
    parser.add_argument("--repeat-masker", type=str)
    parser.add_argument("--chrom-map", type=str)
    parser.add_argument("--bed-file", type=str)
    parser.add_argument("--heatmap-col", default="z_score_median")
    parser.add_argument("--bin-size", type=int, default=100000)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--gap", type=int, default=1_000_000)
    parser.add_argument("--min-maf", type=float, default=0.01)
    parser.add_argument("--min-window", type=int, default=20000)
    parser.add_argument("--max-window", type=int, default=500000)
    parser.add_argument("--percentile", action="store_true")
    parser.add_argument("--fst-top-pct", type=float, default=5.0)
    parser.add_argument("--glm-top-pct", type=float, default=5.0)
    parser.add_argument("--overlap-target-size", type=int, default=250000)
    parser.add_argument("--y-max", type=float, default=None)
    parser.add_argument("--y-min", type=float, default=0)
    parser.add_argument("--ribbon-alpha", type=float, default=0.4)
    parser.add_argument("--smooth-sigma", type=float, default=1.0)
    parser.add_argument("--sem-scale", type=float, default=3.0)
    parser.add_argument("--min-alpha", type=float, default=0.4)
    parser.add_argument("--text-scale", type=float, default=1.0)
    parser.add_argument("--fig-width", type=float, default=7.5)
    parser.add_argument("--fig-height", type=float, default=7.0)
    parser.add_argument("--fig-ratios", nargs='+', type=float, default=[], 
                        help="Custom height ratios. Must match number of tracks if provided.")
    parser.add_argument("--tracks", nargs='+', 
                        default=["heatmap", "glm", "effect", "overlap", "repeats", "genes", "footer"],
                        choices=["heatmap", "glm", "effect", "overlap", "repeats", "genes", "footer"],
                        help="Select tracks to plot and their order.")

    parser.add_argument("--remove-repeats", action="store_true", 
                        help="Remove bins with high repeat density from the plot.")
    parser.add_argument("--repeat-threshold", type=float, default=0.4,
                        help="Threshold (0.0-1.0) for removing repeats.")
    parser.add_argument("--repeat-filter-sigma", type=float, default=100.0,
                        help="Smoothing sigma (in bins) for detecting macro repeat blocks.")
    

    parser.add_argument("--plot-bin-stats", action="store_true",
                        help="Plot median GLM point and spread intervals (quantile lines) for each bin.")
    parser.add_argument("--bin-stats-alpha", type=float, default=1.0,
                        help="Opacity (alpha) of the vertical spread lines. Default: 1.0")
    parser.add_argument("--stats-lower-q", type=float, default=0.25,
                        help="Lower quantile for point spread (default: 0.25 = 25th percentile).")
    parser.add_argument("--stats-upper-q", type=float, default=0.75,
                        help="Upper quantile for point spread (default: 0.75 = 75th percentile).")
    
    parser.add_argument("--no-sem", action="store_true",
                        help="Do not plot the SEM ribbon.")
    parser.add_argument("--no-trend-line", action="store_true",
                        help="Do not plot the smoothed trend line.")
    parser.add_argument("--connect-points", action="store_true",
                        help="Draw a thin line connecting the bin median points.")

    return parser.parse_args()

def main():
    args = parse_args()

    if args.fig_ratios and len(args.fig_ratios) != len(args.tracks):
        sys.exit(f"Error: Provided --fig-ratios ({len(args.fig_ratios)}) does not match number of --tracks ({len(args.tracks)}).")

    output_base = os.path.splitext(args.output)[0]
    svg_out = args.svg_output or f"{output_base}.svg"
    json_out = args.json_output or f"{output_base}.json"

    chroms, offsets, layout, width, fasta_sizes = build_genome_structure(args.heatmap_files, args.top_n, args.gap, args.genome_fasta)
    hm_matrix, sizes, sorted_files = process_heatmap_binned(args.heatmap_files, args.heatmap_col, chroms, offsets, width, args.bin_size, args.percentile, args.min_window, args.max_window, layout)
    chrom_map = load_chrom_map(args.chrom_map)
    rep_x, rep_y = None, None
    if args.repeat_masker and fasta_sizes:
        rep_x, rep_y = process_repeats_binned(args.repeat_masker, chrom_map, chroms, offsets, width, args.bin_size, fasta_sizes)
    
    glm_x, glm_p, glm_sem, glm_e, raw_glm_df, mean_bg, bg_range, dist_stats = process_glm_smoothed_mean(
        args.glm, args.target_pval, chroms, offsets, width, args.bin_size, args.min_maf, rep_y, layout, 
        smooth_sigma=args.smooth_sigma, lower_q=args.stats_lower_q, upper_q=args.stats_upper_q
    )

    over_x, over_y, highlights = process_overlaps_and_candidates(
        raw_glm_df, sorted_files, args.heatmap_col, chroms, offsets, width, args.bin_size, args.fst_top_pct, args.glm_top_pct, args.overlap_target_size, layout
    )
    annot_data, annot_labels = load_bed_annotations(args.bed_file, chroms, offsets)
    
    if args.remove_repeats:
        if rep_y is None:
            sys.exit("Error: --remove-repeats requires --repeat-masker input.")

        hm_matrix, glm_p, glm_sem, glm_e, over_y, rep_y, annot_data, highlights, layout, width, dist_stats = squeeze_dataset(
            hm_matrix, glm_p, glm_sem, glm_e, over_y, rep_y, annot_data, highlights, layout, args.bin_size, args.repeat_threshold, args.repeat_filter_sigma,
            dist_stats=dist_stats
        )
        n_bins = len(glm_p)
        glm_x = np.arange(n_bins) * args.bin_size
        over_x = glm_x 
        rep_x = glm_x  

    render_final(hm_matrix, sizes, glm_x, glm_p, glm_sem, glm_e, mean_bg, bg_range, dist_stats,
                 over_x, over_y, highlights,
                 rep_x, rep_y, annot_data, annot_labels, layout, width, args.bin_size, args, args.output, clean_mode=False)
    
    clean_filename = f"{output_base}_clean.png"
    render_final(hm_matrix, sizes, glm_x, glm_p, glm_sem, glm_e, mean_bg, bg_range, None, # Clean mode doesn't need points
                 over_x, over_y, None, 
                 rep_x, rep_y, annot_data, annot_labels, layout, width, args.bin_size, args, clean_filename, clean_mode=True)

    if svg_out:
        render_final(hm_matrix, sizes, glm_x, glm_p, glm_sem, glm_e, mean_bg, bg_range, dist_stats,
                     over_x, over_y, highlights,
                     rep_x, rep_y, annot_data, annot_labels, layout, width, args.bin_size, args, svg_out, clean_mode=False)

    export_interactive_payload(json_out, chroms, offsets, layout, width, hm_matrix, sizes,
                               glm_x, glm_p, glm_sem, glm_e, over_x, over_y, highlights,
                               rep_x, rep_y, annot_data, args.bin_size, dist_stats)

    export_candidates_tsv(raw_glm_df, chroms, output_base, highlights, offsets)

if __name__ == "__main__":
    main()