import os
import glob

#PILEUP = glob_wildcards(os.path.join(config["input_dir"], "{pileup}.mpileup.gz")).pileup
BAM = glob_wildcards(os.path.join(config["input_dir"], "{bam}.bam")).bam

def fb(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            pool_size = parts[1]
            mapping[key] = pool_size
    return mapping

FB_DICT = fb(config["fb_map"])
FB = list(FB_DICT.keys())

# with open(config["fb_map"], 'r') as file:
#     FB = [line.strip() for line in file if not line.strip().startswith('#') and line.strip()]


############################################################################################################################################################################################################
##### merge wildcard ####
############################################################################################################################################################################################################

def merge_dict(merge_map):
    mapping = {}
    with open(merge_map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            key, value = line.split(': ')
            mapping[key] = value.split('/')
    return mapping

MERGE_DICT = merge_dict(config["merge_map"])
MERGE = list(MERGE_DICT.keys())

def pileup(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            samples = parts[1].split(',')
            mapping[key] = samples
    return mapping

PILEUP_DICT = pileup(config["pileup_map"])
PILEUP = list(PILEUP_DICT.keys())


############################################################################################################################################################################################################
##### FST Wildcards ####
############################################################################################################################################################################################################

# def process_fst(line):
#     items = line.strip().split(',')
#     return f"{items[0]}-{items[1]}"
#     
# def get_comparend1(wildcards):
#     items = wildcards.fst.split('-')
#     return f"{items[0]}:1"
#     
# def get_comparend2(wildcards):
#     items = wildcards.fst.split('-')
#     return f"{items[1]}:1"
# 
# 
# def get_comparend1_name(wildcards):
#     items = wildcards.fst.split('-')
#     return f"{items[0]}"
# 
# def get_comparend2_name(wildcards):
#     items = wildcards.fst.split('-')
#     return f"{items[1]}"
# 
# with open(config["fst_map"], 'r') as file:
#     FST = [process_fst(line) for line in file if not line.strip().startswith('#') and line.strip()]

# def poolfstat(map):
#     mapping = {}
#     with open(map) as f:
#         for line in f:
#             line = line.strip()
#             if line.startswith('#') or not line:  # Skip commented and empty lines
#                 continue
#             parts = line.split(': ')
#             key = parts[0]
#             value = parts[1].split(',')
#             variants = parts[2]
#             mapping[key] = (value, variants)
#     return mapping


def poolfstat(mapfile):
    mapping = {}
    with open(mapfile) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = [p.strip() for p in line.split(':')]
            key = parts[0]
            samples = parts[1].split(',')
            variant = parts[2]
            groups = parts[3]
            mapping[key] = {
                "samples": samples,
                "variant": variant,
                "groups": groups
            }
    return mapping


POOLFSTAT_DICT = poolfstat(config["poolfstat_map"])
POOLFSTAT = list(POOLFSTAT_DICT.keys())


def grenfst(mapfile):
    mapping = {}
    with open(mapfile) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = [p.strip() for p in line.split(':')]
            key = parts[0]
            variants = parts[1]
            window_queue_count = parts[2]
            window_queue_stride = parts[3]
            filter_total_snp_min_frequency = parts[4]
            mapping[key] = (variants, window_queue_count, window_queue_stride, filter_total_snp_min_frequency)
    return mapping

GRENFST_DICT = grenfst(config["grenfst_map"])
GRENFST = list(GRENFST_DICT.keys())


# def grenfst_multiplot(mapfile):
#     mapping = {}
#     with open(mapfile) as f:
#         for line in f:
#             line = line.strip()
#             if line.startswith('#') or not line:
#                 continue
#             parts = [p.strip() for p in line.split(':')]
#             key = parts[0]
#             fst = parts[1]
#             samps = parts[2].strip().split('/')
#             mapping[key] = (fst, samps)
#     return mapping
# 
# GRENFST_MULTIPLOT_DICT = grenfst_multiplot(config["grenfst_multiplot_map"])
# GRENFST_MULTIPLOT = list(GRENFST_MULTIPLOT_DICT.keys())

def grenfst_multiplot(mapfile):
    mapping = {}
    with open(mapfile) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = [p.strip() for p in line.split(':')]
            
            # Key = Analysis Name (e.g., E10fe9wGREN_btwTB)
            key = parts[0]
            
            # Val 0 = VCF Key (e.g., E10fe9wGREN)
            vcf_key = parts[1]
            
            # Val 1 = Target Replicates (e.g., ['T1','T2','T3','T4'])
            target_reps = [x.strip() for x in parts[2].split(',')]
            
            # Val 2 = Reference Replicates (e.g., ['B1','B2','B3','B4'])
            ref_reps = [x.strip() for x in parts[3].split(',')]
            
            # Val 3 = Generations (e.g., ['01','02','06'...])
            gens = [x.strip() for x in parts[4].split(',')]
            
            mapping[key] = (vcf_key, target_reps, ref_reps, gens)
    return mapping

GRENFST_MULTIPLOT_DICT = grenfst_multiplot(config["grenfst_multiplot_map"])
GRENFST_MULTIPLOT = list(GRENFST_MULTIPLOT_DICT.keys())







############################################################################################################################################################################################################
##### CMH Wildcards ####
############################################################################################################################################################################################################


def cmh_samples(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            key, value, comps = line.split(': ')
            mapping[key] = (value.split(','))
    return mapping

def cmh_comps(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            key, value, comps = line.split(': ')
            mapping[key] = (comps.split('/'))
    return mapping


CMH_SAMPLES_DICT = cmh_samples(config["cmh_map"])
CMH_COMPS_DICT = cmh_comps(config["cmh_map"])
CMH = list(CMH_SAMPLES_DICT.keys())


def cmh2(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            data = parts[1]
            trt = parts[2],
            add_freq = parts[3]
            mapping[key] = (data, trt, add_freq)
    return mapping

CMH2_DICT = cmh2(config["cmh2_map"])
CMH2 = list(CMH2_DICT.keys())


############################################################################################################################################################################################################
##### Ind Wildcards ####
############################################################################################################################################################################################################

def indmap(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            variants = parts[1]
            samples = parts[2].split(',')
            mapping[key] = (variants, samples)
    return mapping

IND_DICT = indmap(config["ind_map"])
IND = list(IND_DICT.keys())


############################################################################################################################################################################################################
##### Variants Wildcards ####
############################################################################################################################################################################################################


def poolsnp(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            samples = parts[1]
            params = parts[2].split(',')
            pfilter = parts[3].split(',')
            mapping[key] = (samples, params, pfilter)
    return mapping

POOLSNP_DICT = poolsnp(config["poolsnp_map"])
POOLSNP = list(POOLSNP_DICT.keys())

# def snape(map):
#     with open(map, 'r') as file:
#         return [line.strip() for line in file]

#SNAPE = snape(config["snape_map"])

def variants(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            key, value = line.split(': ')
            mapping[key] = (value.split(','))
    return mapping

VARIANTS_DICT = variants(config["variants_map"])
VARIANTS = list(VARIANTS_DICT.keys())


############################################################################################################################################################################################################
##### fvariants Wildcards ####
############################################################################################################################################################################################################

# def pfilter(map):
#     mapping = {}
#     with open(map) as f:
#         for line in f:
#             line = line.strip()
#             if line.startswith('#') or not line:  # Skip commented and empty lines
#                 continue
#             parts = line.split(': ')
#             key = parts[0]
#             inputs = parts[1].split(',')
#             mapping[key] = inputs
#     return mapping
# 
# PFILTER_DICT = pfilter(config["pfilter_map"])
# PFILTER = list(PFILTER_DICT.keys())


def fvariants(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            inputs = parts[1].split(',')
            variants = parts[2]
            correction = parts[3]
            mapping[key] = (inputs, variants, correction)
    return mapping

FVAR_DICT = fvariants(config["fvariants_map"])
FVAR = list(FVAR_DICT.keys())


############################################################################################################################################################################################################
##### GLM Wildcards ####
############################################################################################################################################################################################################


def glm(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            data = parts[1]
            trt = parts[2],
            add_freq = parts[3]
            contrast_base = parts[4]
            mapping[key] = (data, trt, add_freq, contrast_base)
    return mapping

GLM_DICT = glm(config["glm_map"])
GLM = list(GLM_DICT.keys())


def glmcomp(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            in1 = parts[1].split(',')
            in2 = parts[2].split(',')
            in1name = parts[3]
            in2name = parts[4]
            in3comp = parts[5].split(',')
            mapping[key] = (in1, in2, in1name, in2name, in3comp)
    return mapping


GLMCOMP_DICT = glmcomp(config["glmcomp_map"])
GLMCOMP = list(GLMCOMP_DICT.keys())


def glm_permuted(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            original_glm = parts[1]
            order = parts[2]
            mapping[key] = (original_glm, order)
    return mapping

GLMPERM_DICT = glm_permuted(config["glm_permute_map"])
GLMPERM = list(GLMPERM_DICT.keys())


def fdr_permuted(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            original_glm = parts[1]
            permutation = parts[2].strip().split(',')
            term = parts[3]
            mapping[key] = (original_glm, permutation, term)
    return mapping

FDR_PERMUTED_DICT = fdr_permuted(config["fdr_permuted_map"])
FDR_PERMUTED = list(FDR_PERMUTED_DICT.keys())


def glm_window(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            window_size = parts[1]
            glm = parts[2]
            fdr = parts[3]
            mapping[key] = (window_size, glm, fdr)
    return mapping

W_DICT = glm_window(config["glm_window"])
W = list(W_DICT.keys())

def gglm(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            size = parts[1]
            glm = parts[2]
            fdr = parts[3]
            mapping[key] = (size, glm, fdr)
    return mapping

G_DICT = gglm(config["gglm_map"])
G = list(G_DICT.keys())


def dg(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            delta = parts[1]
            column = parts[2]
            cutoff = parts[3]
            fdr = parts[4]
            glm = parts[5]
            mapping[key] = (delta, column, cutoff, fdr, glm)
    return mapping

DG_DICT = dg(config["dg_map"])
DG = list(DG_DICT.keys())






############################################################################################################################################################################################################
##### DELTA Wildcards ####
############################################################################################################################################################################################################

def delta(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            delta = parts[1]
            freq = parts[2]
            reference = parts[3]
            #delta_files = parts[2].split(',')            
            #concat_names = []
            #for group in delta_files:
            #    samples = group.split('|')
            #    concat_name = ''.join(samples)
            #    concat_names.append(concat_name)
            mapping[key] = (delta, freq, reference)
    return mapping


DELTA_DICT = delta(config["delta_map"])
DELTA = list(DELTA_DICT.keys())

############################################################################################################################################################################################################
##### DELTA Wildcards ####
############################################################################################################################################################################################################

def ly(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            wild_cuttoff = parts[1]
            delta = parts[2]
            mapping[key] = (wild_cuttoff, delta)
    return mapping


LY_DICT = ly(config["lynch_map"])
LY = list(LY_DICT.keys())




############################################################################################################################################################################################################
##### PCA Wildcards ####
############################################################################################################################################################################################################

def pca(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            freq = parts[1]
            wild_sites = parts[2]
            glm_freq = parts[3]
            mapping[key] = (freq, wild_sites, glm_freq)
    return mapping


PCA_DICT = pca(config["pca_map"])
PCA = list(PCA_DICT.keys())



############################################################################################################################################################################################################
##### Freq Wildcards ####
############################################################################################################################################################################################################

def process_fst(line):
    items = line.strip().split(',')
    return f"{items[0]}-{items[1]}"
    
def freq_dict(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            key, value, comp = line.split(': ')
            mapping[key] = (value.split(','), process_fst(comp))
    return mapping

FREQS_DICT = freq_dict(config["freq_map"])
FREQS = list(FREQS_DICT.keys())


def freq_change_dict(map_file):
    mapping = {}
    with open(map_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            key, value = [item.strip() for item in line.split(':')]
            mapping[key] = [v.strip() for v in value.split(',')]
    return mapping

FREQ_CHANGE_DICT = freq_change_dict(config["freq_change_map"])
FREQ_CHANGE = list(FREQ_CHANGE_DICT.keys())


############################################################################################################################################################################################################
##### Poolsim Wildcards ####
############################################################################################################################################################################################################

def poolsim(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            variants = parts[1]
            mapping[key] = variants
    return mapping

POOLSIM_DICT = poolsim(config["poolsim_map"])
POOLSIM = list(POOLSIM_DICT.keys())


############################################################################################################################################################################################################
##### Levene Wildcards ####
############################################################################################################################################################################################################

def levene(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            key, value = line.split(': ')
            mapping[key] = (value.split(','))
    return mapping

LEVENE_DICT = levene(config["levene_map"])
LEVENE = list(LEVENE_DICT.keys())


############################################################################################################################################################################################################
##### TajimasD Wildcards ####
############################################################################################################################################################################################################

def taj(map):
    mapping = {}
    with open(map) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:  # Skip commented and empty lines
                continue
            parts = line.split(': ')
            key = parts[0]
            fvar = parts[1]
            samples = parts[2].split(',')
            params = parts[3].split(',')
            mapping[key] = (fvar, samples, params)
    return mapping

TAJ_DICT = taj(config["taj_map"])
TAJ = list(TAJ_DICT.keys())