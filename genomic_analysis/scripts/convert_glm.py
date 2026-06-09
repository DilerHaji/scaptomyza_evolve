import polars as pl
import gzip
import csv
import argparse
import sys
from typing import List, Dict, Tuple, Optional

vcf_headers = []

def parse_vcf_line(line: str) -> Optional[Tuple[str, int, Dict[str, Tuple[float, int]]]]:
    try:
        parts = line.strip().split('\t')
        if len(parts) < 10:
            return None
            
        chrom = parts[0]
        try:
            pos = int(parts[1])
        except ValueError:
            return None
        
        format_str = parts[8]
        format_fields = format_str.split(':')
        
        try:
            freq_idx = format_fields.index('FREQ')
            dp_idx = format_fields.index('DP')
        except ValueError:
            return None
        
        sample_data = {}
        
        for i, sample_str in enumerate(parts[9:], start=9):
            if i >= len(vcf_headers):
                break
                
            if sample_str.startswith('.'): 
                continue

            sample_values = sample_str.split(':')
            
            if len(sample_values) != len(format_fields):
                continue
            
            freq_str = sample_values[freq_idx]
            
            if ',' in freq_str:
                return None
            
            dp_str = sample_values[dp_idx]
            
            if freq_str == '.' or dp_str == '.':
                continue

            try:
                freq = float(freq_str.rstrip('%')) / 100.0 if '%' in freq_str else float(freq_str)
                dp = int(dp_str)
                
                sample_name = vcf_headers[i]
                sample_data[sample_name] = (freq, dp)
            except ValueError:
                continue
                
        if not sample_data:
            return None
            
        return chrom, pos, sample_data
        
    except Exception:
        return None

def load_metadata(metadata_file: str, treatments: List[str]) -> Dict[str, Dict]:
    metadata = {}
    try:
        with open(metadata_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if treatments is None or row['trt'] in treatments:
                    try:
                        metadata[row['sample']] = {
                            'gen': int(row['gen']),
                            'pop': int(row['pop']),
                            'trt': row['trt']
                        }
                    except (ValueError, KeyError) as e:
                        pass
    except Exception as e:
        sys.exit(1)

    return metadata

def process_vcf(vcf_file: str, metadata_file: str, output_file: str, treatments: List[str]):
    metadata = load_metadata(metadata_file, treatments)
    if not metadata:
        sys.exit(1)

    cols = {
        'CHROM': [],
        'POS': [],
        '1.FREQ': [],
        'neff': [],
        'gen': [],
        'pop': [],
        'trt': [],
        'SourceFile': []
    }

    skipped_variants = 0
    total_variants = 0
    processed_sites = 0
    
    global vcf_headers
    
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    
    try:
        with open_func(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    if line.startswith('#CHROM'):
                        vcf_headers = line.strip().split('\t')
                    continue

                if not vcf_headers:
                    continue

                total_variants += 1

                parsed = parse_vcf_line(line)
                
                if parsed is None:
                    skipped_variants += 1
                    continue
                
                chrom, pos, sample_dict = parsed
                

                site_has_data = False
                for sample_name, (freq, dp) in sample_dict.items():
                    if sample_name in metadata:
                        meta = metadata[sample_name]

                        cols['CHROM'].append(chrom)
                        cols['POS'].append(pos)
                        cols['1.FREQ'].append(freq)
                        cols['neff'].append(dp)
                        cols['gen'].append(meta['gen'])
                        cols['pop'].append(meta['pop'])
                        cols['trt'].append(meta['trt'])
                        cols['SourceFile'].append(sample_name)
                        
                        site_has_data = True
                
                if site_has_data:
                    processed_sites += 1

    except Exception as e:
        sys.exit(1)

    try:
        df = pl.DataFrame(cols)
        
        if df.height == 0:
            with open(output_file, 'w') as f:
                f.write(','.join(cols.keys()) + '\n')
            return

        df = df.sort(['CHROM', 'POS', 'SourceFile'])

        df.write_csv(output_file)

    except Exception as e:
        sys.exit(1)

def parse_treatments(treatments_str: str) -> List[str]:
    if not treatments_str:
        return None
    return [t.strip() for t in treatments_str.split(',')]

def main():
    parser = argparse.ArgumentParser(description='Convert VCF to GLM CSV')
    parser.add_argument('-v', '--vcf', required=True)
    parser.add_argument('-m', '--metadata', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-t', '--treatments', help='Comma-separated treatments (e.g., "B,T")')

    args = parser.parse_args()
    treatments = parse_treatments(args.treatments)
    
    process_vcf(args.vcf, args.metadata, args.output, treatments)

if __name__ == "__main__":
    main()