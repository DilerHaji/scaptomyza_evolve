import polars as pl
import gzip
import argparse
from pathlib import Path

def read_correction_file(file_path):
    return pl.read_csv(file_path)

def process_vcf(vcf_path, correction_data, output_path, correction_column):
    
    correction_samples = set(correction_data['SourceFile'].unique().to_list())
    
    total_variants = 0
    biallelic_variants = 0
    corrected_variants = 0
    multiallelic_variants = 0
    set_to_missing = 0
    
    with gzip.open(vcf_path, 'rt') as vcf, gzip.open(output_path, 'wt') as out:
        for line in vcf:
            if line.startswith('##'):
                out.write(line)
                continue
            
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
                sample_to_index = {sample: i+9 for i, sample in enumerate(header[9:])}
                out.write(line)  # Write the original header without modification
                continue
            
            total_variants += 1
            
            fields = line.strip().split('\t')
            chrom, pos, id_, ref, alt, qual, filter_, info = fields[:8]
            format_ = fields[8]
            
            if len(alt.split(',')) > 1:
                multiallelic_variants += 1
                out.write(line)  # Write multiallelic variants without modification
                continue
            
            biallelic_variants += 1
            pos = int(pos)
            
            correction = correction_data.filter(
                (pl.col('CHROM') == chrom) & 
                (pl.col('POS') == pos)
            )

            new_fields = fields[:9]  # Keep the first 9 columns
            corrected = False
            
            for sample, index in sample_to_index.items():
                
                if sample in correction_samples:
                    sample_correction = correction.filter(pl.col('SourceFile') == sample)

                    if len(sample_correction) > 0:
                        corrected_dp = sample_correction[str(correction_column)].item()
                        
                        genotype_data = fields[index].split(':')
                        
                        if len(genotype_data) >= 5:
                            try:
                                old_gt, old_rd, old_ad, old_dp, old_freq = genotype_data[:5]
                                new_dp = round(corrected_dp)
                                new_rd = round(new_dp * (1 - float(old_freq)))
                                new_ad = round(new_dp * float(old_freq))
                                new_genotype_data = [
                                    old_gt,
                                    str(new_rd),
                                    str(new_ad),
                                    str(new_dp),
                                    str(old_freq) if new_dp > 0 else "0.00"
                                ]
                                new_fields.append(':'.join(new_genotype_data))
                                corrected = True
                            except ValueError:
                                new_fields.append('./.:.:.:.:.') # If any value is '.', set to missing
                                set_to_missing += 1
                        else:
                            new_fields.append('./.:.:.:.:.') # Set to missing for samples with abnormal genotype structure 
                            set_to_missing += 1
                    else:
                        new_fields.append(fields[index]) # Set to missing for samples not in SourceFile column of correction file
                        set_to_missing += 1
                else:
                    new_fields.append(fields[index])  # Set to missing for samples not in SourceFile column of correction file
                    set_to_missing += 1
            
            if corrected:
                corrected_variants += 1
            out.write('\t'.join(new_fields) + '\n')
    
def main():
    parser = argparse.ArgumentParser(description='Process VCF file with corrections, keeping only biallelic sites and samples matching SourceFile in correction files.')
    parser.add_argument('--vcf', required=True, help='Input VCF file path')
    parser.add_argument('--output', required=True, help='Output VCF file path')
    parser.add_argument('--correction_files', nargs='+', required=True, help='Correction CSV file paths')
    parser.add_argument('--correction_column', required=True, help='Column name in correction files to use for DP correction')
    
    args = parser.parse_args()

    correction_data = pl.concat([read_correction_file(file) for file in args.correction_files])
    
    process_vcf(args.vcf, correction_data, args.output, args.correction_column)

if __name__ == '__main__':
    main()