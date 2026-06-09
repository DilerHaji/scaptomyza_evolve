import csv
import argparse
import os

def read_metadata(metadata_file):
    metadata = {}
    with open(metadata_file, mode='r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)  # Skip the header
        for row in reader:
            if row:  # Check if row is not empty
                filename, gen, trt, pop = row[0], row[1], row[2], row[3]
                metadata[filename] = {'gen': gen, 'trt': trt, 'pop': pop}
    return metadata

def process_files(files, metadata, output_file):
    with open(output_file, 'w') as outfile:
        writer = csv.writer(outfile)
        header = ['CHROM', 'POS', 'REF', 'ALT', '1.REF_CNT', '1.ALT_CNT', '1.COV', '1.FREQ', 'Czech2023', 'Bergland2014', 'gen', 'trt', 'pop', 'SourceFile']
        writer.writerow(header)
        
        for file in files:
            filename = os.path.splitext(os.path.basename(file))[0]
            if filename in metadata:  # Check if filename is in metadata
                gen = metadata[filename]['gen']
                print(gen)
                pop = metadata[filename]['pop']
                print(pop)
                trt = metadata[filename]['trt']
                print(trt)
            else:
                gen = pop = trt = 'Unknown'  # Default value if not found
                
            with open(file, 'r') as infile:
                for line in infile:
                    #print(line)
                    parts = line.strip().split('\t')
                    #print(parts)
                    if parts and len(parts) > 9:  # Ensure there are enough parts to parse
                        try:
                            # The RD, AD, DP, FREQ values are in the 10th part, after splitting by ':' and then by ','
                            genotype_info = parts[9].split(':')
                            #print(genotype_info)
                            #print(len(genotype_info) >= 5)
                            if len(genotype_info) >= 5:  # Check that we have enough parts after splitting
                                rd, ad, dp, freq = genotype_info[1], genotype_info[2], genotype_info[3], genotype_info[4]
#                                 print(rd)
#                                 print(ad)
#                                 print(dp)
#                                 print(freq)
                                chrom, pos, ref, alt = parts[0], parts[1], parts[3], parts[4]
#                                 print(chrom)
#                                 print(pos)
#                                 print(ref)
#                                 print(alt)
                                row = [chrom, pos, ref, alt, rd, ad, dp, freq] + parts[10:] + [gen, trt, pop, filename]
#                                 print(row)
                                writer.writerow(row)
                        except IndexError as e:
                            print(f"Error processing line in file {file}: {line}. Error: {e}")

def main():
    parser = argparse.ArgumentParser(description='Process and reformat genomic data files.')
    parser.add_argument('files', nargs='+', help='File(s) to process')
    parser.add_argument('--metadata', required=True, help='CSV file containing metadata')
    parser.add_argument('--output', default='output.csv', help='Output file name')
    args = parser.parse_args()

    metadata = read_metadata(args.metadata)
    process_files(args.files, metadata, args.output)
    print(f'Processed files. Output saved to {args.output}')

if __name__ == '__main__':
    main()
