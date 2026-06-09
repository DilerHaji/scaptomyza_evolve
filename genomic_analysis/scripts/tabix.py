import argparse
import subprocess
import os

def split_variants_file(variants_file, chunks):
    with open(variants_file, 'r') as f:
        lines = f.readlines()
    
    chunk_size = len(lines) // chunks
    chunks = [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]
    
    chunk_files = []
    for i, chunk in enumerate(chunks):
        chunk_file = f"{variants_file}.chunk_{i}"
        with open(chunk_file, 'w') as f:
            f.writelines(chunk)
        chunk_files.append(chunk_file)
    
    return chunk_files

def process_chunk(input_vcf, variants_chunk, output_prefix, chunk_index):
    output_file = f"{output_prefix}_chunk_{chunk_index}.vcf"
    command = f"tabix -f -R {variants_chunk} {input_vcf} > {output_file}"
    subprocess.run(command, shell=True, check=True)
    return output_file

def concatenate_vcfs(chunk_files, output_file):
    with open(output_file, 'w') as outfile:
        for chunk_file in chunk_files:
            with open(chunk_file, 'r') as infile:
                outfile.write(infile.read())

def main(input_vcf, variants_file, output_prefix, chunks):
    # Split variants file
    chunk_files = split_variants_file(variants_file, chunks)
    
    # Process each chunk
    output_chunks = []
    for i, chunk_file in enumerate(chunk_files):
        output_chunk = process_chunk(input_vcf, chunk_file, output_prefix, i)
        output_chunks.append(output_chunk)
    
    # Concatenate results
    final_output = f"{output_prefix}_body.vcf"
    concatenate_vcfs(output_chunks, final_output)
    
    # Clean up temporary files
    for file in chunk_files + output_chunks:
        os.remove(file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF file in chunks using tabix")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("variants_file", help="Variants file for filtering")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--chunks", type=int, default=10, help="Number of chunks to split the variants file into")
    
    args = parser.parse_args()
    
    main(args.input_vcf, args.variants_file, args.output_prefix, args.chunks)