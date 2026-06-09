import sys

def extract_chromosome_names_and_sizes(fasta_file):
    chromosome_sizes = {}
    current_chrom = None
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                current_chrom = line[1:].strip()
                chromosome_sizes[current_chrom] = 0
            else:
                if current_chrom:
                    chromosome_sizes[current_chrom] += len(line.strip())
    return chromosome_sizes

def get_top_n_chromosomes(chromosome_sizes, n):
    return sorted(chromosome_sizes, key=chromosome_sizes.get, reverse=True)[:n]

def write_to_file(file_path, data):
    with open(file_path, 'w') as file:
        for item in data:
            file.write(f"{item}\n")

if __name__ == "__main__":

    fasta_file_path = sys.argv[1]
    output_file_path = sys.argv[3]
    topn = sys.argv[2]
    print(topn)

    chromosome_sizes = extract_chromosome_names_and_sizes(fasta_file_path)
    top_chromosomes = get_top_n_chromosomes(chromosome_sizes, n=int(topn))
    write_to_file(output_file_path, top_chromosomes)
