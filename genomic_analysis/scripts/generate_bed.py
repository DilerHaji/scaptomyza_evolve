import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ref', type=str, required=True, help='Path reference fasta')
parser.add_argument('--output', type=str, required=True, help='name of output bed file')
parser.add_argument('--chunk_size', type=str, required=True, help='chunk size in bases')
args = parser.parse_args()


def generate_bed_file(ref, output, chunk_size):
    from Bio import SeqIO
    
    chunk_size = int(chunk_size)

    with open(output, 'w') as bed_file:
        for record in SeqIO.parse(ref, "fasta"):
            seq_length = len(record.seq)
            for start in range(0, seq_length, chunk_size):
                end = min(start + chunk_size, seq_length)
                # Adjust the end point to be non-inclusive
                if end != seq_length:
                    end -= 1
                bed_file.write(f"{record.id}\t{start}\t{end}\n")


generate_bed_file(args.ref, args.output, args.chunk_size)