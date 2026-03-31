import argparse
import os
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tabulate import tabulate
import pandas as pd

def print_intro():
    print("\n-- Assembly Chopper --")
    print("Chop multi-fasta sequences into fixed-length fragments based on TSV config.\n")

def process_fasta_file(fasta_file):
    return list(SeqIO.parse(fasta_file, "fasta"))

def chop_linear(seq_record, fragment_length):
    fragments = []
    seq_len = len(seq_record.seq)
    for start in range(0, seq_len - fragment_length + 1, fragment_length):
        end = start + fragment_length
        fragments.append((seq_record.seq[start:end], start, end))
    return fragments

def chop_circular(seq_record, fragment_length):
    # For simplicity, treat circular like linear, ignore remainder
    return chop_linear(seq_record, fragment_length)

def write_fragments(fragments, output_fasta, output_tsv, parent_file, contig, start_id):
    frag_id = start_id
    with open(output_fasta, "a") as fasta_handle, open(output_tsv, "a", newline='') as tsv_handle:
        tsv_writer = csv.writer(tsv_handle, delimiter='\t')
        if os.stat(output_tsv).st_size == 0:
            tsv_writer.writerow(["fragment_id", "parent_file", "contig", "start", "end", "length"])
        for seq, start, end in fragments:
            record = SeqRecord(seq, id=f"frag_{frag_id}", description="")
            SeqIO.write(record, fasta_handle, "fasta")
            tsv_writer.writerow([record.id, parent_file, contig, start, end, len(seq)])
            frag_id += 1
    return frag_id

def parse_arguments():
    parser = argparse.ArgumentParser(description="Chop assemblies into fixed-length fragments based on TSV config")
    parser.add_argument('-i', '--input_directory', required=True, help='Input directory with fasta files')
    parser.add_argument('-o', '--output_directory', required=True, help='Output directory')
    parser.add_argument('--tsv_file', required=True, help='TSV config file specifying fasta files and circularity')
    parser.add_argument('--fragment_length', type=int, required=True, help='Length of fragments to generate')
    parser.add_argument('--prefix', default='fragments', help='Prefix for output files')
    return parser.parse_args()

def print_summary_table(summary_data):
    headers = ["Fasta File", "Total Fragments", "Fragment Length"]
    print("\nSummary of Chopping:")
    print(tabulate(summary_data, headers=headers, tablefmt="simple"))
    print("\n")

def main():
    args = parse_arguments()
    print_intro()

    os.makedirs(args.output_directory, exist_ok=True)
    output_fasta = os.path.join(args.output_directory, f"{args.prefix}.fragments.fasta")
    output_tsv = os.path.join(args.output_directory, f"{args.prefix}.fragments.tsv")

    open(output_tsv, "w").close()

    tsv_config = pd.read_csv(args.tsv_file, sep='\t', skipinitialspace=True)
    
    frag_id_counter = 1
    summary_data = []

    for idx, row in tsv_config.iterrows():
        fasta_file = os.path.join(args.input_directory, row['fasta_file'])
        is_circular = bool(row['is_circular'])
        parent_name = row['fasta_file']
        
        sequences = process_fasta_file(fasta_file)
        total_fragments = 0

        for seq_record in sequences:
            if is_circular:
                fragments = chop_circular(seq_record, args.fragment_length)
            else:
                fragments = chop_linear(seq_record, args.fragment_length)
            frag_id_counter = write_fragments(fragments, output_fasta, output_tsv, parent_name, seq_record.id, frag_id_counter)
            total_fragments += len(fragments)
        summary_data.append([parent_name, total_fragments, args.fragment_length])

    print_summary_table(summary_data)
    print(f"Fragments written to {output_fasta}")
    print(f"Fragment metadata written to {output_tsv}")
    print("\n--- Done ---")

if __name__ == "__main__":
    main()
