import argparse
import time
import numpy as np
import plotly.graph_objects as go
from scipy.stats import norm, skewnorm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random
import os
import pandas as pd
import csv
from tabulate import tabulate

def print_intro():
    print("\n")
    print("-- Generate and visualise synthetic reads --")
    print("This script generates synthetic reads from mutliple multi-fasta files under a specified length distribution and introduced error rates.")
    print("\n")

def print_settings_info(args):
    print(f"-- loading parsed arguments:")
    print(f"  I/O")
    print(f"    > input path: {args.input_directory}")
    print(f"    > output path: {args.output_directory}")
    print(f"    > tsv file: {args.tsv_file}")
    print(f"    > prefix: {args.prefix}")
    print(f"  read length distribution")
    print(f"    > peak of size distribution: {args.peak}")
    print(f"    > SD of size Distribution: {args.std_dev}")
    print(f"    > Skewness: {args.skewness}")
    print(f"  error rate distribution")
    print(f"    > mean error rate: {args.error_mean}")
    print(f"    > SD error rate: {args.error_std_dev}")
    print("\n")
    
def generate_perfect_distribution(peak, std_dev, num_points=1000):
    x_vals = np.linspace(peak - 4 * std_dev, peak + 4 * std_dev, num_points)
    y_vals = norm.pdf(x_vals, loc=peak, scale=std_dev)
    y_vals /= np.sum(y_vals)
    return x_vals, y_vals

def generate_perfect_skewed_distribution(peak, std_dev, skewness, num_points=1000):
    print(f"-- calculating perfect skewed gaussian distrubution for read lengths")
    x_vals = np.linspace(peak - 4 * std_dev, peak + 4 * std_dev, num_points)
    y_vals = skewnorm.pdf(x_vals, a=skewness, loc=peak, scale=std_dev)
    y_vals /= np.sum(y_vals)
    return x_vals, y_vals

def generate_random_distribution(x_vals, y_vals, num_reads):
    random_sizes = np.random.choice(x_vals, size=num_reads, p=y_vals)
    return random_sizes

def generate_error_rates(num_reads, mean, std_dev):
    print(f"-- calculating error rate distribution")
    error_rates = np.random.normal(mean, std_dev, num_reads)
    error_rates = np.clip(error_rates, 0, 1)
    return error_rates

def visualize_distribution(fragment_lengths, perfect_distribution, output_directory, output_file_prefix, prefix, error_mean, error_std_dev):
    print(f"-- visualising length distribution of extracted fragments")
    hist_trace = go.Histogram(x=fragment_lengths, nbinsx=150, name='extracted fragments', opacity=0.75, marker=dict(color='grey', line=dict(color='white', width=1)))
    line_trace = go.Scatter(x=perfect_distribution[0], y=perfect_distribution[1], mode='lines', line=dict(color='rgba(0,0,0,0)'), name='Perfect Distribution', showlegend=False)
    subplot = go.Figure([go.Scatter(x=perfect_distribution[0], y=perfect_distribution[1], mode='lines', name='Perfect Distribution')])

    subplot.update_layout(title='Perfect Distribution', xaxis=dict(title='Read Size'), yaxis=dict(title='Count'), yaxis2=dict(title='Count', overlaying='y', side='right'), showlegend=True)

    layout = go.Layout(
        title='Length Distribution',
        xaxis=dict(title='Read Length (bp)', showline=True, linewidth=1, linecolor='black', ticks='outside'),
        yaxis=dict(title='Count', showline=True, linewidth=1, linecolor='black',ticks='outside', rangemode="tozero"),
        yaxis2=dict(title='Extracted Fragments', overlaying='y', side='right', showline=True, linewidth=1, linecolor='black',ticks='outside',rangemode="tozero"),
        showlegend=True,
        paper_bgcolor='white',
        plot_bgcolor='white',
        legend=dict(orientation="h", y=1.15)
    )

    fig = go.Figure(data=[hist_trace, line_trace], layout=layout)
    fig.add_trace(go.Scatter(x=subplot.data[0].x, y=subplot.data[0].y, mode='lines', name='Perfect Distribution', yaxis='y2', line=dict(color='orange')))
    fig.update_layout(xaxis2=dict(title='Read Length (bp)', overlaying='x', side='bottom'), yaxis2=dict(title='Count', overlaying='y', side='right'))

    output_image_file = os.path.join(output_directory, f"{prefix}.{error_mean}_{error_std_dev}.{output_file_prefix}_size_distribution.png")
    fig.write_image(output_image_file)
    print(f"    size distribution plot saved as {output_image_file}")

def visualize_read_quality(fragments_with_errors, output_directory, prefix, error_mean, error_std_dev):
    print(f"-- visualising introduced error rates vs read length")
    lengths = [frag[0] for frag in fragments_with_errors]
    error_rates = [frag[1] for frag in fragments_with_errors]
    q_values = [-10 * np.log10(er) if er > 0 else 50 for er in error_rates]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=lengths, y=error_rates, mode='markers', name='Error Rate'))
    fig.update_layout(title='Read Length vs Error Rate', xaxis_title='Read Length (bp)', yaxis_title='Error Rate')
    output_image_file = os.path.join(output_directory, f"{prefix}.{error_mean}_{error_std_dev}.read_quality_error_rate.png")
    fig.write_image(output_image_file)
    print(f"    read quality plot (error rate) saved as {output_image_file}")

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=lengths, y=q_values, mode='markers', name='Q Value'))
    fig.update_layout(title='Read Length vs Q Value', xaxis_title='Read Length (bp)', yaxis_title='Q Value')
    output_image_file = os.path.join(output_directory, f"{prefix}.{error_mean}_{error_std_dev}.read_quality_q_value.png")
    fig.write_image(output_image_file)
    print(f"    read quality plot (Q value) saved as {output_image_file}")

def read_tsv_file(tsv_file):
    df = pd.read_csv(tsv_file, sep='\t', skipinitialspace=True)
    return df

def introduce_errors(sequence, error_rate):
    bases = ['A', 'C', 'G', 'T']
    num_errors = int(len(sequence) * error_rate)
    error_positions = random.sample(range(len(sequence)), num_errors)

    sequence_list = list(sequence)
    for pos in error_positions:
        current_base = sequence_list[pos]
        new_base = random.choice([base for base in bases if base != current_base])
        sequence_list[pos] = new_base

    return ''.join(sequence_list)


def process_fasta_file(fasta_file, is_circular):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    if is_circular and len(sequences) == 1:
        sequence = sequences[0].seq
        concatenated_sequence = sequence + sequence
        concatenated_record = SeqRecord(concatenated_sequence, id=sequences[0].id, description=sequences[0].description + " (concatenated)")
        return [concatenated_record]
    return sequences

#def calculate_required_reads(seq_length, read_length, coverage):
#    return int((coverage * seq_length) / read_length)

def calculate_required_reads(seq_length, read_length, coverage):
    return int(np.ceil((coverage * seq_length) / read_length))
    
def extract_fragment_from_circular(seq_record, size):
    """Extract a fragment of `size` from a circular sequence.
    Returns the fragment, and normalized start/end positions within original sequence length.
    """
    seq_len = len(seq_record.seq) // 2  # because we concatenate once
    start_pos = random.randint(0, seq_len - 1)
    end_pos = (start_pos + size) % seq_len
    fragment = seq_record.seq[start_pos:start_pos + size]
    normalized_start = start_pos % seq_len
    normalized_end = (start_pos + size) % seq_len

    return fragment, normalized_start, normalized_end
    

def extract_fragments(input_files, random_sizes, output_fasta, output_fastq, output_tsv, coverages, is_circular_flags, error_rates, prefix):
    print(f"-- extracting fragments and introducing errors (if any)")
    records = []
    for file, is_circular in zip(input_files, is_circular_flags):
        records.extend(process_fasta_file(file, is_circular))

    summary_data = []
    read_count = 0
    fragments_with_errors = []

    with open(output_fasta, "w") as output_handle, open(output_fastq, "w") as fastq_handle, open(output_tsv, "w", newline='') as tsv_handle:
        tsv_writer = csv.writer(tsv_handle, delimiter='\t')
        tsv_writer.writerow(["sequence_id", "origin", "contig", "length", "start_position", "end_position", "error_rate", "q_value"])

        for file, coverage in zip(input_files, coverages):
            file_records = process_fasta_file(file, is_circular_flags[input_files.index(file)])
            reads_for_file = 0
            for seq_record in file_records:
                #num_reads = calculate_required_reads(len(seq_record.seq) // 3, np.mean(random_sizes), coverage)
                num_reads = calculate_required_reads(len(seq_record.seq), np.mean(random_sizes), coverage)
                for _ in range(num_reads):
                    if read_count >= len(random_sizes):
                        break
                    size = int(random.choice(random_sizes))
                    #if is_circular_flags[input_files.index(file)] and size > len(seq_record.seq) // 3:
                    if is_circular_flags[input_files.index(file)] and size > len(seq_record.seq):
                        fragment, start, end = extract_fragment_from_circular(seq_record, size)
                    else:
                        # Standard extraction, non-circular
                        if size > len(seq_record.seq):
                            # Skip this read if the requested size is invalid
                            print(f"Warning: skipped fragment (size {size}) > sequence length ({len(seq_record.seq)}) in {seq_record.id}")
                            continue
                        start = random.randint(0, len(seq_record.seq) - size)
                        fragment = seq_record.seq[start:start + size]
                        end = (start + size) % len(seq_record.seq)
                        #end = (start + size) % (len(seq_record.seq) // 3)

                    error_rate = error_rates[read_count]
                    fragment_with_errors = introduce_errors(fragment, error_rate)

                    fragment_record = SeqRecord(Seq(fragment_with_errors), id=f"read_{read_count+1}", description=f"length={size} origin={os.path.basename(file)} contig={seq_record.id}")
                    SeqIO.write(fragment_record, output_handle, "fasta")

                    q_value = -10 * np.log10(error_rate) if error_rate > 0 else 50
                    tsv_writer.writerow([fragment_record.id, os.path.basename(file), seq_record.id, size, start, end, error_rate, q_value])
                    fragments_with_errors.append((size, error_rate))

                    # Convert error rate to Q-score, clip to max 50
                    q_value = -10 * np.log10(error_rate) if error_rate > 0 else 50
                    q_value = min(int(q_value), 50)
                    quality_scores = [q_value] * len(fragment)

                    fastq_record = SeqRecord(
                        Seq(fragment_with_errors),
                        id=fragment_record.id,
                        description="",
                        letter_annotations={"phred_quality": quality_scores}
                    )
                    SeqIO.write(fastq_record, fastq_handle, "fastq")

                    read_count += 1
                    reads_for_file += 1

            summary_data.append([os.path.basename(file), coverage, is_circular_flags[input_files.index(file)], reads_for_file])

    print(f"    generated reads saved as fasta: {output_fasta}")
    print(f"    generated reads saved as fastq: {output_fastq}")
    print(f"    read origin and lengths saved as {output_tsv}")

    print_summary_table(summary_data)

    return fragments_with_errors

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate and visualize synthetic reads from multiple multi-fasta files. Optionally introducing sequencing errors.')
    parser.add_argument('-i', '--input_directory', required=True, help='Input directory containing multi-fasta files')
    parser.add_argument('-o', '--output_directory', required=True, help='Output directory')
    parser.add_argument('--tsv_file', required=True, help='Input TSV file containing fasta file names, coverage levels, and circular element flags')
    parser.add_argument('--peak', type=float, default=13000, help='Peak of the size distribution (gaussian) (default: 13000)')
    parser.add_argument('--std_dev', type=float, default=5000, help='SD of size distribution (gaussian) (default: 1000)')
    parser.add_argument('--skewness', type=float, default=3.5, help='Skew of the size distribution (default: 3.5)')
    parser.add_argument('--error_mean', type=float, default=0.001, help='Mean error rate (default: 0.001)')
    parser.add_argument('--error_std_dev', type=float, default=0.0005, help='SD for error rates (default: 0.0005)')
    parser.add_argument('--prefix', default='output', help='Prefix for output files')
    return parser.parse_args()

def print_summary_table(summary_data):
    headers = ["Fasta File", "Coverage", "Circular", "Reads Extracted"]
    print("\nSummary of Read Extraction:")
    print(tabulate(summary_data, headers=headers, tablefmt="simple"))
    print("\n")

def compute_n50(lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    running_sum = 0
    for length in sorted_lengths:
        running_sum += length
        if running_sum >= total / 2:
            return length


def format_bp(num):
    if num >= 1e6:
        return f"{num:,} ({num / 1e6:.2f} Mb)"
    elif num >= 1e3:
        return f"{num:,} ({num / 1e3:.3f} Kb)"
    else:
        return f"{num:,} bp"
        
                    
def main():
    start_time = time.time()
    args = parse_arguments()
    print_intro()
    print_settings_info(args)
    tsv_data = read_tsv_file(args.tsv_file)

    input_files = [os.path.join(args.input_directory, file) for file in tsv_data['fasta_file'].tolist()]

    coverages = tsv_data['coverage'].tolist()
    is_circular_flags = tsv_data['is_circular'].tolist()

    perfect_distribution = generate_perfect_skewed_distribution(args.peak, args.std_dev, skewness=args.skewness)
    average_read_length = np.mean(perfect_distribution[0])

    # calculate total reads needed based on coverage
    #total_reads_needed = sum([calculate_required_reads(len(rec.seq), average_read_length, cov) for file, cov in zip(input_files, coverages) for rec in process_fasta_file(file, is_circular_flags[input_files.index(file)])])

    all_records = {
        file: process_fasta_file(file, is_circular_flags[input_files.index(file)])
        for file in input_files
    }

    total_reads_needed = sum([
        calculate_required_reads(len(rec.seq), average_read_length, cov)
        for file, cov in zip(input_files, coverages)
        for rec in all_records[file]
    ])

    # generate random distribution and error rates for the total number of reads needed
    random_distribution = generate_random_distribution(perfect_distribution[0], perfect_distribution[1], total_reads_needed)
    error_rates = generate_error_rates(total_reads_needed, args.error_mean, args.error_std_dev)

	# set up output files
    output_fasta = os.path.join(args.output_directory, f"{args.prefix}.{args.error_mean}_{args.error_std_dev}.generated_reads.fasta")
    output_fastq = os.path.join(args.output_directory, f"{args.prefix}.{args.error_mean}_{args.error_std_dev}.generated_reads.fastq")
    output_tsv = os.path.join(args.output_directory, f"{args.prefix}.{args.error_mean}_{args.error_std_dev}.generated_reads.tsv")

    # extract fragments with lengths under given distribution as well as errors
    fragments_with_errors = extract_fragments(input_files, random_distribution, output_fasta,output_fastq , output_tsv, coverages, is_circular_flags, error_rates, args.prefix)

    read_lengths = [length for length, _ in fragments_with_errors]
    total_reads = len(read_lengths)
    total_bases = sum(read_lengths)
    n50 = compute_n50(read_lengths)

    print("\nDataset Summary Statistics:")
    print(f"  Total Reads: {total_reads:,}")
    print(f"  Total Bases: {format_bp(total_bases)}")
    print(f"  N50: {format_bp(n50)}")    
    # Visualise
    #visualize_distribution(random_distribution, perfect_distribution, args.output_directory, "extracted_fragments", args.prefix, args.error_mean, args.error_std_dev)
    #visualize_read_quality(fragments_with_errors, args.output_directory, args.prefix, args.error_mean, args.error_std_dev)
    print("\n--- Done ---")
    end_time = time.time()
    print(f"Total runtime: {end_time - start_time:.2f} seconds\n")

if __name__ == "__main__":
    main()

