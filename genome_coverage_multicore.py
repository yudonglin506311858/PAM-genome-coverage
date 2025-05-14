import argparse
import re
from Bio.Seq import Seq
from multiprocessing import Pool, cpu_count
import os
import sys
from collections import defaultdict

# IUPAC nucleotide code mapping
IUPAC_CODES = {
    'R': '[AG]',
    'Y': '[CT]',
    'M': '[AC]',
    'K': '[GT]',
    'S': '[GC]',
    'W': '[AT]',
    'H': '[ATC]',
    'B': '[GTC]',
    'V': '[GAC]',
    'D': '[GAT]',
    'N': '[ACGT]',
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T'
}

def pattern_to_regex(pattern):
    """Convert a pattern with IUPAC codes to a regular expression."""
    regex = []
    for char in pattern.upper():
        if char in IUPAC_CODES:
            regex.append(IUPAC_CODES[char])
        else:
            raise ValueError(f"Invalid IUPAC code: {char}")
    return re.compile(''.join(regex))

def get_reverse_complement(seq):
    """Get reverse complement of a DNA sequence."""
    return str(Seq(seq).reverse_complement())

def generate_all_possible_sequences(pattern):
    """Generate all possible DNA sequences that match the IUPAC pattern."""
    from itertools import product
    
    bases_map = {
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T']
    }
    
    # Get all possible base combinations for the pattern
    possible_bases = [bases_map[char] for char in pattern]
    # Generate all combinations
    return [''.join(comb) for comb in product(*possible_bases)]

def process_chromosome(args):
    """Process a single chromosome sequence."""
    chrom_name, sequence, pattern = args
    sequence = sequence.upper()
    pattern_len = len(pattern)
    matches = []
    
    # Generate all possible forward sequences that match the pattern
    possible_sequences = generate_all_possible_sequences(pattern)
    
    # Generate all possible reverse complement sequences
    rc_possible_sequences = {get_reverse_complement(s) for s in possible_sequences}
    
    # Check every possible position in the sequence
    for i in range(len(sequence) - pattern_len + 1):
        window = sequence[i:i+pattern_len]
        
        # Check forward strand
        if window in possible_sequences:
            matches.append((window, chrom_name, i, i+pattern_len, '+'))
        
        # Check reverse strand
        if window in rc_possible_sequences:
            rc_seq = get_reverse_complement(window)
            matches.append((rc_seq, chrom_name, i, i+pattern_len, '-'))
    
    return matches, len(sequence)

def split_genome(file_path):
    """Split genome file into chromosomes for parallel processing."""
    chromosomes = []
    current_chr = None
    current_seq = []
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_chr is not None:
                    chromosomes.append((current_chr, ''.join(current_seq)))
                current_chr = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
    
        if current_chr is not None:
            chromosomes.append((current_chr, ''.join(current_seq)))
    
    return chromosomes

def count_sequences(genome_file, pattern, threads):
    """Count sequences using multiple threads."""
    chromosomes = split_genome(genome_file)
    total_bases = 0
    all_matches = []
    
    # Prepare arguments for parallel processing
    tasks = [(chrom, seq, pattern) for chrom, seq in chromosomes]
    
    with Pool(threads) as pool:
        results = pool.map(process_chromosome, tasks)
    
    for matches, chrom_len in results:
        all_matches.extend(matches)
        total_bases += chrom_len
    
    return all_matches, total_bases

def write_matches(output_file, matches):
    """Write matches to output file."""
    with open(output_file, 'w') as out:
        # Write header
        out.write("Sequence\tChromosome\tStart\tEnd\tStrand\n")
        
        # Write matches
        for match in matches:
            seq, chrom, start, end, strand = match
            out.write(f"{seq}\t{chrom}\t{start}\t{end}\t{strand}\n")

def append_statistics(output_file, matches, total_bases, pattern_len):
    """Append statistics to the output file."""
    # Calculate unique match positions (considering both position and strand)
    unique_positions = set()
    position_counts = defaultdict(int)
    covered_bases = set()      
    for seq, chrom, start, end, strand in matches:
        unique_positions.add((chrom, start, end, strand))
        position_counts[(chrom, start, end)] += 1
        covered_bases.update(range(start, end + 1))
    
    total_unique_matches = len(unique_positions)
    total_possible = 2 * (total_bases - pattern_len + 1)  # Both strands
    coverage_percentage = (len(unique_positions) / total_bases) * 100 if total_bases > 0 else 0
    
    # Additional statistics
    double_strand_matches = sum(1 for count in position_counts.values() if count > 1)
    
    with open(output_file, 'a') as out:
        out.write("\n# Summary Statistics\n")
        out.write(f"# Total raw matches: {len(matches)}\n")
        out.write(f"# Unique positions: {total_unique_matches}\n")
        out.write(f"# Positions matched on both strands: {double_strand_matches}\n")
        out.write(f"# Total possible positions: {total_possible}\n")
        out.write(f"# Sequence coverage: {coverage_percentage:.6f}%\n")  
        out.write(f"# Pattern length: {pattern_len}\n")
        out.write(f"# Total genome bases: {total_bases}\n")

def main():
    parser = argparse.ArgumentParser(description='Count specific sequences in a genome.')
    parser.add_argument('-genome', required=True, help='Path to the genome FASTA file')
    parser.add_argument('-pam', required=True, help='PAM sequence to search for (IUPAC codes allowed)')
    parser.add_argument('-o', required=True, help='Output file path')
    parser.add_argument('-t', '--threads', type=int, default=max(1, cpu_count() - 1),
                       help='Number of threads to use (default: all available -1)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.o), exist_ok=True) if os.path.dirname(args.o) else None
    
    # Process genome
    matches, total_bases = count_sequences(args.genome, args.pam, args.threads)
    
    # Write output
    write_matches(args.o, matches)
    append_statistics(args.o, matches, total_bases, len(args.pam))
    
    print(f"Processing complete. Results saved to {args.o}", file=sys.stderr)

if __name__ == '__main__':
    main()
