from frequency_finder import compute_protein_level_frequencies
from translation import filter_and_translate_orfs
from file_reading import read_fasta_and_find_orfs
from distance_calculator import compute_distance_matrix, save_phylip_matrix
from table_exporter import save_frequency_table
from pathlib import Path

DATA_DIR = Path("data")
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)

def process_genome(fasta_path):

    #1st and 2nd
    orfs = read_fasta_and_find_orfs(fasta_path)

    #3rd and 4th
    translated = filter_and_translate_orfs(orfs, min_len=100)
    print(f"   {len(translated)} good ORFs")

    if not translated:
        return None, None
    
    #5th
    avg_codon_freq, avg_dicodon_freq = compute_protein_level_frequencies(translated)

    return avg_codon_freq, avg_dicodon_freq

if __name__ == "__main__":
    fasta_files = sorted(DATA_DIR.glob("*.fasta"))

    codon_results = {}
    dicodon_results = {}

    for fasta_file in fasta_files:
        genome_name = fasta_file.stem
        print(f"Processing {genome_name}")
        codon_freq, dicodon_freq = process_genome(fasta_file)
        if codon_freq and dicodon_freq:
            codon_results[genome_name] = codon_freq
            dicodon_results[genome_name] = dicodon_freq

    save_frequency_table(codon_results, OUTPUT_DIR / "codon_frequencies.csv")
    save_frequency_table(dicodon_results, OUTPUT_DIR / "dicodon_frequencies.csv")

    #6th
    codon_matrix = compute_distance_matrix(codon_results)
    dicodon_matrix = compute_distance_matrix(dicodon_results)

    save_phylip_matrix(codon_matrix, OUTPUT_DIR / "codon_distance_matrix.txt")
    save_phylip_matrix(dicodon_matrix, OUTPUT_DIR / "dicodon_distance_matrix.txt")

    print("Done. Check dicodon_distance_matrix.txt and codon_distance_matrix.txt")