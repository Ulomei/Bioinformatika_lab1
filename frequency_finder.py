from collections import Counter
from itertools import product

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
ALL_DICODONS = ["".join(p) for p in product(AMINO_ACIDS, repeat=2)]

def compute_protein_level_frequencies(translated):
    total_aa_counts = Counter()
    total_dicodon_counts = Counter()
    total_aa = 0
    total_dicodons = 0

    for _, _, _, _, _, protein in translated:
        protein = protein.replace("*", "")
        total_aa_counts.update(protein)
        total_aa += len(protein)

        dicodons = [protein[i:i+2] for i in range(len(protein) - 1)]
        total_dicodon_counts.update(dicodons)
        total_dicodons += max(len(protein) - 1, 0)

    avg_codon_freq = {
        aa: total_aa_counts[aa] / total_aa if total_aa > 0 else 0.0
        for aa in AMINO_ACIDS
    }

    avg_dicodon_freq = {
        pair: total_dicodon_counts[pair] / total_dicodons if total_dicodons > 0 else 0.0
        for pair in ALL_DICODONS
    }

    return avg_codon_freq, avg_dicodon_freq