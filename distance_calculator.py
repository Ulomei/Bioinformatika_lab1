import math
import itertools

def euclidean_distance(freqs1: dict, freqs2: dict) -> float:
    keys = sorted(freqs1.keys())
    return math.sqrt(sum((freqs1[k] - freqs2[k]) ** 2 for k in keys))

def kl_divergence(p: list, q: list) -> float:
    epsilon = 1e-12 
    return sum(pi * math.log((pi + epsilon) / (qi + epsilon)) for pi, qi in zip(p, q))

def jensen_shannon_distance(freqs1: dict, freqs2: dict) -> float:
    keys = sorted(freqs1.keys())
    p = [freqs1[k] for k in keys]
    q = [freqs2[k] for k in keys]

    m = [(pi + qi) / 2 for pi, qi in zip(p, q)]

    jsd = (kl_divergence(p, m) + kl_divergence(q, m)) / 2
    return math.sqrt(jsd)

def compute_distance_matrix(freq_dicts: dict) -> dict:
    genomes = list(freq_dicts.keys())
    matrix = {g1: {} for g1 in genomes}

    for g1, g2 in itertools.combinations(genomes, 2):
        dist = jensen_shannon_distance(freq_dicts[g1], freq_dicts[g2]) #Change the dunction to euclidean_distance() to check how the matrix looks using euclidean distance.
        matrix[g1][g2] = dist
        matrix[g2][g1] = dist

    for g in genomes:
        matrix[g][g] = 0.0

    return matrix

def save_phylip_matrix(matrix: dict, output_file: str):
    genomes = list(matrix.keys())
    n = len(genomes)
    with open(output_file, "w") as f:
        f.write(f"{n}\n")
        for g in genomes:
            row = " ".join(f"{matrix[g][g2]:.3f}" for g2 in genomes)
            f.write(f"{g[:10]:<10} {row}\n")