from Bio import SeqIO
from orf_finder import find_orfs_in_sequence

def read_fasta_and_find_orfs(fasta_file: str):
    results = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        orfs = find_orfs_in_sequence(seq)
        results.extend(orfs)
    return results