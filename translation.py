from Bio.Seq import Seq

def filter_and_translate_orfs(orfs, min_len):
    proteins = []
    for frame, strand, start, stop, dna_seq in orfs:
        if len(dna_seq) >= min_len:
            protein = str(Seq(dna_seq).translate(table="Standard"))
            proteins.append((frame, strand, start, stop, dna_seq, protein))
    return proteins