from Bio.Seq import Seq

START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}

def find_orfs_in_sequence(seq: str):
    orfs = []

    for frame in range(3):
        orfs.extend(find_orfs_in_frame(seq, frame, strand="+"))

    rev_seq = str(Seq(seq).reverse_complement())
    for frame in range(3):
        orfs.extend(find_orfs_in_frame(rev_seq, frame, strand="-"))

    return orfs


def find_orfs_in_frame(seq: str, frame: int, strand: str):
    orfs = []
    seq_len = len(seq)
    codons = [seq[i:i+3] for i in range(frame, seq_len, 3) if len(seq[i:i+3]) == 3]

    start_positions = []
    for i, codon in enumerate(codons):
        if codon == START_CODON and not start_positions:
            start_positions.append(i)
        elif codon in STOP_CODONS:
            for start_i in start_positions:
                start = frame + start_i * 3
                stop = frame + i * 3 + 3
                orf_seq = seq[start:stop]
                orfs.append((frame, strand, start, stop, orf_seq))
                start_positions = []
    return orfs