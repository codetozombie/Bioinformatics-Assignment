import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import config
from src.utils import save_score_to_file

def calculate_identity(aligned_seq1, aligned_seq2):
    """
    Calculates the percentage identity between two aligned sequences of the same length.
    """
    if len(aligned_seq1) != len(aligned_seq2):
        raise ValueError("Aligned sequences must be of the same length.")
    
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
    length = len(aligned_seq1)
    
    identity = (matches / length) * 100
    return round(identity, 2)

if __name__ == "__main__":
    # Example usage (assuming alignments are loaded)
    seq1 = "MVHLTPEEKSAVTALWGKVNV--DEVGGEALGRLL"
    seq2 = "MVHLTDAEKAAVNGLWGKVNPDDPEVGGEALGRLL"
    ident = calculate_identity(seq1, seq2)
    print(f"Sequence Identity: {ident}%")
    save_score_to_file("identity_metrics.txt", f"Identity: {ident}%\n")