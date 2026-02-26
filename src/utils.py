import os
import sys
from Bio import SeqIO

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import config

def load_fasta(filepath):
    """Loads a FASTA file and returns a list of SeqRecord objects."""
    if not os.path.exists(filepath):
        print(f"Error: {filepath} not found.")
        return []
    return list(SeqIO.parse(filepath, "fasta"))

def save_score_to_file(filename, data_string):
    """Saves text/metrics to the results directory."""
    path = os.path.join(config.METRICS_DIR, filename)
    with open(path, "w") as f:
        f.write(data_string)
    print(f"Metrics saved to {path}")