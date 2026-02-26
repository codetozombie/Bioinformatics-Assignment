import os
from src import data_collection
from src import preprocessing
from src import pairwise_alignment
from src import msa
from src import consensus
from Bio import SeqIO
import config

def run_pipeline():
    print("="*50)
    print("DCIT 411: BIOINFORMATICS PIPELINE INITIALIZED")
    print("="*50)
    
    # 1. Fetch Data
    accessions = ["NP_000549.1", "NP_001009042.1", "NP_032246.2", "NP_036706.1"]
    raw_file = "hemoglobin_raw.fasta"
    data_collection.download_sequences(accessions, raw_file)
    print("-" * 50)
    
    # 2. Preprocess Data
    clean_file = "hemoglobin_clean.fasta"
    preprocessing.clean_fasta(raw_file, clean_file)
    print("-" * 50)
    
    # 3. Pairwise Alignment (Using first two sequences)
    clean_path = os.path.join(config.PROCESSED_DATA_DIR, clean_file)
    if os.path.exists(clean_path):
        records = list(SeqIO.parse(clean_path, "fasta"))
        if len(records) >= 2:
            pairwise_alignment.perform_alignment(records[0], records[2], mode='global')
            pairwise_alignment.perform_alignment(records[0], records[2], mode='local')
    print("-" * 50)
    
    # 4. Multiple Sequence Alignment (MSA) using MUSCLE
    # Ensure you have 'muscle' (or 'muscle.exe' on Windows) in your root folder
    executable_path = "muscle" if os.name != 'nt' else "muscle.exe"
    msa.perform_msa_muscle(clean_file, executable_path)
    print("-" * 50)
    
    # 5. Advanced Topic: Consensus Generation
    consensus.generate_consensus("hemoglobin_msa.fasta")
    print("=" * 50)
    print("PIPELINE EXECUTION COMPLETE.")
    print("=" * 50)

if __name__ == "__main__":
    run_pipeline()