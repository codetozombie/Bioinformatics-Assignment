import os
import sys
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

# Add the parent directory to sys.path to import config.py
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import config

def perform_alignment(seq_record1, seq_record2, mode='global'):
    """
    Performs pairwise sequence alignment using dynamic programming.
    mode='global' applies Needleman-Wunsch.
    mode='local' applies Smith-Waterman.
    """
    aligner = PairwiseAligner()
    aligner.mode = mode
    
    # Milestone Requirement: Use a substitution matrix (BLOSUM62 for proteins)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    
    # Milestone Requirement: Experiment with gap penalties
    # An affine gap penalty: higher cost to open a gap, lower cost to extend it
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    
    print(f"\n--- Running {mode.capitalize()} Alignment ---")
    print(f"Aligning: {seq_record1.id} vs {seq_record2.id}")
    
    # Perform the alignment
    alignments = aligner.align(seq_record1.seq, seq_record2.seq)
    
    # The aligner returns a generator of optimal alignments. We take the first one.
    best_alignment = alignments[0]
    
    print(f"Alignment Score: {best_alignment.score}")
    
    # Print a snippet of the actual alignment string (first 80 characters for readability)
    print("Alignment Preview:")
    alignment_lines = str(best_alignment).split('\n')
    for line in alignment_lines:
        print(line[:80])
        
    return best_alignment

if __name__ == "__main__":
    input_path = os.path.join(config.PROCESSED_DATA_DIR, "hemoglobin_clean.fasta")
    
    try:
        # Load the cleaned sequences
        records = list(SeqIO.parse(input_path, "fasta"))
        
        if len(records) >= 2:
            # Let's compare Human (index 0) and Mouse (index 2) Hemoglobin
            human_seq = records[0]
            mouse_seq = records[2]
            
            # Run Needleman-Wunsch
            
            global_aln = perform_alignment(human_seq, mouse_seq, mode='global')
            
            # Run Smith-Waterman
            
            local_aln = perform_alignment(human_seq, mouse_seq, mode='local')
            
            # Save the global alignment output to the results folder
            output_file = os.path.join(config.ALIGNMENTS_DIR, "human_vs_mouse_global.txt")
            with open(output_file, "w") as f:
                f.write(str(global_aln))
            print(f"\nGlobal alignment saved to {output_file}")
            
        else:
            print("Not enough sequences found for pairwise alignment. Check preprocessing.")
            
    except FileNotFoundError:
        print(f"Error: {input_path} not found. Please run preprocessing.py first.")