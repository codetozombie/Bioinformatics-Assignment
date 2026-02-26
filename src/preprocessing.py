import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Add the parent directory to sys.path to import config.py
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import config

def clean_fasta(input_filename, output_filename):
    """
    Reads raw FASTA files, cleans the sequences (removes gaps and stop codons), 
    and saves the normalized sequences to the processed directory.
    """
    input_path = os.path.join(config.RAW_DATA_DIR, input_filename)
    output_path = os.path.join(config.PROCESSED_DATA_DIR, output_filename)
    
    print(f"Reading raw sequences from {input_filename}...")
    
    cleaned_records = []
    try:
        for record in SeqIO.parse(input_path, "fasta"):
            # Convert sequence to string and ensure uppercase
            seq_str = str(record.seq).upper()
            
            # Remove trailing stop codons (*) and any pre-existing gaps (-)
            seq_str = seq_str.replace("*", "").replace("-", "")
            
            # Create a new cleaned SeqRecord
            cleaned_record = SeqRecord(
                Seq(seq_str),
                id=record.id,
                description=record.description + " | CLEANED"
            )
            cleaned_records.append(cleaned_record)
            
        # Write the cleaned records to the processed directory
        SeqIO.write(cleaned_records, output_path, "fasta")
        print(f"Successfully cleaned and saved {len(cleaned_records)} sequences to:\n{output_path}")
        
    except FileNotFoundError:
        print(f"Error: Could not find {input_path}. Ensure data_collection.py ran successfully.")
    except Exception as e:
        print(f"An error occurred during preprocessing: {e}")

if __name__ == "__main__":
    # Run the cleaning function on the file we downloaded in Step 2
    clean_fasta("hemoglobin_raw.fasta", "hemoglobin_clean.fasta")