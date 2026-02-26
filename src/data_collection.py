import os
import sys
from Bio import Entrez, SeqIO

# Add the parent directory to sys.path so we can import config.py
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import config

def download_sequences(accession_list, output_filename):
    """
    Fetches biological sequences from the NCBI Protein database.
    """
    # Set the email from your config to comply with NCBI guidelines
    Entrez.email = config.ENTREZ_EMAIL
    
    output_path = os.path.join(config.RAW_DATA_DIR, output_filename)
    print(f"Fetching {len(accession_list)} sequences from NCBI...")
    
    try:
        # efetch queries the database and returns the data in FASTA format
        handle = Entrez.efetch(
            db="protein", 
            id=accession_list, 
            rettype="fasta", 
            retmode="text"
        )
        
        # Parse the returned data into SeqRecord objects
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        
        # Write the records to your data/raw/ directory
        SeqIO.write(records, output_path, "fasta")
        print(f"Successfully saved {len(records)} sequences to:\n{output_path}")
        
        return output_path
        
    except Exception as e:
        print(f"An error occurred during data collection: {e}")
        return None

if __name__ == "__main__":
    # Accession numbers for Hemoglobin alpha subunit across four species
    # NP_000549.1 (Human), NP_001009042.1 (Chimp), NP_032246.2 (Mouse), NP_036706.1 (Rat)
    hemoglobin_accessions = [
        "NP_000549.1", 
        "NP_001009042.1", 
        "NP_032246.2", 
        "NP_036706.1"
    ]
    
    # Run the download function
    download_sequences(hemoglobin_accessions, "hemoglobin_raw.fasta")