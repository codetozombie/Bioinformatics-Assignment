import os
import sys
import re
import logging
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Add the parent directory to sys.path to import config.py
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import config

# ── Logging setup (falls back to print if logging fails) ──────────────────
try:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    log = logging.getLogger(__name__)
except Exception:
    log = None  # Fallback to print statements

def _log(msg, level="info"):
    """Internal logger wrapper that falls back to print."""
    if log and hasattr(log, level):
        getattr(log, level)(msg)
    else:
        print(f"[{level.upper()}] {msg}")

# =============================================================
# SECTION 1: Individual Cleaning Functions (Compatible Core)
# =============================================================

def remove_gaps(seq_str: str) -> str:
    """Remove all gap characters ('-' and '.') from a sequence string."""
    cleaned = seq_str.replace("-", "").replace(".", "")
    return cleaned

def to_uppercase(seq_str: str) -> str:
    """Convert sequence to uppercase for consistent processing."""
    return seq_str.upper()

def remove_whitespace(seq_str: str) -> str:
    """Remove all whitespace characters (spaces, tabs, newlines)."""
    return re.sub(r"\s+", "", seq_str)

def handle_ambiguous_dna(seq_str: str, replacement: str = "N") -> str:
    """Replace non-standard DNA characters while keeping IUPAC codes."""
    valid = set("ACGTNRYSWKMBDHV")
    return "".join(ch if ch in valid else replacement for ch in seq_str)

def handle_ambiguous_protein(seq_str: str, replacement: str = "X") -> str:
    """Replace non-standard amino acid characters; remove stop codons (*)."""
    seq_str = seq_str.replace("*", "")  # Remove stop codons first
    valid = set("ACDEFGHIKLMNPQRSTVWYBOUXZ")
    return "".join(ch if ch in valid else replacement for ch in seq_str)

def detect_sequence_type(seq_str: str) -> str:
    """
    Heuristically detect whether a sequence is DNA, RNA, or Protein.
    Returns 'DNA', 'RNA', or 'Protein'.
    """
    seq_upper = seq_str.upper()
    total = len(seq_upper)
    if total == 0:
        return "Protein"  # Default fallback
    
    counts = Counter(seq_upper)
    dna_bases = sum(counts.get(b, 0) for b in "ACGTN")
    rna_bases = sum(counts.get(b, 0) for b in "ACGUN")
    
    if dna_bases / total > 0.85:
        return "DNA"
    if rna_bases / total > 0.85 and counts.get("U", 0) > 0:
        return "RNA"
    return "Protein"

# =============================================================
# SECTION 2: Main Preprocessing Pipeline (Drop-in Compatible)
# =============================================================

def preprocess_sequence(seq_str: str, seq_type: str = "auto") -> str:
    """
    Apply full preprocessing to a sequence string.
    Maintains compatibility with simple version's output format.
    """
    # Step 1 & 2: Normalize
    seq_str = to_uppercase(seq_str)
    seq_str = remove_whitespace(seq_str)
    
    # Step 3: Remove gaps (matches simple behavior)
    seq_str = remove_gaps(seq_str)
    
    # Step 4: Detect type if auto
    if seq_type == "auto":
        seq_type = detect_sequence_type(seq_str)
    
    # Step 5: Handle ambiguous characters based on type
    if seq_type in ("DNA", "RNA"):
        seq_str = handle_ambiguous_dna(seq_str)
    elif seq_type == "Protein":
        seq_str = handle_ambiguous_protein(seq_str)
    
    return seq_str

def clean_fasta(input_filename, output_filename):
    """
    Reads raw FASTA files, cleans the sequences (removes gaps and stop codons), 
    and saves the normalized sequences to the processed directory.
    
    COMPATIBLE SIGNATURE: Same parameters and behavior as simple version.
    """
    input_path = os.path.join(config.RAW_DATA_DIR, input_filename)
    output_path = os.path.join(config.PROCESSED_DATA_DIR, output_filename)
    
    print(f"Reading raw sequences from {input_filename}...")
    
    cleaned_records = []
    try:
        # Ensure output directory exists (robustness improvement)
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        
        for record in SeqIO.parse(input_path, "fasta"):
            # Convert sequence to string and preprocess with auto-detection
            seq_str = str(record.seq)
            seq_type = detect_sequence_type(seq_str)
            cleaned_seq = preprocess_sequence(seq_str, seq_type=seq_type)
            
            # Create a new cleaned SeqRecord (matches simple output format)
            cleaned_record = SeqRecord(
                Seq(cleaned_seq),
                id=record.id,
                description=record.description + " | CLEANED"
            )
            cleaned_records.append(cleaned_record)
            
        # Write the cleaned records to the processed directory
        count = SeqIO.write(cleaned_records, output_path, "fasta")
        print(f"Successfully cleaned and saved {len(cleaned_records)} sequences to:\n{output_path}")
        
    except FileNotFoundError:
        print(f"Error: Could not find {input_path}. Ensure data_collection.py ran successfully.")
    except Exception as e:
        print(f"An error occurred during preprocessing: {e}")

if __name__ == "__main__":
    # Run the cleaning function on the file we downloaded in Step 2
    # COMPATIBLE: Same call signature as simple version
    clean_fasta("hemoglobin_raw.fasta", "hemoglobin_clean.fasta")
