import os
import sys
import subprocess
from Bio import AlignIO

# Add the parent directory to sys.path to import config.py
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import config

def perform_msa_muscle(input_filename, muscle_exe="muscle"):
    """
    Executes MUSCLE for Multiple Sequence Alignment using Python's subprocess.
    
    Args:
        input_filename (str): The name of the cleaned FASTA file.
        muscle_exe (str): The path to the MUSCLE executable. 
    """
    input_path = os.path.join(config.PROCESSED_DATA_DIR, input_filename)
    
    # MUSCLE typically outputs in FASTA format, which Biopython parses easily
    output_path = os.path.join(config.ALIGNMENTS_DIR, "hemoglobin_msa.fasta")
    
    if not os.path.exists(input_path):
        print(f"Error: Input file {input_path} not found. Run preprocessing.py first.")
        return None

    print(f"Initiating MUSCLE Multiple Sequence Alignment on {input_filename}...")
    
    # Construct the command-line arguments for MUSCLE v5
    # If using an older MUSCLE v3, the flags are "-in" and "-out" instead of "-align" and "-output"
    command = [
        muscle_exe, 
        "-align", input_path, 
        "-output", output_path
    ]
    
    print(f"Executing system command: {' '.join(command)}")
    
    try:
        # subprocess.run is the modern standard for executing external CLI tools in Python
        result = subprocess.run(
            command, 
            capture_output=True, 
            text=True, 
            check=True
        )
        
        if os.path.exists(output_path):
            print(f"\nMSA successfully completed. Alignment saved to:\n{output_path}")
            
            
            
            # Parse the alignment back into Python
            alignment = AlignIO.read(output_path, "fasta")
            
            print(f"\n--- MSA Summary ---")
            print(f"Number of sequences aligned: {len(alignment)}")
            print(f"Alignment length (including gaps): {alignment.get_alignment_length()}")
            
            # Print a quick preview of the first sequence's alignment
            print(f"Preview ({alignment[0].id}): {str(alignment[0].seq)[:50]}...")
            
            return alignment
            
    except subprocess.CalledProcessError as e:
        print(f"\nCRITICAL ERROR: MUSCLE execution failed.")
        print(f"Standard Error Output:\n{e.stderr}")
        print("\nTroubleshooting: Ensure the MUSCLE executable is in your project folder or system PATH, and check if you are using MUSCLE v3 or v5 (flags differ).")
        return None
    except FileNotFoundError:
        print(f"\nCRITICAL ERROR: MUSCLE executable '{muscle_exe}' not found.")
        print("Please download MUSCLE and place the executable in your project root.")
        return None

if __name__ == "__main__":
    # If on Windows, you might need to change this to "muscle.exe"
    # If you put it in the root folder, use "./muscle" (Linux/Mac) or "muscle.exe" (Windows)
    executable_path = "muscle.exe" 
    
    msa_result = perform_msa_muscle("hemoglobin_clean.fasta", executable_path)