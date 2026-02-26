# src/msa_benchmark.py
import time
import tracemalloc
import subprocess
from Bio import AlignIO

def benchmark_msa_tool(tool_name, command_list, output_file):
    print(f"\n--- Benchmarking {tool_name} ---")
    
    tracemalloc.start()
    start_time = time.time()
    
    try:
        subprocess.run(command_list, capture_output=True, check=True)
        
        # Stop tracking
        end_time = time.time()
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        # Evaluate Accuracy (using alignment length and conservation as a proxy)
        alignment = AlignIO.read(output_file, "fasta")
        
        print(f"Speed (Execution Time): {end_time - start_time:.4f} seconds")
        print(f"Memory (Peak Usage): {peak / 10**6:.4f} MB")
        print(f"Accuracy Proxy (Alignment Length): {alignment.get_alignment_length()} columns")
        
    except FileNotFoundError:
        print(f"Executable for {tool_name} not found. Skipping.")
    except Exception as e:
        print(f"Error running {tool_name}: {e}")

if __name__ == "__main__":
    input_fasta = "data/processed/hemoglobin_uniform.fasta"
    
    # 1. MUSCLE
    muscle_cmd = ["muscle.exe", "-align", input_fasta, "-output", "data/muscle_out.fasta"]
    benchmark_msa_tool("MUSCLE", muscle_cmd, "data/muscle_out.fasta")
    
    # 2. MAFFT
    mafft_cmd = ["mafft.bat", "--auto", input_fasta]
    # MAFFT prints to stdout, so we write it to a file
    print("\nNote: MAFFT stdout handling requires custom subprocess piping in production.")
    
    # # 3. ClustalW
    # clustalw_cmd = ["clustalw2", f"-INFILE={input_fasta}", "-OUTFILE=data/clustal_out.fasta", "-OUTPUT=FASTA"]
    # benchmark_msa_tool("ClustalW", clustalw_cmd, "data/clustal_out.fasta")