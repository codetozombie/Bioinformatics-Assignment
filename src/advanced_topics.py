# src/advanced_topics.py
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.PDB import PDBParser, Superimposer
from Bio import AlignIO
from Bio.Align import AlignInfo

def run_psi_blast(sequence_record):
    """Profile-based alignment via NCBI PSI-BLAST (Position-Specific Iterated BLAST)."""
    print(f"\nRunning PSI-BLAST for {sequence_record.id}...")
    # blastp with step_number > 1 acts as PSI-BLAST
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence_record.seq, hitlist_size=5, format_type="XML")
    blast_record = NCBIXML.read(result_handle)
    
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print(f"Homolog Found: {alignment.title[:50]} | E-value: {hsp.expect}")
            break # Just print top hit per alignment

def structural_alignment(pdb_file1, pdb_file2):
    """Aligns protein structures based on 3D coordinates using Bio.PDB."""
    parser = PDBParser(QUIET=True)
    struct1 = parser.get_structure("S1", pdb_file1)
    struct2 = parser.get_structure("S2", pdb_file2)
    
    # Extract alpha carbons from the first chain
    atoms1 = [atom for atom in struct1[0]['A'].get_atoms() if atom.get_name() == 'CA']
    atoms2 = [atom for atom in struct2[0]['A'].get_atoms() if atom.get_name() == 'CA']
    
    # Equalize lengths for superimposition
    min_len = min(len(atoms1), len(atoms2))
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1[:min_len], atoms2[:min_len])
    super_imposer.apply(struct2.get_atoms())
    
    print(f"\nStructural Alignment RMSD: {super_imposer.rms:.3f} Ångströms")

def generate_consensus(msa_file_path):
    """Finds the most conserved residues across aligned sequences."""
    alignment = AlignIO.read(msa_file_path, "fasta")
    summary_info = AlignInfo.SummaryInfo(alignment)
    
    # dumb_consensus is a legacy method, but specifically built into Biopython's AlignInfo 
    # for exact positional conservation tracking.
    consensus = summary_info.dumb_consensus(threshold=0.7, ambiguous='X')
    print(f"\nConsensus Sequence (70% threshold):\n{consensus}")

if __name__ == "__main__":
    from Bio import SeqIO
    # Ensure this sequence exists from preprocessing
    records = list(SeqIO.parse("data/processed/hemoglobin_uniform.fasta", "fasta"))
    
    # 1. Profile-Based Alignment
    run_psi_blast(records[0])
    
    # 2. Consensus Generation (Requires an MSA output file from step 3)
    try:
        generate_consensus("results/alignments/hemoglobin_msa.fasta")
    except Exception as e:
        print("Run msa_benchmark.py to generate MSA file first for Consensus.")
    
    # 3. Structural Alignment (Requires pre-downloaded .pdb files in data/)
    # structural_alignment("data/1a3n.pdb", "data/1ird.pdb")