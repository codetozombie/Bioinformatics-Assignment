import os
import sys
import urllib.request
from Bio.PDB import PDBParser, Superimposer

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import config

def download_pdb(pdb_id, output_dir):
    """Downloads a PDB file from the RCSB Protein Data Bank."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    filepath = os.path.join(output_dir, f"{pdb_id}.pdb")
    if not os.path.exists(filepath):
        print(f"Downloading {pdb_id}...")
        urllib.request.urlretrieve(url, filepath)
    return filepath

def align_structures(pdb_id1, pdb_id2):
    """Aligns two 3D protein structures and calculates the RMSD."""
    # Ensure structures directory exists
    struct_dir = config.STRUCTURES_DIR
    
    file1 = download_pdb(pdb_id1, struct_dir)
    file2 = download_pdb(pdb_id2, struct_dir)
    
    parser = PDBParser(QUIET=True)
    structure1 = parser.get_structure(pdb_id1, file1)
    structure2 = parser.get_structure(pdb_id2, file2)
    
    # Extract alpha carbon (CA) atoms from the first model and first chain
    atoms1 = [atom for atom in structure1[0]['A'].get_atoms() if atom.get_name() == 'CA']
    atoms2 = [atom for atom in structure2[0]['A'].get_atoms() if atom.get_name() == 'CA']
    
    # Ensure we only align matching lengths
    min_len = min(len(atoms1), len(atoms2))
    atoms1 = atoms1[:min_len]
    atoms2 = atoms2[:min_len]
    
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    super_imposer.apply(structure2.get_atoms())
    
    print(f"--- 3D Structural Alignment ---")
    print(f"PDB {pdb_id1} vs {pdb_id2}")
    print(f"RMSD: {round(super_imposer.rms, 3)} Ångströms")
    return super_imposer.rms

if __name__ == "__main__":
    # Example: Human Hemoglobin vs. Mouse Hemoglobin structures
    align_structures("1A3N", "1IRD")