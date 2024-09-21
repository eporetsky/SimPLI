import argparse
import MDAnalysis as mda
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import matplotlib.pyplot as plt
import numpy as np
from io import BytesIO
from PIL import Image
import networkx as nx

def load_protein(pdb_path):
    """Load the protein structure from a PDB file."""
    return mda.Universe(pdb_path)

def load_ligand(sdf_path):
    """Load the ligand structure from an SDF file."""
    mol = Chem.SDMolSupplier(sdf_path, removeHs=False)[0]
    return mol

def calculate_interactions(protein, ligand, distance_cutoff=4.0):
    """Calculate interactions between protein and ligand atoms within a specified cutoff distance."""
    protein_atoms = protein.select_atoms('all')
    ligand_conf = ligand.GetConformer()
    interactions = []
    
    for atom in protein_atoms:
        protein_pos = atom.position
        for i in range(ligand.GetNumAtoms()):
            ligand_pos = ligand_conf.GetAtomPosition(i)
            distance = np.linalg.norm(protein_pos - np.array([ligand_pos.x, ligand_pos.y, ligand_pos.z]))
            
            if distance < distance_cutoff:
                interactions.append((atom, i, distance))
    
    return interactions

def get_one_letter_code(resname):
    """Convert three-letter amino acid codes to one-letter codes."""
    amino_acid_codes = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    return amino_acid_codes.get(resname, resname)

def plot_interactions_with_networkx(protein, ligand, interactions):
    """Plot interactions between a protein and ligand using NetworkX for visualization."""
    AllChem.Compute2DCoords(ligand)
    drawer = Draw.MolDraw2DCairo(500, 500)
    drawer.DrawMolecule(ligand)
    drawer.FinishDrawing()

    ligand_img_data = drawer.GetDrawingText()
    ligand_img = Image.open(BytesIO(ligand_img_data))
    
    fig, ax = plt.subplots()
    ax.imshow(ligand_img, extent=(0, 500, 0, 500))

    ligand_coords = []
    atom_colors = []
    for i in range(ligand.GetNumAtoms()):
        pos = drawer.GetDrawCoords(i)
        ligand_coords.append((pos.x, pos.y))
        
        atom = ligand.GetAtomWithIdx(i)
        atom_type = atom.GetSymbol()
        if atom_type == 'O':
            atom_colors.append('red')
        elif atom_type == 'N':
            atom_colors.append('blue')
        elif atom_type == 'C':
            atom_colors.append('black')
        else:
            atom_colors.append('gray')

    ligand_coords = np.array(ligand_coords)
    ligand_coords[:, 1] = 500 - ligand_coords[:, 1]

    unique_interactions = {}
    for atom, i, _ in interactions:
        residue_key = (atom.resname, atom.resid)
        if residue_key not in unique_interactions:
            unique_interactions[residue_key] = []
        unique_interactions[residue_key].append(i)

    G = nx.Graph()

    # Add nodes for ligand atoms
    for idx, (x, y) in enumerate(ligand_coords):
        G.add_node(f"Ligand_{idx}", pos=(x, y), type='ligand', color=atom_colors[idx])

    # Add nodes for unique protein residues
    for residue_key in unique_interactions:
        one_letter_code = get_one_letter_code(residue_key[0])
        G.add_node(residue_key, pos=None, type='protein', label=f'{one_letter_code}{residue_key[1]}')

    # Add edges between ligand atoms and residues
    for (residue_key, ligand_indices) in unique_interactions.items():
        for i in ligand_indices:
            G.add_edge(residue_key, f"Ligand_{i}")

    # Use a layout algorithm to position nodes
    pos = nx.spring_layout(G, 
                           pos={n: d['pos'] for n, d in G.nodes(data=True) if d['pos'] is not None}, 
                           fixed=[f"Ligand_{i}" for i in range(ligand.GetNumAtoms())],
                           k=40, iterations=100
                          )

    # Draw nodes and edges
    nx.draw(
        G, pos, ax=ax,
        with_labels=False,
        node_size=100,
        node_color=[d['color'] if d['type'] == 'ligand' else 'lightblue' for n, d in G.nodes(data=True)],
        font_weight='bold'
    )

    for node, (x, y) in pos.items():
        if G.nodes[node]['type'] == 'protein':
            label = G.nodes[node]['label']
            ax.text(x, y, label, ha='center', va='center', fontweight='bold', fontsize=10, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

    ax.axis('off')
    plt.tight_layout()
    plt.savefig('simpli.png')

def main():
    parser = argparse.ArgumentParser(description='Visualize protein-ligand interactions with NetworkX.')
    parser.add_argument('-p', '--protein', required=True, help='Path to the protein PDB file.')
    parser.add_argument('-l', '--ligand', required=True, help='Path to the ligand SDF file.')
    parser.add_argument('-d', '--distance', type=float, default=4.0, help='Distance cutoff for interactions (default: 4.0 Å).')

    args = parser.parse_args()

    protein = load_protein(args.protein)
    ligand = load_ligand(args.ligand)
    interactions = calculate_interactions(protein, ligand, args.distance)
    plot_interactions_with_networkx(protein, ligand, interactions)

if __name__ == "__main__":
    main()