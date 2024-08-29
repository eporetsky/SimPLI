import sys
import MDAnalysis as mda
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import matplotlib.pyplot as plt
import numpy as np
from io import BytesIO
from PIL import Image
from scipy.optimize import minimize

def load_protein(pdb_path):
    return mda.Universe(pdb_path)

def load_ligand(sdf_path):
    mol = Chem.SDMolSupplier(sdf_path, removeHs=False)[0]
    return mol

def calculate_interactions(protein, ligand, distance_cutoff=4.0):
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

def optimize_amino_acid_positions(unique_interactions, ligand_coords, num_trials=10):
    best_result = None
    best_total = float('inf')

    def total_distance(angles):
        total = 0
        positions = []
        for angle in angles:
            pos = np.array([250 + 200 * np.cos(angle), 250 + 200 * np.sin(angle)])
            positions.append(pos)
        
        # Check for overlap and adjust positions
        for i, pos1 in enumerate(positions):
            for j, pos2 in enumerate(positions):
                if i != j and np.linalg.norm(pos1 - pos2) < 80:  # 80 is the width of the rectangle
                    total += 500  # Penalize overlap heavily

        for (residue_key, ligand_indices), pos in zip(unique_interactions.items(), positions):
            for i in ligand_indices:
                ligand_pos = ligand_coords[i]
                total += np.linalg.norm(pos - ligand_pos)
                # Penalize if too close to ligand atoms
                if np.linalg.norm(pos - ligand_pos) < 80:  # Arbitrary threshold for "too close"
                    total += 500  # Penalize being too close to ligand atoms

        return total

    for _ in range(num_trials):
        initial_angles = np.random.uniform(0, 2*np.pi, len(unique_interactions))
        bounds = [(0, 2*np.pi) for _ in range(len(unique_interactions))]
        result = minimize(total_distance, initial_angles, bounds=bounds, method='L-BFGS-B')
        
        if result.fun < best_total:
            best_total = result.fun
            best_result = result.x

    return best_result

def plot_interactions(protein, ligand, interactions):
    # Generate a 2D depiction of the ligand
    AllChem.Compute2DCoords(ligand)
    drawer = Draw.MolDraw2DCairo(500, 500)
    drawer.DrawMolecule(ligand)
    drawer.FinishDrawing()

    # Convert the drawing to an image
    ligand_img_data = drawer.GetDrawingText()
    ligand_img = Image.open(BytesIO(ligand_img_data))

    # Set up a figure with matching pixel dimensions
    fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
    ax.imshow(ligand_img, extent=(0, 500, 0, 500))

    # Use RDKit to get scaled 2D coordinates of ligand atoms
    ligand_coords = []
    for i in range(ligand.GetNumAtoms()):
        pos = drawer.GetDrawCoords(i)
        ligand_coords.append((pos.x, pos.y))
    ligand_coords = np.array(ligand_coords)

    # Flip the y-coordinates around the x-axis through the center of the image
    ligand_coords[:, 1] = 500 - ligand_coords[:, 1]

    unique_interactions = {}
    for atom, i, _ in interactions:
        residue_key = (atom.resname, atom.resid)
        if residue_key not in unique_interactions:
            unique_interactions[residue_key] = []
        unique_interactions[residue_key].append(i)

    optimized_angles = optimize_amino_acid_positions(unique_interactions, ligand_coords)
    colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_interactions)))

    # Draw interaction lines first
    for idx, ((resname, resid), ligand_indices) in enumerate(unique_interactions.items()):
        angle = optimized_angles[idx]
        protein_pos = (250 + 200 * np.cos(angle), 250 + 200 * np.sin(angle))
        
        for i in ligand_indices:
            ligand_pos = ligand_coords[i]
            ax.plot([protein_pos[0], ligand_pos[0]], [protein_pos[1], ligand_pos[1]], 
                    color=colors[idx], linestyle='--', linewidth=2, alpha=0.7)

    # Draw rectangles and labels
    for idx, ((resname, resid), ligand_indices) in enumerate(unique_interactions.items()):
        angle = optimized_angles[idx]
        protein_pos = (250 + 200 * np.cos(angle), 250 + 200 * np.sin(angle))
        
        # Draw a rectangle to represent the amino acid label
        rect = plt.Rectangle((protein_pos[0]-40, protein_pos[1]-15), 80, 30, 
                             facecolor=colors[idx], edgecolor='black', alpha=0.7)
        ax.add_patch(rect)
        
        # Add the amino acid label
        ax.text(protein_pos[0], protein_pos[1], f"{resname}{resid}", 
                ha='center', va='center', fontweight='bold', fontsize=10)

    ax.set_xlim(0, 500)
    ax.set_ylim(0, 500)
    ax.set_aspect('equal', 'box')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig('simpli.png')

def main(protein_pdb_path, ligand_sdf_path, distance_cutoff):
    protein = load_protein(protein_pdb_path)
    ligand = load_ligand(ligand_sdf_path)
    interactions = calculate_interactions(protein, ligand, distance_cutoff)
    plot_interactions(protein, ligand, interactions)

# Replace 'protein.pdb' and 'ligand.sdf' with your file paths
protein_pdb_path = sys.argv[1]
ligand_sdf_path = sys.argv[2]
distance_cutoff = float(sys.argv[3])

main(protein_pdb_path, ligand_sdf_path, distance_cutoff)
