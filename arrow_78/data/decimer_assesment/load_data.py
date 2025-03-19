from rdkit import Chem
from rdkit.Chem import Draw
import os

# Make a list with the first 1000 lines of the molecules.txt file
molecules = []
with open('molecules_1000_random.txt', 'r') as file:
    # Get all lines in the file
    for i in range(1000):
        molecules.append(file.readline().strip())
        
# Store images from each molecule into directory "molecule_images"
if not os.path.exists('molecule_images'):
    os.makedirs('molecule_images')
    
for i, molecule in enumerate(molecules):
    mol = Chem.MolFromSmiles(molecule)
    if mol is not None:
        Draw.MolToFile(mol, f'molecule_images/molecule_{i}.png')
    else:
        print(f"Could not convert molecule {i} to RDKit Mol object")
        
print(f"{len(molecules)} images have been saved in the directory molecule_images")

# Terminal comands: python load_data.py
# Output: 1000 images have been saved in the directory molecule_images