
import subprocess
import numpy as np
from PIL import Image
from rdkit import Chem, RDLogger


def run_mol2smiles_dec(img_path):
    """
    Run DECIMER SMILES translation in a subprocess.
    """
    command = f"conda run --name DECIMER python functions/mol2smiles_dec.py"
    
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        return None
    else:
        print("Error:", result.stderr.strip())

# EXAMPLE USAGE

image_path = "mol_images/test2911.png"
run_mol2smiles_dec(image_path)