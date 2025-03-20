import sys
from rdkit import Chem


def canonicalize_smiles(smiles_list):
    """Canonicalize a list of SMILES strings."""
    canonical_smiles = []
    invalids = 0 # Count invalid SMILES
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi.strip())
        if mol:
            canonical_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))
        else:
            canonical_smiles.append(None)  # Invalid SMILES
            invalids += 1
    print(f"Invalid SMILES: {invalids}")
    return canonical_smiles

def compare_smiles(predicted_smiles, true_smiles):
    """Compare predicted SMILES against true SMILES and calculate accuracy."""
    correct = 0
    
    with open(predicted_smiles, 'r') as f:
        predicted = f.readlines()
    with open(true_smiles, 'r') as f:
        true = f.readlines()
    
    # Canonicalize both predicted and true SMILES
    canonical_predicted = canonicalize_smiles(predicted)
    canonical_true = canonicalize_smiles(true)
    
    # Compare canonical SMILES
    print(f"The lengths of both lists of canonicalized SMILES are: {len(canonical_predicted)} & {len(canonical_true)}")
    print(f"First 5 canonicalized SMILES from predicted list: {canonical_predicted[:5]}")
    for pred_smi in canonical_predicted:
        if pred_smi in canonical_true:
            correct += 1

    accuracy = correct / len(canonical_predicted) if len(canonical_predicted) > 0 else 0
    
    return correct, accuracy


def main():
    if len(sys.argv) != 3:
        print("Usage: compare_smiles.py predicted_smiles.txt true_smiles.txt")
        sys.exit(1)

    predicted_smiles = sys.argv[1]
    true_smiles = sys.argv[2]

    correct, accuracy = compare_smiles(predicted_smiles, true_smiles)
    print(f"Correct predictions: {correct}")
    print(f"Accuracy: {accuracy:.2%}")
    
# Terminal commands: python compare_smiles.py predicted_smiles.txt true_smiles.txt

if __name__ == '__main__':
    main()