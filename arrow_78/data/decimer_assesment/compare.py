import sys

# Compare predicted_smiles.txt against molecules_1000.txt
# and output the number of correct predictions and % accuracy
'''
def compare_smiles(predicted_smiles, true_smiles):
    # Run each line of predicted_smiles against all lines of true_smiles
    # and count the number of correct predictions
    correct = 0
    
    with open(predicted_smiles, 'r') as f:
        predicted = f.readlines()
    with open(true_smiles, 'r') as f:
        true = f.readlines()

    for _ in predicted.lower():
        if _ in true.lower():
            correct += 1
    accuracy = correct / len(predicted)
    
    return correct, accuracy
'''
def compare_smiles(predicted_smiles, true_smiles):
    # Run each line of predicted_smiles against all lines of true_smiles
    # and count the number of correct predictions
    correct = 0
    
    with open(predicted_smiles, 'r') as f:
        predicted = [line.strip().lower() for line in f]  # Convert each line to lowercase and strip whitespace
    with open(true_smiles, 'r') as f:
        true = [line.strip().lower() for line in f]  # Convert each line to lowercase and strip whitespace

    for pred in predicted:
        if pred in true:
            correct += 1
    accuracy = correct / len(predicted)
    
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
    
# Terminal comands: python compare_smiles.py predicted_smiles.txt true_smiles.txt

if __name__ == '__main__':
    main()