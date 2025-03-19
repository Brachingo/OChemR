import sys
import cv2
import os

print("\n Importing DECIMER AI might take some seconds... \n")
from DECIMER import predict_SMILES # Import SMILES prediction function from DECIMER.ai

# Store in each line of a new file the SMILES code of the molecules in the images, Store the smiles in a file

predicted_smiles = []

# Iterate over all the images in the directory molecule_images and predict the SMILES code, store the SMILES code in a list
for i, file in enumerate(os.listdir('molecule_images')):
    img = cv2.imread(f"molecule_images/{file}")
    if img is not None:
        smiles = predict_SMILES(img)
        # Put the SMILES code in the file directly
        predicted_smiles.append(smiles)
        if len(predicted_smiles) == 250:
            print(f"Processed 25% of the images")
        elif len(predicted_smiles) == 500:
            print(f"Processed 50% of the images")
        elif len(predicted_smiles) == 750:
            print(f"Processed 75% of the images")
    else:
        print(f"Could not read image {file}")

print("\nAll images have been processed\nNow writing the SMILES codes to a file...\n")
        
# Write the SMILES code in a new file
with open('predicted_smiles_random.txt', 'w') as file:
    for smiles in predicted_smiles:
        file.write(smiles + '\n')
        
# Put the smiles codes directly into the predicted¡_smiles.txt file
        
print("SMILES codes have been written to predicted_smiles.txt")

# Restore stdout to print the final result only
sys.stdout = sys.__stdout__

# Print the number of SMILES codes predicted
print(f"\n{len(predicted_smiles)} SMILES codes have been predicted")
