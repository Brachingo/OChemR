import sys
import cv2

print("\n Importing DECIMER AI might take some seconds... \n")
from DECIMER import predict_SMILES # Import SMILES prediction function from DECIMER.ai
print("\n\nReading image...")

def img_to_smiles_decimer(image):
    try:

        cv2.imwrite("tmp.png", image) # Save the image
        print("Getting SMILES...")
        return predict_SMILES("tmp.png")  # Predict the SMILES string
    except Exception:
        return ""

# Restore stdout to print the final result only
sys.stdout = sys.__stdout__

# Load image and print only the result
#mol2 = "images/Englerin-A_3"
mol2 = "mol_images/Englerin-A_312"

img = cv2.imread("../"+mol2+".png")
if img is not None:
    print(f"\nThe SMILES code for {mol2} is:\n{img_to_smiles_decimer(img)}")
else:
    print("\nError: Image not found!")
