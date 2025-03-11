import os
import sys
parts = os.path.abspath(__file__).split(os.sep)
# import importlib
import re
import glob

# Check for the decimer package
try:
    print("\nInitializing DECIMER. This can take about 30 seconds...")
    from DECIMER import predict_SMILES
except ImportError:
    print("""\nDECIMER not found. Please install the package before running this script.
        # Create a conda wonderland
        conda create --name DECIMER python=3.10.0 -y
        conda activate DECIMER

        # Equip yourself with DECIMER
        pip install decimer""")
    
    sys.exit(1)

# Select the folder containing the images
print("\nSelect the folder containing the images in PNG or JPG format.")

folder_path = ["../../images/val"]
if not folder_path:
    print("\nExiting...")
    sys.exit(1)

filelist = glob.glob(os.path.join(folder_path, '**/*.png'), recursive=True) + \
           glob.glob(os.path.join(folder_path, '**/*.jpg'), recursive=True) + \
           glob.glob(os.path.join(folder_path, '**/*.jpeg'), recursive=True)

out_dir = os.path.join(folder_path, 'out')
os.makedirs(out_dir, exist_ok=True)

base_outfile_name = 'smiles_out'
outfile_name = f"{base_outfile_name}.csv"
counter = 1
while os.path.exists(os.path.join(out_dir, outfile_name)):
    outfile_name = f"{base_outfile_name}_{counter}.csv"
    counter += 1

with open(os.path.join(out_dir, outfile_name), 'w') as f:
    f.write("Image,SMILES\n")
    
    for filename in filelist:
        # Unleash the power of DECIMER
        img_file = os.path.basename(filename)
        SMILES = predict_SMILES(filename)
        print(f"Image: {img_file} | Decoded SMILES: {SMILES}")

        # Write the SMILES to a csv file with the Image name as one column and SMILES as another column
    
        f.write(f"{img_file},{SMILES}\n")

if filelist:
    print("\nAll images processed. Check the 'out' directory for the SMILES csv file.")