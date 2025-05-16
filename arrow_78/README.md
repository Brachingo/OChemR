# OChemR project

## Generate Training data set
### Target Structure:
- images/
    - custom_train.json     # JSONs for DETR training.
    - custom_val.json       #
    - json_annotest.json    #
    - train2017/
    - val2017/
    - test2017/
    - labelled_train2017/
    
##### Step 1: Set proper paths:
 - folder: arrow_78
 - script : params.json
        - train_path = "location where train images will be" ex: ../images/train2017/
        - labelled_path = "location where train images will be" ex: ../images/labelled2017/



##### Step 2: Set proper params: - generate 60k images. Check: info_params.py for more inf.
 - folder: arrowDS/params/
 - script: params.json
        -   {
            "dataset_params": { 
                "train_path" : "<train_path>",
                    "labelled_path" : "<labelled_path>",
                    "img_width" : [512,650,800,1024, 800, 1024],
                "img_height" : [512,650,800,1024, 800, 1024],
                "molecules_sizes" : [3,6,8,10, 8, 10],
                "molecules_rotations" : [0,30,100,330],
                "num_molecules_per_reaction" : 12,
                "num_reactions_per_epoch" : 10000,
                "epochs" : 6}
            }
 - By running 10000 reactions per epoch, we will get 60k images in 6 epochs.
 - For validation set, we will use 15% of the training set (1500 reactions per epoch). 9k images.
 - For test set, we only do 200 reactions per epoch. 1200 images.

##### Step 3: Run.
 - folder: arrow_78
        - terminal: bash < arrow.lsf # will run main.py with parameters.

##### Step 4: Check progress.
 - folder: arrow_78/logs/currentDatefolder/
        - code experiment_log.txt

##### Step 5: Loop for validation images.
 - Remember to change "train_path" in all scripts to validation path. ex: images/val/
 - Remember to also change number of reactions per epoch to 1500, and in file `labellingImage.py` the number of reactions that will be labelled (line 56)

##### Step 6: Create .json files to run DETR.
 - folder: DatasetsConversion/_Yolo2COCO/
 - script: main.py  # set proper categories
        - classes = [
                "molec",
                "arrow",
                "text",
                "plus"
            ]

 - TRAIN:
 - script: run.lsf
        - --path location where train images are. ex: ../../images/train2017/
        - --output json_annotation.json
 - terminal: bash < run.lsf

 - VAL:
 - script: runval.lsf
        - --path location where validation images are. ex: ../../images/val2017/
        - --output json_annotval.json
 - terminal: bash < runval.lsf

 - TEST:
 - Add images to test folder if not done yet.
 - script: runtest.lsf
        - --path location where validation images are. ex: ../../images/test2017/
        - --output json_annotest.json
 - terminal: bash < runtest.lsf


##### Step 7: Move json file to images/ folder:
 - folder: DatasetConversion/_Yolo2COCO/output/
        - terminal: mv json* "path to images/annotations/ folder"

## DETR Custom Dataset Training:
 - ##### Step 0: Clone DETR: 
 - terminal: git clone https://github.com/woctezuma/finetune-detr.git

 In order to get a better training than the one in the original repo, we will modify some files.

 - To train the model on previous weights, we have to follow the steps from the original repo. But instead of training the model in the Jupyter Notbeook, we train the model by running  `python train_d2.py` in the terminal. A personalized script meant to improve the training process is provided in the folder `detr/`.

 - To run inference on the model, we will use the script `inference_D2.py` in the folder `detr/detectron_2/detr/`. This script is meant to run inference on the model and save the output in a JSON file. We can run inferene on our test data, or on our real-world data. I personally recommend to run inference on `vizz_detectron2.ipynb`, as it gives a clear output on all the testing data we have, with also the posibility to sae the inference on real-world data in an folder as .png files.


 *STILL HAVE TO APPPLY CHANGES TO train_d2.py AND inference_D2.py*