# Copyright (c) 2022 rxn4chemistry - Mark Martori Lopez

import argparse
from email.mime import image
import os
import sys
import distutils.core
import torch
from pathlib import Path
import json
import os
import glob
import tensorflow as tf
import time
import matplotlib.pyplot as plt
import matplotlib
print("Primary imports done.")

p = os.path.abspath('../')
sys.path.append(p)


# ViT - Detect Objects in image (Molecules - Arrows/+ - Text)
from Backend.functions.detr_inference import detr_inference
from detectron_2.detr.models import build_model
print("Imports done for DETR inference.")
# OCR - Get text from detected Text bboxes from ViT.

from paddleocr import PaddleOCR
ocr = PaddleOCR(use_angle_cls=True, lang="en")  # Initialize PaddleOCR

print("Imports done for OCR inference.")

# Molvec - Translate Molecules to SMILES
from Backend.functions.mol2smiles_dec import img_to_smiles_decimer
import cv2

# Read Direction of images
from Backend.functions.findArrows import get_arrow_direction
from Backend.functions.finDirection import findClosestType,groupMolecules,orderArrows
from Backend.functions.order_annot_reaction import Order_Reaction, Annotate_Reaction
import json
print("Imports done for direction inference.")
# Output
from Backend.functions.storeresults import storeResults

#p = os.path.abspath('detr/detectron2/detr')
#sys.path.append(p)
# Some basic setup:
# Setup detectron2 logger
import detectron2
from detectron2.utils.logger import setup_logger
setup_logger()
print("Logger setup done.")
# import some common libraries

# import some common detectron2 utilities
from detectron2 import model_zoo
from detectron2.engine import DefaultPredictor
from detectron2.config import get_cfg
from detectron2.utils.visualizer import Visualizer
from detectron2.data import MetadataCatalog, DatasetCatalog
from detectron2.data.datasets import register_coco_instances
print("Imports done for detectron2.")

from d2.detr import add_detr_config
print("Imports done for detectron2 training.")

# current path
print(f"Current path for COCO annotation = {os.getcwd()}")
register_coco_instances("custom_train",
                        {},
                        "../images/annotations/custom_train.json",
                        "../images/train2017/")

register_coco_instances("custom_val",
                        {},
                        "../images/annotations/custom_val.json",
                        "../images/val2017/")



#### ARGUMENTS
def get_args_parser():
    parser = argparse.ArgumentParser("Arrow project - Input images to extract SMILES.", add_help = False)

    #### --------  Input
    parser.add_argument("--data_path", default = "images/", type = str)

    #### --------   Vision Transformer:
    parser.add_argument('--lr', default=1e-4, type=float)
    parser.add_argument('--lr_backbone', default=1e-5, type=float)
    parser.add_argument('--batch_size', default=6, type=int)
    parser.add_argument('--weight_decay', default=1e-4, type=float)
    parser.add_argument('--epochs', default=300, type=int)
    parser.add_argument('--lr_drop', default=200, type=int)
    parser.add_argument('--clip_max_norm', default=0.1, type=float,
                        help='gradient clipping max norm')

    # Model parameters
    parser.add_argument('--frozen_weights', type=str, default=None,
                        help="Path to the pretrained model. If set, only the mask head will be trained")
    # * Backbone
    parser.add_argument('--backbone', default='resnet50', type=str,
                        help="Name of the convolutional backbone to use")
    parser.add_argument('--dilation', action='store_true',
                        help="If true, we replace stride with dilation in the last convolutional block (DC5)")
    parser.add_argument('--position_embedding', default='sine', type=str, choices=('sine', 'learned'),
                        help="Type of positional embedding to use on top of the image features")

    # * Transformer
    parser.add_argument('--enc_layers', default=6, type=int,
                        help="Number of encoding layers in the transformer")
    parser.add_argument('--dec_layers', default=6, type=int,
                        help="Number of decoding layers in the transformer")
    parser.add_argument('--dim_feedforward', default=2048, type=int,
                        help="Intermediate size of the feedforward layers in the transformer blocks")
    parser.add_argument('--hidden_dim', default=256, type=int,
                        help="Size of the embeddings (dimension of the transformer)")
    parser.add_argument('--dropout', default=0.1, type=float,
                        help="Dropout applied in the transformer")
    parser.add_argument('--nheads', default=8, type=int,
                        help="Number of attention heads inside the transformer's attentions")
    parser.add_argument('--num_queries', default=100, type=int,
                        help="Number of query slots")
    parser.add_argument('--num_classes', type=int, default=5, 
                        help='Number of object classes (including no-object)')
    parser.add_argument('--pre_norm', action='store_true')

    # * Segmentation
    parser.add_argument('--masks', action='store_true',
                        help="Train segmentation head if the flag is provided")

    # # Loss
    parser.add_argument('--no_aux_loss', dest='aux_loss', action='store_false',
                        help="Disables auxiliary decoding losses (loss at each layer)")
    # * Matcher
    parser.add_argument('--set_cost_class', default=1, type=float,
                        help="Class coefficient in the matching cost")
    parser.add_argument('--set_cost_bbox', default=5, type=float,
                        help="L1 box coefficient in the matching cost")
    parser.add_argument('--set_cost_giou', default=2, type=float,
                        help="giou box coefficient in the matching cost")
    # * Loss coefficients
    parser.add_argument('--mask_loss_coef', default=1, type=float)
    parser.add_argument('--dice_loss_coef', default=1, type=float)
    parser.add_argument('--bbox_loss_coef', default=5, type=float)
    parser.add_argument('--giou_loss_coef', default=2, type=float)
    parser.add_argument('--eos_coef', default=0.1, type=float,
                        help="Relative classification weight of the no-object class")

    # dataset parameters
    parser.add_argument('--dataset_file', default='arrow')
    parser.add_argument('--data_panoptic_path', type=str)
    parser.add_argument('--remove_difficult', action='store_true')

    parser.add_argument('--device_detr', default='cpu',
                        help='device to use for training / testing')
    parser.add_argument('--resume', default = 'detr/output/checkpoint.pth', help='Get model weights from checkpoint')

    parser.add_argument('--thresh', default=0.95, type=float)
    #### ----- End ViT. ---------------

    # Device
    parser.add_argument("--device", default = "cpu",
                        help="device used for loading model and inference.")

    # Debugging
    parser.add_argument("--debugging", default = False, help = "Get images with bboxes croped, arrows start/end points.")

    # Output
    parser.add_argument("--output_dir", default ="output/", help= " Path where to store SMILES and/or Images detection results.")
    parser.add_argument("--detection", default = False,
                        help = "Store images with bounding boxes detected in output folder.")
    parser.add_argument('--detection_dir',
                        help='path where to save the results, empty for no saving')
    return parser

### MAIN
def main(args_detr):
    """
    1 - Feed ViT:
        Output : dict
                detr_outputs : {img_name : [labels, bboxes] , next_img_name : [labels, bboxes]}
    
    2 - Find direction in images.

    3 - OCR text, OMR molecules, OSR arrows.

    4 - Join results -> SMILES.
    """

    device_detr = torch.device(args_detr.device_detr)
    data_path = args_detr.data_path
    debugging = args_detr.debugging                                # Debugging ? 

    ### 0: Read images - check extension. --------------
    img_files = []
    types = ('*.png', '*.jpg','*.jpeg') 
    for files in types:
        img_files.extend(glob.glob(data_path + files))

    if os.path.exists(args_detr.output_dir+"smiles.txt"):
        os.system("rm "+args_detr.output_dir+"smiles.txt")

    print(f"Analyzing files = {img_files}.\n")
    
    print(f"Building DETR and OCR models ...")
    ### 1: Feed images to Vision Transformer and get outputs. -----------------
    if args_detr.detection_dir: # Save images with bboxes and attention.
        detection_path = args_detr.output_dir+args_detr.detection_dir
        Path(detection_path).mkdir(parents=True, exist_ok=True)
    else:                                                          # Only get output of detr.
        detection_path = "none"

    threshold = args_detr.thresh
    model, _, _ = build_model(args_detr) # Get DETR model with current args.
    
    # Print my current path
    print(f"Current path for model loading = {os.getcwd()}")
    
    
    # Load model weights from model_final.pth
    cfg = get_cfg()
    add_detr_config(cfg)
    cfg.merge_from_file("../detectron_2/detr/d2/configs/detr_256_6_6_torchvision.yaml")
    cfg.MODEL.ROI_HEADS.SCORE_THRESH_TEST = 0.5 # Set threshold for this model
    cfg.DATASETS.TRAIN = ("custom_train",)
    cfg.DATASETS.TEST = ("custom_val",)
    cfg.MODEL.WEIGHTS = '../detectron_2/output/model_final.pth' # Set path model .pth
    cfg.MODEL.DETR.NUM_CLASSES = 4
    cfg.MODEL.ROI_HEADS.NUM_CLASSES = 4
    cfg.MODEL.ROI_HEADS.SCORE_THRESH_TEST = 0.95 # Confidence threshold for this model
    model = DefaultPredictor(cfg)
        
    
    ## 2: Load OCR for text_imgs  - EasyOCR --------------------
    #reader = easyocr.Reader(['en']) #, gpu= gpu # this needs to run only once to load the model into memory 
    # reader was used then using EasyOCR to get text from images. Now we use PaddleOCR.

    print("DETR Inference ...")
    # Run detr inference
    detr_outputs = detr_inference(img_files, model, device_detr, detection_path, args_detr.thresh)

    ## 3: Translate all bboxes found by ViT:
    for it, (image_path, outputs) in enumerate(list(detr_outputs.items())):
        timest = time.time()                                   # Time it!
        img = cv2.imread(image_path)                           # Read image
        print(f"Image {image_path} generating SMILES ...")
        height,width,_ = img.shape
        filename = os.path.basename(image_path)
        labels = outputs[0]                                    # Get Labels and bboxes...
        boxes = outputs[1]                                     # ...detected of current image.

        boxes_dict = {}                                        # {0 : bbox, 1 : bbox ...}
        SMILES_dict = {}                                       # {0 : SMILES, 1 : text, 2 : start/end arrow ...}
        
        for step,(x,y,w,h) in enumerate(boxes):                # Iterate by bbox
            boxes_dict[step] = [int(x),int(y),int(w),int(h)]   # Store coordinates of the bbox. Key = index.
            label = labels[step]                               # Store label of that bbox.
            x,y,w,h = int(x),int(y),int(w),int(h)              # Convert bbox coords. to int.
            padd = 5                                           # If possible, crop bbox with a 5+padding.
            if x- padd >= 0 and y - padd >= 0 and w + padd <= width and h + padd <= height:
                new_img = img[y - padd:h + padd, x - padd :w + padd]
            else:                                              # If not possible, crop the proper bbox.
                new_img = img[y:h, x:w]  
                                                               # TRANSLATION starts:
            if label == 0: # Get SMILES                        # If it is a molecule...
                smiles = img_to_smiles_decimer(new_img)        # Run DECIMER AI to get SMILES
                if smiles == "":                               # If it can't translate
                    SMILES_dict[step] = "molecule"+str(step)   # If Molvec can't translate, we store molecule+keyidx.
                else: 
                    SMILES_dict[step] = smiles                 # Store SMILES
                                                  # DEBUGGING
                molname = "mol_images/"+filename.replace(".",str(step)+".")
                cv2.imwrite(molname, new_img)


            elif label == 1:                                   # If bbox is an arrow:
                startpoint,endpoint = get_arrow_direction(new_img,image_path,step, debugging)
                SMILES_dict[step] = [startpoint,endpoint]      # Store Start/End points of arrow.
                if debugging:
                    print(f"Arrow startpoint = {startpoint} and endpoint = {endpoint}.\n")

                arrowname = "arrows/"+filename.replace(".",str(step)+".")
                cv2.imwrite(arrowname,new_img)
                
            elif label == 2:                                   # Text - OCR - PADDLEOCR
                textname = "text_images/"+filename.replace(".",str(step)+".")
                cv2.imwrite(textname,new_img)
                # PaddleOCR
                # PaddleOCR
                result = ocr.ocr(textname, cls=True)                
                # Extract text from PaddleOCR result
                extracted_text = []
                for line in result:
                    #Skip this part when you get error TypeError: 'NoneType' object is not iterable
                    if line is None:
                        continue
                    for word_info in line:                        
                        extracted_text.append(word_info[1][0])  # word_info[1][0] contains detected text
                        
                try:
                    SMILES_dict[step] = extracted_text         # Store text.
                except:
                    SMILES_dict[step] = " "

            else: # + symbol.
                SMILES_dict[step] = "."

        if debugging:
            print(f" Image = {filename} : \n")
            print(f"{SMILES_dict}")
            print(f"Labels:\n {labels}")
            print(f"bboxes:\n{boxes_dict}")


        ## 4: --------------- Find overall reaction direction and output full path SMILES: --------------------
        #
        # HOW: - By using arrows start/end points, we can calculate mindistance to all molecules. Then create path.
        #
        # Problem:
        #      1 - If start/end points of arrow are not the correct ones, the reaction will not make sense at the end.
        # Solutions:
        #      1 - Improve sp,ep arrow retrieval and image retrieval (cleaner images).
        #       or, we can store only SMILES that take part in reaction. (w/o direction)
        
        # Group molecules that are summing up and MODIFY dictionaries by deleting the key of the first molecule of the sum, and joining its SMILES to the second molecule:
        boxes_dict,SMILES_dict,labels = groupMolecules(boxes_dict,SMILES_dict, labels)

        unorderedReaction = {}                                 # Unordered Final Dictionary to dump on Json.
        if debugging: print(f"Bboxes_dict after groupmols = {boxes_dict.values()}\n Labels = {labels}")
        for key in boxes_dict:                                 # Iterate in boxes_dict keys, in case we deleted a molecule index when joining molecules that are summing.
            if labels[key] == 1:                               # and check label from labels list.
                sp,ep = SMILES_dict[key]                       # Get startpoint, endpoint of arrow:
                pos_prev,middle,pos_post = findClosestType(sp,ep, boxes_dict, labels,key, typeobj = 0, threshold = 85)
                textinf = []                                   # Returns list of molecules indexes that are closer to arrowidx than threhsold.           
                if len(middle) > 0:                            # Ex: mol1_index, text1_index, mol3_index
                    for a in middle: textinf.append(a)         # Store text information of the arrow.
                prevmol = SMILES_dict[pos_prev]                # Get SMILES molecule previous to the arrow.
                postmol = SMILES_dict[pos_post]                # Get SMILES molecule to which the arrow points.
                                                               # Find text (typeobj = 2) closest to arrow:
                pos_prev,middle,pos_post = findClosestType(sp,ep, boxes_dict, labels,key, typeobj = 2, threshold = 70) # returns list of text associated with arrow.
                if len(middle) > 0:
                    for a in middle: textinf.append(a)
                textinf = [SMILES_dict[txt] for txt in textinf]
                unorderedReaction["arrow"+str(key)] = {
                                                "prev_mol":prevmol, # Create Arrow structure in main dictionary.
                                                "text"    :textinf, # Add text of current process/arrow.
                                                "post_mol":postmol  # Molecule to which the arrow points.
                                                }
                # Save the outputs in a dictionary.
                """
                dict = {
                    "Arrow ID": key,
                    "Previous Molecule": prevmol,
                    "Text": textinf,
                    "Post Molecule": postmol
                }
                # Print the dictionary into a JSON File
                with open(args_detr.output_dir+filename.replace(".png", ".json"), "a") as f:
                    json.dump(dict, f, indent=4)
                """
                    
                if debugging:
                    print(f"For arrow = {key}:\n")
                    print(f"Prev mol = {prevmol}")
                    print(f"Middle text = {textinf}")
                    print(f"Post mol = {postmol}")
                    
        final_ordered_reaction = Order_Reaction(unorderedReaction)
        
        # Annotate reaction with the correct information.
        Annotate_Reaction(final_ordered_reaction, args_detr.output_dir, filename)
        
                
        #final_ordered_reaction, SMILESresult = orderArrows(unorderedReaction)    # Order unordered dictionary.
        #print(f"Unordered reaction = {unorderedReaction}")
        #print(f"Final ordered reaction = {final_ordered_reaction}")
        #if final_ordered_reaction:
            #storeResults(final_ordered_reaction, filename, args_detr.output_dir) # Store results as Json

        #    with open(args_detr.output_dir+"smiles.txt","a") as f:
        #        f.write(f"Filename = {filename}\nSMILES = {SMILESresult[0:-2]}\n\n")
        timefin = time.time()
        if debugging:
            print(f"Total time to translate, find direction and assemble = {abs(timefin - timest)}s.")


if __name__ == '__main__':
    parser_detr = argparse.ArgumentParser("ARROW Project - Retrieve SMILES from images.", parents=[get_args_parser()])
    args_detr = parser_detr.parse_args()
    if args_detr.output_dir:
            Path(args_detr.output_dir).mkdir(parents=True, exist_ok=True)
           
    main(args_detr)    

