# Copyright (c) 2022 rxn4chemistry - Mark Martori Lopez

import os
import sys
import cv2
from PIL import Image
import random
import torch
import numpy as np
import matplotlib.pyplot as plt
p = os.path.abspath('../detr/detr/')
sys.path.append(p)
#from datasets.coco import make_coco_transforms # Useful before

## Rescaling/Filtering bboxes ----------
def box_cxcywh_to_xyxy(x):
    x_c, y_c, w, h = x.unbind(1)
    b = [(x_c - 0.5 * w), (y_c - 0.5 * h),
         (x_c + 0.5 * w), (y_c + 0.5 * h)]
    return torch.stack(b, dim=1)

def rescale_bboxes(out_bbox, size):
    img_w, img_h = size
    b = box_cxcywh_to_xyxy(out_bbox)
    b = b * torch.tensor([img_w, img_h,
                          img_w, img_h
                          ], dtype=torch.float32)
    return b

def filter_boxes(scores, boxes, confidence=0.9, apply_nms=True, iou=0.5):
    keep = scores.max(-1).values > confidence                           # If one category is scored > threshold
    scores, boxes = scores[keep], boxes[keep]                           # keep that 'detection'
    return scores, boxes
##    ------------      -----------------------------


## PLOT function -------------
def plot_one_box(x, img, color=None, label=None, line_thickness=1):
    tl = 2                                                              # line_thickness or round(0.002 * (img.shape[0] + img.shape[1]) / 2) + 1  # line/font thickness
    color = color or [random.randint(0, 255) for _ in range(3)]
    c1, c2 = (int(x[0]), int(x[1])), (int(x[2]), int(x[3]))
    cv2.rectangle(img, c1, c2, color, thickness=tl, lineType=cv2.LINE_AA)

    if label:
        tf = tl                                                         # font thickness
        t_size = cv2.getTextSize(label, 0, fontScale=tl / 4, thickness=tf)[0]
        c2 = c1[0] + t_size[0], c1[1] - t_size[1] - 3
        cv2.rectangle(img, c1, c2, (45, 230, 230), -1, cv2.LINE_AA)     # filled
        cv2.putText(img, label, (c1[0], c1[1] - 2), 0, tl / 4, [0, 0, 0], thickness=tf, lineType=cv2.LINE_AA)
## -------------       -------------------------


CLASSES = [
    'N/A','molec','arrow','text','letter','cbox','plus' 
]

# colors for visualization
COLORS = [[255,255,255],[119, 179, 50], [204, 83, 23], [236, 177, 32],
          [126, 47, 142], [0, 114, 178], [77, 190, 238]]

CLASSES = ['molecule', 'arrow', 'text', 'plus']  # Update based on your dataset


# Replace the detr_inference function with:
@torch.no_grad()
def detr_inference(images_path, predictor, device, output_path, threshold):
    outputs_dict = {}
    for img_sample in images_path:
        filename = os.path.basename(img_sample)
        im = cv2.imread(img_sample)
        outputs = predictor(im)  # Use detectron2's predictor
        instances = outputs["instances"]
        
        # Filter by confidence
        keep = instances.scores > threshold
        scores = instances.scores[keep].cpu().numpy()
        boxes = instances.pred_boxes[keep].tensor.cpu().numpy()
        labels = instances.pred_classes[keep].cpu().numpy()
        
        # Convert boxes to [xmin, ymin, xmax, ymax] format
        boxes = [[xmin, ymin, xmax, ymax] for xmin, ymin, xmax, ymax in boxes]
        
        # Store results
        outputs_dict[img_sample] = [labels.tolist(), boxes]
        
        # Visualization (optional)
        """v = Visualizer(im[:, :, ::-1], MetadataCatalog.get(cfg.DATASETS.TRAIN[0]), scale=1.2)
        out = v.draw_instance_predictions(instances[keep].to("cpu"))
        cv2.imwrite(os.path.join("inference", filename), out.get_image()[:, :, ::-1])"""
            
        # Save the inference image
        #cv2.imwrite(os.path.join("inference", filename), im)    
            
    return outputs_dict