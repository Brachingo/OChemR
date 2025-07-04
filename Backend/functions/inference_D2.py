# Copyright (c) 2022 rxn4chemistry - Mark Martori Lopez

import os
import cv2
import glob
from PIL import Image
import time
import random

# From Facebook's detectron2:
from d2.detr import add_detr_config
from detectron2 import model_zoo
from detectron2.engine import DefaultPredictor
from detectron2.config import get_cfg
from detectron2.utils.visualizer import Visualizer
from detectron2.data import MetadataCatalog
from detectron2.data.datasets import register_coco_instances
from detectron2.utils.visualizer import ColorMode
# -------
# Make real directory
if not os.path.exists('images/real/real'):
    os.makedirs('images/real/real')
#register_coco_instances("my_dataset_test", {}, "../images/json_annotest.json", "../images/test/")
# Try on new (Real World) images:
register_coco_instances("my_dataset_test", {}, "../../../real_data/json_annotreal.json", "../../../real_data/real")

## Prepare Config File:
cfg = get_cfg()
add_detr_config(cfg)
cfg.merge_from_file("d2/configs/detr_256_6_6_torchvision.yaml")
cfg.MODEL.WEIGHTS = "outputs_4/model_final.pth" # os.path.join(cfg.OUTPUT_DIR, "model_final.pth")
cfg.DATASETS.TEST = ("my_dataset_test", )
cfg.MODEL.DETR.NUM_CLASSES = 4
cfg.DATALOADER.NUM_WORKERS = 2
cfg.SOLVER.IMS_PER_BATCH = 2
cfg.SOLVER.BASE_LR = 0.00025  # pick a good LR
cfg.SOLVER.MAX_ITER = 30000    # 300 iterations seems good enough for this toy dataset; you may need to train longer for a practical dataset
cfg.MODEL.ROI_HEADS.BATCH_SIZE_PER_IMAGE = 512   # faster, and good enough for this toy dataset (default: 512)
cfg.MODEL.ROI_HEADS.NUM_CLASSES = 4  # only has one class (ballon)
predictor = DefaultPredictor(cfg)
test_metadata = MetadataCatalog.get("my_dataset_test")


## Define Bounding Box Creation Function + Labels:
def plot_one_box(x, img, color=None, label=None, line_thickness=0.5):
    tl = 2 # line_thickness or round(0.002 * (img.shape[0] + img.shape[1]) / 2) + 1  # line/font thickness
    # color = color or [random.randint(0, 255) for _ in range(3)]
    c1, c2 = (int(x[0]), int(x[1])), (int(x[2]), int(x[3]))
    cv2.rectangle(img, c1, c2, color, thickness=tl, lineType=cv2.LINE_AA)

    if label:
        tf = tl #max(tl - 1, 1)  # font thickness
        t_size = cv2.getTextSize(label, 0, fontScale=tl / 4, thickness=tf)[0]
        c2 = c1[0] + t_size[0], c1[1] - t_size[1] - 3
        cv2.rectangle(img, c1, c2, (45, 230, 230), -1, cv2.LINE_AA)  # filled
        cv2.putText(img, label, (c1[0], c1[1] - 2), 0, tl / 4, (0,0,0), thickness=tf, lineType=cv2.LINE_AA)

# Colors for visualization
COLORS = [[119, 179, 50], [255, 0, 0], [0, 0, 255],
          [0, 255, 0], [0, 114, 178], [186, 0, 128]]

CLASSES = [
    'N/A','molec','arrow','text','letter','cbox','plus' 
]


## ----------------- Run Test: Get inference time.
threshold = 0.95
tested_images = ['../../../real_data/real/*png'] # 'images/val/*png'
for paths in tested_images:
    for it,imageName in enumerate(glob.glob(paths)):
        filename = os.path.basename(imageName)
        im = cv2.imread(imageName)
        start_t = time.perf_counter()
        outputs = predictor(im)
        end_t = time.perf_counter()
        print(f"Inference time for img {filename} = {abs(end_t - start_t)} seconds.")
        v = Visualizer(im[:, :, ::-1],
                        metadata=test_metadata, 
                        scale=0.8
                        )
        outputs["instances"] = outputs["instances"].to("cpu")
        keep = outputs["instances"].scores > threshold
        boxes = outputs['instances'][keep].pred_boxes
        boxes = boxes.to("cpu")

        scores = outputs['instances'][keep].pred_classes
        scores = scores.to("cpu")
        scores = scores.tolist()

        # confidence
        conf = outputs['instances'][keep].scores
        conf = conf.to("cpu")
        conf = conf.tolist()
        
        for it,box in enumerate(boxes):
            box = box.tolist()
            label = scores[it]
            confidence = int(conf[it] * 100)
            text = f"{label} {confidence:.2f}"
            color = COLORS[label]
            plot_one_box(box, im, color, label=text)

        cv2.imwrite('images/real/real/'+filename, im)
