# Copyright (c) 2022 rxn4chemistry - Mark Martori Lopez

# DEPENDENCIES
import os
import detectron2
from detectron2.utils.logger import setup_logger
setup_logger()

# import some common detectron2 utilities
from detectron2 import model_zoo
#from detectron2.engine import DefaultPredictor
from detectron2.config import get_cfg
#from detectron2.utils.visualizer import Visualizer
from detectron2.data import MetadataCatalog
from detectron2.data.catalog import DatasetCatalog
# -------

## Prepare Custom Dataset
from detectron2.data.datasets import register_coco_instances
register_coco_instances("my_dataset_train", {}, "../../images/annotations/custom_train.json", "../../images/train2017")
register_coco_instances("my_dataset_val", {}, "../../images/annotations/custom_val.json", "../../images/val2017")
# register_coco_instances("my_dataset_test", {}, "images/json_annotest.json", "images/test")

# Prepare metrics
metricsPath = "output_silva_custom/metrics_json.json"
metrics2Path = "output_silva_custom/metrics.json"
if os.path.exists(metricsPath):
        os.remove(metricsPath)
if os.path.exists(metrics2Path):
        os.remove(metrics2Path)
# Create Output Directory
outputPath = "output_silva/"
os.makedirs(outputPath, exist_ok=True)

## Training Config
cfg = get_cfg()
cfg.merge_from_file(model_zoo.get_config_file("COCO-Detection/faster_rcnn_R_50_FPN_3x.yaml"))
cfg.DATASETS.TRAIN = ("my_dataset_train",)
cfg.DATASETS.TEST = ("my_dataset_val",)

cfg.DATALOADER.NUM_WORKERS = 4
cfg.OUTPUT_DIR = outputPath
cfg.MODEL.WEIGHTS = model_zoo.get_checkpoint_url("COCO-Detection/faster_rcnn_R_50_FPN_3x.yaml") # Let training initialize from model zoo
cfg.SOLVER.IMS_PER_BATCH = 4
cfg.SOLVER.BASE_LR = 1e-4 
cfg.SOLVER.WARMUP_ITERS = 1000
cfg.SOLVER.MAX_ITER = 45000 # 24100 = 51k images. 
cfg.SOLVER.STEPS = (15000, 22500)
cfg.SOLVER.GAMMA = 0.05
cfg.MODEL.ROI_HEADS.BATCH_SIZE_PER_IMAGE = 64
cfg.MODEL.ROI_HEADS.NUM_CLASSES = 4
cfg.TEST.EVAL_PERIOD = 4000 # Must be atleast 4000 8000


## Create Trainer Model from Facebook's Detectron2
from detectron2.engine import DefaultTrainer
from detectron2.evaluation import COCOEvaluator

class CocoTrainer(DefaultTrainer): # COCO Dataset trained with model in CFG.
    @classmethod
    def build_evaluator(cls, cfg, dataset_name, output_folder=None):
        if output_folder is None:
                os.makedirs("coco_eval", exist_ok=True)
                output_folder = "coco_eval"
        return COCOEvaluator(dataset_name, cfg, False, output_folder)

## --- TRAIN
import os
os.makedirs(cfg.OUTPUT_DIR, exist_ok=True)
trainer = CocoTrainer(cfg)
trainer.resume_or_load(resume = False)
trainer.train()
