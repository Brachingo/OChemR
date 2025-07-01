import os
from detectron2.engine import DefaultTrainer
from detectron2.config import get_cfg
from detectron2.data.datasets import register_coco_instances
from detectron2.data import MetadataCatalog, DatasetCatalog
from detectron2.model_zoo import get_checkpoint_url, get_config_file

# Step 1: Register your dataset
register_coco_instances(
    "train_dataset", {},
    "../images/json_annotation.json",
    "../images/train"
)

register_coco_instances(
    "val_dataset", {}, # Name and metadata
    "../images/json_annotval.json", # Path to json annotation file
    "../images/val" # Path to images
)

outputPath = "output_silva_custom"
os.makedirs(outputPath, exist_ok=True)

metricsPath = "output_silva_custom/metrics_json.json"
metrics2Path = "output_silva_custom/metrics.json"
if os.path.exists(metricsPath):
    os.remove(metricsPath)
if os.path.exists(metrics2Path):
    os.remove(metrics2Path)

# Step 2: Prepare the configuration
# Load Metadata
metadata = MetadataCatalog.get("train_dataset")

# Step 3: Configure the model
cfg = get_cfg()
cfg.merge_from_file(get_config_file("COCO-Detection/faster_rcnn_R_50_FPN_3x.yaml"))  # Use Faster R-CNN config
cfg.DATASETS.TRAIN = ("train_dataset",)
cfg.DATASETS.TEST = ("val_dataset",)
cfg.DATALOADER.NUM_WORKERS = 6  # Number of data loading threads
cfg.MODEL.WEIGHTS = get_checkpoint_url("COCO-Detection/faster_rcnn_R_50_FPN_3x.yaml")  # Pretrained weights
cfg.SOLVER.IMS_PER_BATCH = 4  # Images per batch
cfg.SOLVER.BASE_LR = 1e-4  # Learning rate
cfg.SOLVER.MAX_ITER = 23000  # Adjust based on dataset size
cfg.MODEL.ROI_HEADS.BATCH_SIZE_PER_IMAGE = 128  # Adjust based on GPU memory
cfg.MODEL.ROI_HEADS.NUM_CLASSES = 4  # Number of classes (1 for "scratch")
cfg.OUTPUT_DIR = outputPath  # Output directory

# Step 4: Create the trainer
trainer = DefaultTrainer(cfg)
trainer.resume_or_load(resume=False)
#trainer.train()

from detectron2.evaluation import COCOEvaluator, inference_on_dataset
from detectron2.data import build_detection_test_loader

# Load the trained model
cfg.MODEL.WEIGHTS = os.path.join(cfg.OUTPUT_DIR, "model_final.pth")  # Path to final weights

cfg.MODEL.ROI_HEADS.SCORE_THRESH_TEST = 0.5  # Set a threshold for predictions

# Initialize the evaluator
evaluator = COCOEvaluator("scratch_val", cfg, False, output_dir=cfg.OUTPUT_DIR)
val_dataset
# Build the data loader for validation set
val_loader = build_detection_test_loader(cfg, "val_dataset")

# Perform evaluation
results = inference_on_dataset(trainer.model, val_loader, evaluator)
print(results)
