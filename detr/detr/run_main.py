!python main.py \
  --dataset_file "custom" \
  --coco_path "../../images/" \
  --output_dir "outputs" \
  --resume "detr-r50_no-class-head.pth" \
  --num_classes $num_classes \
  --epochs 30