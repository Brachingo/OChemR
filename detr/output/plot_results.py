from models.detr import build_detr  # Adjust this import based on your DETR implementation
import torch

checkpoint = torch.load("output/checkpoint.pth", map_location="cpu", weights_only=False)
print(checkpoint.keys())  # Check available keys


# Create the model (ensure the num_classes matches your dataset)
model = build_detr(num_classes=91)  # Change 91 if using a custom dataset

# Load weights into the model
model.load_state_dict(checkpoint["model"])

# Set the model to evaluation mode
model.eval()
