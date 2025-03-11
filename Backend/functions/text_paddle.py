from paddleocr import PaddleOCR
import cv2

# Only use for visualization (OCR_PADDLE WONDERLAND)
#from matplotlib import pyplot as plt 
#import numpy as np

# Initialize PaddleOCR
ocr = PaddleOCR(use_angle_cls=True, lang='en')  # Set use_gpu=False if GPU is not available

# Path to the image
img_path = '../images/test29.png'
img = cv2.imread(img_path)

# Run the OCR
result = ocr.ocr(img, cls=True)

# Print the results
for idx, line in enumerate(result):
    print(f"Line {idx + 1}:")
    for word in line:
        text = word[1][0]  # Extracted text
        confidence = word[1][1]  # Confidence score
        box = word[0]  # Bounding box coordinates
        print(f"  Text: {text}, Confidence: {confidence:.4f}, Bounding Box: {box}")

# Visualize the image with bounding boxes ||| ONLY USE ON OCR_PADDLE WONDERLAND
"""
img_rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)  # Convert to RGB for matplotlib
plt.figure(figsize=(10, 10))
plt.imshow(img_rgb)

# Draw bounding boxes on the image
for line in result:
    for word in line:
        box = word[0]  # Bounding box coordinates
        # Convert box coordinates to ints
        box = [(int(point[0]), int(point[1])) for point in box]
        # Draw the bounding box
        cv2.polylines(img_rgb, [np.array(box)], isClosed=True, color=(0, 255, 0), thickness=2)
        # Annotate the text
        cv2.putText(img_rgb, word[1][0], (box[0][0], box[0][1]), cv2.FONT_HERSHEY_SIMPLEX, 0.8, (255, 0, 0), 1)

# Show the image with bounding boxes and text
plt.imshow(img_rgb)
plt.axis('off')  # Hide axes
plt.show()
"""