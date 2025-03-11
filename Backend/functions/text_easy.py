import easyocr
import cv2

def detect_text(image_path):
    reader = easyocr.Reader(['en'])
    result = reader.readtext(image_path)
    return result

if __name__ == "__main__":
    image_path = '../images/test29.png'
    detected_text = detect_text(image_path)
    for (bbox, text, prob) in detected_text:
        print(f"Detected text: {text} (Confidence: {prob:.2f})")