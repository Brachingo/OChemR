
from doctr.models import ocr_predictor

model = ocr_predictor(pretrained=True)
print(ocr_predictor.__module__)
from doctr.io import DocumentFile
# Load a PDF document
doc = DocumentFile.from_images('../images/Englerin-A_3.png')
#doc = DocumentFile.from_images('../images/test29.png')
# Run the model
result = model(doc)
print(result)
result.show()
