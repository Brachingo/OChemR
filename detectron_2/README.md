## DETR with D2 Wrapper
### Trainign:

    - Examine the Jupyter notebook `detr_finetune.ipynb` for loading all necessary steps and auxiliary functions DETR model training needs.
    - To train the model, I prefer using a separate file, you can create it yourself easily, which is then executed in the terminal, or you can execute it in a screen if you are working on external servers.
    - Run the trainng.
    - The model will be trained on the synthetic dataset created in the previous step, and evaluated on the validation set. The model will be saved in the output folder of choice.

### Evaluation
- Done in the training while its running. The evaluation is done every x iterations. The evaluation is done on the validation set. 

### Inference

    - The inference is done in the file `vizz_detectron2.ipynb`. This Jupyter notebook will load the trained model and run inference on the test set, or any choice of dataset.
    - It is an easy and flexible way to run inference on any dataset, as it allows you to change the dataset and the model easily. While also rapidly eaxmining the results.
    - Run the last cell to save the inference of the dataset of choice.