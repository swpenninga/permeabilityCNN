# permeabilityCNN
## Code for generating 3D porous images and for training a CNN through Keras.

Workflow is usually as follows:
* `generate_data.m` lets you define size, porosity, etc. of the geometries, exports to `dataset.mat` in the same folder.

* `permeability.m` uses `dataset.mat` and calculates permeability, exports porosity `porosities.mat` and permeabilities `permeabilities.mat`

* You can now use your own matlab script (or later, a python script) to add symmetric permutations as datapoints.

* Open `dataprocessing.ipynb` to convert the matlab files to numpy arrays (this enlarges the files 5-fold because matlab's compression is lost). This file should be easy to configure to your system (ram size and disk space etc.)

* `keras.ipynb` imports the numpy notebooks in a loop-structure that should be changed to correspond to the datasets you created in the previous step. 

* to run `keras.ipynb` you have to make sure that the first cell shows `Num GPUs Available:  1` to recognize the GPU, the setup described there is most likely different for your PC and you should use the links down below.

https://www.tensorflow.org/install/source#gpu

https://www.tensorflow.org/install/gpu


The results of the report for the CNN are given in `predicted.csv` and `truth.csv`.
