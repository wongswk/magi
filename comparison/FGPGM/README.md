# Fast Gaussian Process Based Gradient Matching

Implementation of FGPGM in python 2.7, based on the paper "Fast Gaussian Process Based Gradient Matching for Parameter Identification in Systems of Nonlinear ODEs" by Philippe Wenk, Alkis Gotovos, Stefan Bauer, Nico Gorbach, Andreas Krause and Joachim M. Buhmann. If you use this code for further research, please consider citing our paper (http://arxiv.org/abs/1804.04378
).

## Out of the Box Experiments

The code must be run from a separate folder where the scripts can store .csv files and plots. Currently, there are two experiments implemented that should run out of the box. These experiments can be found in the folder "mainFiles". Each experiment can be run by executing the following steps

1) Create a folder for the scripts to save their output.
2) Add the root of the whole directory to your PYTHONPATH
3) Run "createExperiment.py" from the folder created in step 1 to create simulated data of the respective experiment.
4) Run "inferHyperparams.py" to infer all the hyperparameters later used by FGPGM.
5) Run "doFGPGM.py" to infer the parameters of the ODEs

## Implementing new Experiments

To test the algorithm on a different system of ODEs, create a new instance of the Experiment class located in the Experiment.py file. The provided experiments are saved in the Experiments folder. After creating a new experiment, appropriate main files need to be created. To train one set of hyperparameters for all states, the LotkaVolterra main files can be used as guidance. To train one set for each state dimension, the ProteinTransduction main files should be used.

## Implementing a new Kernel

To implement a new kernel, create a new instance of the Kernel class located in the Kernel.py file. The provided kernels are stored in the Kernels folder.
