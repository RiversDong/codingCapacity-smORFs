# sORF
Python scripts to extract sequence-derived features and perform RF-based validation and prediction

## featureCombineTrain.py 
Get all features of train dataset

## train5fold.py
Grid search based on five-fold validation in Prokaryotic species so that I can validate the best parameters (tree number, depth, sample number for splitting) of RF-based model for the prediction

## modelPredict.py 
Training and cross-species prediction and obtain the performance

## featureCombine.py 
Get all the features of testing datasets

## main.py 
Grid search based on five-fold validation in eukaryotic species so that I can validate the best parameter of RF model for the prediction

## tryAll.py
Cross-species prediction (One to one)

## gridSearchMixed.py
Merge all kind of features of all species and perform cross-species prediction (Four to one)

## featurePackage.py
Including some functions to calculate the sequence-derived features, which is called by the above scripts
