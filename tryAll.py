from joblib import load, dump
import pandas as pd
from sklearn.metrics import accuracy_score, recall_score, precision_score, f1_score, matthews_corrcoef, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import os

def mergeFeature(modelFeature1, modelFeature2):
    Feature1 = pd.read_csv(modelFeature1)
    Feature1_data = Feature1.iloc[:, 1:]
    Feature1['Labels'] = (Feature1['Labels'] == "coding").astype(int)
    Feature2 = pd.read_csv(modelFeature2)
    Feature2_data = Feature2.iloc[:,1:]
    X = pd.concat([Feature1_data, Feature2_data], axis=1)
    Y = Feature1["Labels"]
    return X, Y

path = "/home/chuand/small_orf/data/test_feature"
apecies2parameters = {"Ara":"190:15:10", "Hum":"200:19:20", "Mou":"200:19:10", "Pro":"200:19:10"}
out = open("/home/chuand/small_orf/result/eachToEach", "w")
for i in apecies2parameters:
    modelSpecies = i
    #modelf1 = os.path.join(path, i + "_zcurveCodonNAR.csv")
    #modelf2 = os.path.join(path, i + "_kmer.csv")
    #X,Y = mergeFeature(modelf1, modelf2)
    i_path = os.path.join(path, i + ".zckn")
    Feature1 = pd.read_csv(i_path)
    Feature1['Labels'] = (Feature1['Labels'] == "coding").astype(int)
    Y = Feature1["Labels"]
    X = Feature1.iloc[:, 1:]
    modelPara = apecies2parameters[i].split(":")
    n_estimators = int(modelPara[0])
    max_depth = int(modelPara[1])
    min_samples_split = int(modelPara[2])
    rf=RandomForestClassifier(random_state=1,\
            n_estimators=n_estimators,\
            max_depth=max_depth,\
            min_samples_split=min_samples_split)
    rf.fit(X,Y)
    testSpecies = []
    for j in apecies2parameters:
        if modelSpecies != j:
            testSpecies.append(j)
    for j in testSpecies:
        #predictf1 = os.path.join(path, j + "_zcurveCodonNAR.csv")
        #predictf2 = os.path.join(path, j + "_kmer.csv")
        #X_j, Y_j = mergeFeature(predictf1, predictf2)
        predictf = os.path.join(path, j + ".zckn")
        predict_data = pd.read_csv(predictf)
        predict_data["Labels"] = (predict_data['Labels'] == "coding").astype(int)
        Y_j = predict_data["Labels"]
        X_j = predict_data.iloc[:,1:]
        Y_prep = rf.predict(X_j)
        acc=round(accuracy_score(Y_j, Y_prep),4)
        recall=round(recall_score(Y_j, Y_prep),4)
        precision=round(precision_score(Y_j, Y_prep),4)
        f1=round(f1_score(Y_j, Y_prep),4)
        mcc = round(matthews_corrcoef(Y_j, Y_prep),4)
        out.write(i+"\t"+j+"\t")
        out.write("%f\t%f\t%f\t%f\t%f\n"%(acc, recall, precision, f1, mcc))
out.close()

