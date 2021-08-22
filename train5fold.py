from sklearn.model_selection import GridSearchCV,StratifiedKFold
from sklearn.metrics import make_scorer, accuracy_score, recall_score, precision_score, f1_score, matthews_corrcoef, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import numpy as np
import os
import sys

path = "/home/chuand/small_orf/data/train_feature/"
files = os.listdir(path)
scoring = {'Accuracy':make_scorer(accuracy_score), 'Recall':make_scorer(recall_score), 'Precision':make_scorer(precision_score), "F1":make_scorer(f1_score), "MCC":make_scorer(matthews_corrcoef), "AUC":make_scorer(roc_auc_score)}    
for i in files:
    OUT=open("/home/chuand/small_orf/result/5fold"+ "_" + i, "w")
    feature = os.path.join(path, i)
    data = pd.read_csv(feature, sep=',')
    data['Labels'] = (data['Labels'] == "coding").astype(int)
    X=data.iloc[:,1:]
    Y=data['Labels']
    X_train,y_train = X ,Y
    rf=RandomForestClassifier(random_state=1)
    parameters = {'n_estimators':range(50,201,10),'max_depth':range(3,20,2),'min_samples_split':range(10,301,10)}
    #parameters = {'n_estimators':range(50,201,100),'max_depth':range(3,20,10),'min_samples_split':range(10,301,100)}
    kflod = StratifiedKFold(n_splits=5, random_state=1, shuffle=True)
    grid = GridSearchCV(rf,parameters, scoring=scoring, refit='F1', cv=kflod.split(X_train,y_train), n_jobs=10)
    grid.fit(X_train,y_train)
    results=grid.cv_results_
    #print(results.keys())
    accs=results["mean_test_Accuracy"]
    acc_std=results["std_test_Accuracy"]
    recalls=results["mean_test_Recall"]
    recall_std=results["std_test_Recall"]
    precisions=results["mean_test_Precision"]
    precision_std=results["std_test_Precision"]
    f1_score=results["mean_test_F1"]
    f1_score_std=results["std_test_F1"]
    mcc_score = results["mean_test_MCC"]
    mcc_std = results["std_test_MCC"]
    auc_score = results["mean_test_AUC"]
    auc_std = results["std_test_AUC"]
    params=results["params"]
    
    index_range=range(0, len(accs))
    for j in index_range:
        tree_number=str(params[j]["n_estimators"])
        depth=str(params[j]["max_depth"])
        min_samples_split=str(params[j]["min_samples_split"])
        acc=str(round(accs[j],4))
        recall=str(round(recalls[j],4))
        precision=str(round(precisions[j],4))
        f1=str(round(f1_score[j],4))
        mcc = str(round(mcc_score[j],4))
        auc = str(round(auc_score[j],4))
        info=["{0}:{1}:{2}".format(tree_number,depth,min_samples_split),acc,recall,precision,f1,mcc, auc]
        OUT.write(",".join(info)+"\n")
    del rf
    del grid
    OUT.close()


