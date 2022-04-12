from joblib import load, dump
import pandas as pd
from sklearn.metrics import accuracy_score, recall_score, precision_score, f1_score, matthews_corrcoef, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import os

trainpath = "/home/chuand/small_orf/data/train_feature/"
testpath = "/home/chuand/small_orf/data/test_feature/"

fileToBest = {}
fileToBest["c"] = "190:15:10"
fileToBest["ck"] = "160:13:10"
fileToBest["ckn"] = "200:19:10"
fileToBest["cn"] = "130:19:10"
fileToBest["k"] = "190:17:10"
fileToBest["kn"] = "170:13:20"
fileToBest["n"] = "190:11:10"
fileToBest["z"] = "170:13:10"
fileToBest["zc"] = "50:19:10"
fileToBest["zck"] = "170:19:10"
fileToBest["zckn"] = "200:15:10"
fileToBest["zcn"] = "150:13:10"
fileToBest["zk"] = "190:15:10"
fileToBest["zkn"] = "190:13:10"
fileToBest["zn"] = "60:15:20"

test_list = ["Ara", "Hum", "Mou", "Pro", "Bac"]
OUT = open("/home/chuand/small_orf/result/trainToOther", "w")
for i in fileToBest:
    parameters = fileToBest[i].split(":")
    n_estimators = int(parameters[0])
    max_depth = int(parameters[1])
    min_samples_split = int(parameters[2])
    f = os.path.join(trainpath, "train."+i)
    benchmark=pd.read_csv(f, sep=',')
    benchmark['Labels'] = (benchmark['Labels'] == "coding").astype(int)
    rows, columns = benchmark.shape
    X=benchmark.iloc[:,1:]
    Y=benchmark['Labels']
    X_train,y_train = X, Y
    rf=RandomForestClassifier(random_state=1, n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split)
    rf.fit(X_train, y_train)
    
    for each_test in test_list:
        test_f = os.path.join(testpath, each_test + "."  + i)
        test = pd.read_csv(test_f, sep=',')
        test['Labels'] = (test['Labels'] == "coding").astype(int)
        X=test.iloc[:,1:]
        Y=test['Labels']
        X_test,y_test = X ,Y
        y_prediction = rf.predict(X_test)
        acc=accuracy_score(y_test, y_prediction)
        recall=recall_score(y_test, y_prediction)
        precision=precision_score(y_test, y_prediction)
        f1=f1_score(y_test, y_prediction)
        mcc = matthews_corrcoef(y_test, y_prediction)
        outinfo = [i, each_test, ":".join(parameters), str(round(acc, 4)), str(round(recall, 4)), str(round(precision, 4)),str(round(f1, 4)), str(round(mcc, 4))]
        OUT.write("\t".join(outinfo) + "\n")
    del rf
OUT.close()
