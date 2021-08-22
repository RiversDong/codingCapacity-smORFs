import sys
import os
import pandas as pd
from sklearn.model_selection import GridSearchCV,StratifiedKFold
from sklearn.metrics import make_scorer, accuracy_score, recall_score, precision_score, f1_score, matthews_corrcoef, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import numpy as np

scoring = {'Accuracy':make_scorer(accuracy_score), 'Recall':make_scorer(recall_score), 'Precision':make_scorer(precision_score), "F1":make_scorer(f1_score), "MCC":make_scorer(matthews_corrcoef), "AUC":make_scorer(roc_auc_score)}    
def trainrf(X, Y, OUT):
    OUT = open(OUT, "w")
    X_train,y_train = X ,Y
    rf=RandomForestClassifier(random_state=1)
    parameters = {'n_estimators':range(50,201,10),'max_depth':range(3,20,2),'min_samples_split':range(10,301,10)}
    #parameters = {'n_estimators':range(50,201,50),'max_depth':range(3,20,2),'min_samples_split':range(10,301,50)}
    kflod = StratifiedKFold(n_splits=5, random_state=1, shuffle=True)
    grid = GridSearchCV(rf,parameters, scoring=scoring, refit='F1', cv=kflod.split(X_train,y_train), n_jobs=10)
    grid.fit(X_train,y_train)
    results=grid.cv_results_
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
    best_parameter = []
    tmp_mcc = 0
    best_prediction = []
    index_range=range(0, len(accs))
    for j in index_range:
        tree_number=str(params[j]["n_estimators"])
        depth=str(params[j]["max_depth"])
        min_samples_split=str(params[j]["min_samples_split"])
        acc = round(accs[j],4)
        recall = round(recalls[j],4)
        precision = round(precisions[j],4)
        f1= round(f1_score[j],4)
        mcc = round(mcc_score[j],4)
        auc = round(auc_score[j],4)
        info=["{0}:{1}:{2}".format(tree_number,depth,min_samples_split), str(acc), str(recall), str(precision), str(f1), str(mcc), str(auc)]
        OUT.write("\t".join(info) + "\n")
        if mcc > tmp_mcc:
            best_parameter = [tree_number,depth,min_samples_split]
            tmp_mcc = mcc
            best_prediction = [acc, recall, precision, f1, mcc, auc]
    OUT.close()
    del rf
    return best_parameter, best_prediction

'''
os.system("python .\\bin1\\feature.py")
# grid search finding the best parameters
files = os.listdir("D:\\small_orf\\data\\train_feature")
files = ["codon.csv", "codonNAR.csv", "zcurveCodon.csv", "zcurveCodonNAR.csv", "zcurveNAR.csv"]
for i in files:
    os.system("python .\\bin1\\train5fold.py {0}".format(i))
    
os.system("python .\\bin1\\figure.py")
os.system("python .\\bin1\\featureInd.py")
os.system("python .\\bin1\\modelPredict.py > D:\\small_orf\\result\ind.res")

# noncoding_Gene.fasta 和 Bac_noncoding.fasta 是一个数据集
NAR单独立
测试集 预测 测试集 热图
5fold_zcurve.csv 2415 200:13:10 0.8092
5fold_codon.csv 1949 180:11:20 0.7913
5fold_zcurveCodon.csv 1930 150:11:10 0.8164
5fold_zcurveNAR.csv 1949 180:11:20 0.8011
5fold_codonNAR.csv 1924 90:11:10 0.7898
5fold_zcurveCodonNAR.csv 3845 100:19:10 0.8232
'''
test = ["Ara", "Hum", "Mou", "Pro"]
outpath = "/home/chuand/small_orf/data/test_feature/"
#OUT = open("/home/chuand/small_orf/result/eachSpecies.cv", "w")
for i in test:
    zcknf = os.path.join(outpath, i + ".zckn")
    zckn = pd.read_csv(zcknf, sep = ',')
    zckn_data = zckn.iloc[:, 1:]
    zckn['Labels'] = (zckn['Labels'] == "coding").astype(int)
    Y = zckn['Labels']
    X = zckn_data
    result5fold = "zckn" + i
    best_parameter, best_prediction = trainrf(X,Y,result5fold)
    print(i, best_parameter, best_prediction)

    

