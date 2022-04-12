import os
import numpy as np
import sys
from featurePackage import *
from Bio import SeqIO
import pandas as pd
from pandas import Series
from sklearn.model_selection import GridSearchCV,StratifiedKFold
from sklearn.metrics import make_scorer, accuracy_score, recall_score, precision_score, f1_score, matthews_corrcoef, roc_auc_score
from sklearn.ensemble import RandomForestClassifier

def mergeSpecesAndFeature():
    path = "/home/chuand/small_orf/data/trainAndTest/test/"
    species = ["Bac", "Pro"]
    positives = []
    negatives = []
    for i in species:
        p = os.path.join(path, i + "_coding.fasta")
        n = os.path.join(path, i + "_nocoding.fasta")
        positives.append(p)
        negatives.append(n)
    mixed_features = pd.DataFrame()
    ids = []
    lables = []
    files = []
    for j in positives:
        base_j = os.path.basename(j)
        features = []
        py2 = "/home/chuand/.bin/miniconda3/envs/py2/bin/python2.7"
        Integrated_Hexamer = "/home/chuand/small_orf/bin1/CPPred-sORF/Hexamer/Integrated_Hexamer.txv"
        Integrated = "/home/chuand/small_orf/bin1/CPPred-sORF/Model/Integrated.range"
        Integrated_model = "/home/chuand/small_orf/bin1/CPPred-sORF/Model/Integrated.model"
        cmd = "{0} /home/chuand/small_orf/bin1/CPPred-sORF/CPPred-sORF.py -i {1} -hex {2} -r {3} -mol {4} -spe Integrated -o sORF.csv".format(py2, j, Integrated_Hexamer, Integrated, Integrated_model)
        os.system(cmd)
        nar = pd.read_csv("sORF.csv", sep = "\t")
        nar_data = nar.iloc[:,1:-2]
        records = SeqIO.parse(j, "fasta")
        for i in records:
            i_id = str(i.id)
            ids.append(i_id)
            i_seq = str(i.seq)
            feature = []
            zvariable = zcurve(i_seq)
            kmer_freq = kmer(i_seq, 4)
            codon = codonFreq(i_seq)
            T_ratio = Tfreq(i_seq)
            feature.extend(zvariable)
            feature.extend(codon)
            feature.extend(kmer_freq)
            feature.append(T_ratio)
            features.append(feature)
            lables.append("coding")
            files.append(base_j)
        features = np.array(features)
        j_df = pd.DataFrame(features)
        j_nar = pd.concat([j_df, nar_data], axis=1)
        mixed_features = pd.concat([mixed_features, j_nar], axis=0)
        #break
    for j in negatives:
        base_j = os.path.basename(j)
        features = []
        py2 = "/home/chuand/.bin/miniconda3/envs/py2/bin/python2.7"
        Integrated_Hexamer = "/home/chuand/small_orf/bin1/CPPred-sORF/Hexamer/Integrated_Hexamer.txv"
        Integrated = "/home/chuand/small_orf/bin1/CPPred-sORF/Model/Integrated.range"
        Integrated_model = "/home/chuand/small_orf/bin1/CPPred-sORF/Model/Integrated.model"
        cmd = "{0} /home/chuand/small_orf/bin1/CPPred-sORF/CPPred-sORF.py -i {1} -hex {2} -r {3} -mol {4} -spe Integrated -o sORF.csv".format(py2, j, Integrated_Hexamer, Integrated, Integrated_model)
        os.system(cmd)
        nar = pd.read_csv("sORF.csv", sep = "\t")
        nar_data = nar.iloc[:,1:-2]
        records = SeqIO.parse(j, "fasta")
        for i in records:
            i_id = str(i.id)
            ids.append(i_id)
            i_seq = str(i.seq)
            feature = []
            zvariable = zcurve(i_seq)
            kmer_freq = kmer(i_seq, 4)
            codon = codonFreq(i_seq)
            T_ratio = Tfreq(i_seq)
            feature.extend(zvariable)
            feature.extend(codon)
            feature.extend(kmer_freq)
            feature.append(T_ratio)
            features.append(feature)
            lables.append("noncoding")
            files.append(base_j)
        features = np.array(features)
        j_df = pd.DataFrame(features)
        j_nar = pd.concat([j_df, nar_data], axis=1)
        mixed_features = pd.concat([mixed_features, j_nar], axis=0)
        #break
    mixed_features["ID"] = ids
    mixed_features["file"] = files
    mixed_features.index = Series(lables)
    mixed_features.index.name = "Labels"
    mixed_features.to_csv("/home/chuand/small_orf/data/Pro-6318Bac.csv")

scoring = {'Accuracy':make_scorer(accuracy_score), 'Recall':make_scorer(recall_score), 'Precision':make_scorer(precision_score), "F1":make_scorer(f1_score), "MCC":make_scorer(matthews_corrcoef), "AUC":make_scorer(roc_auc_score)}
def trainrf(X, Y, OUT):
    OUT = open(OUT, "w")
    X_train,y_train = X ,Y
    rf=RandomForestClassifier(random_state=1)
    parameters = {'n_estimators':range(50,201,10),'max_depth':range(3,20,2),'min_samples_split':range(10,301,10)}
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
def fiveFold():
    f = "/home/chuand/small_orf/data/AllSpeciesMixedFeature.csv"
    zckn = pd.read_csv(f, sep = ',')
    zckn['Labels'] = (zckn['Labels'] == "coding").astype(int)
    Y = zckn['Labels']
    X = zckn.iloc[:,1:-2]
    result5fold = "allSpecies"
    best_parameter, best_prediction = trainrf(X,Y,result5fold)
    # nohup python gridSearchMixed.py > ../result/mixedSpeciesFeature.cv&
    print(best_parameter, best_prediction)

def prediction():
    n_estimators = 170 
    max_depth = 19
    min_samples_split = 20
    species = ["Pro", "Bac"]
    pfiles = [i+"_coding.fasta" for i in species]
    nfiles = [i+"_nocoding.fasta" for i in species]
    data = pd.read_csv("/home/chuand/small_orf/data/AllSpeciesMixedFeature.csv")
    data['Labels'] = (data['Labels'] == "coding").astype(int)
    Y_train = data['Labels']
    X_train = data.iloc[:,1:-2]
    rf=RandomForestClassifier(random_state=1, n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split)
    rf.fit(X_train, Y_train)
    predict_data = pd.read_csv("/home/chuand/small_orf/data/Pro-6318Bac.csv")
    predict_data['Labels'] = (predict_data['Labels'] == "coding").astype(int)
    for i in range(0, 2):
        c = pfiles[i]
        nc = nfiles[i]
        prediction = predict_data[(predict_data["file"] == nc) | (predict_data["file"] == c)]
        Y_ = prediction["Labels"]
        X_prediction = prediction.iloc[:,1:-2]
        print(X_prediction.shape)
        Y_prep = rf.predict(X_prediction)
        acc=round(accuracy_score(Y_, Y_prep),4)
        recall=round(recall_score(Y_, Y_prep),4)
        precision=round(precision_score(Y_, Y_prep),4)
        f1=round(f1_score(Y_, Y_prep),4)
        mcc = round(matthews_corrcoef(Y_, Y_prep),4)
        print(c.split("_")[0], acc, recall, precision, f1, mcc)

#mergeSpecesAndFeature()    
prediction()



















