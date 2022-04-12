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
    positives = []
    negatives = []
    train_p1 = "/home/chuand/small_orf/data/trainAndTest/train/coding_sORF_seq_uniq_v1.fa"
    train_n1 = "/home/chuand/small_orf/data/trainAndTest/train/noncoding_sORF_seq_uniq_v1.fa"
    positives.append(train_p1)
    negatives.append(train_n1)
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
    mixed_features.to_csv("/home/chuand/small_orf/data/pro_train.csv")

if __name__ == "__main__":
    mergeSpecesAndFeature()
