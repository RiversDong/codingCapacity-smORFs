from Bio import SeqIO
import numpy as np
import pandas as pd
from pandas import Series
import sys
from featurePackage import *
import os
from itertools import combinations

def GetFeature(infile, label):
    # features from NAR
    py2 = "/home/chuand/.bin/miniconda3/envs/py2/bin/python2.7"
    Integrated_Hexamer = "/home/chuand/small_orf/bin1/CPPred-sORF/Hexamer/Integrated_Hexamer.txv"
    Integrated = "/home/chuand/small_orf/bin1/CPPred-sORF/Model/Integrated.range"
    Integrated_model = "/home/chuand/small_orf/bin1/CPPred-sORF/Model/Integrated.model"
    cmd = "{0} /home/chuand/small_orf/bin1/CPPred-sORF/CPPred-sORF.py -i {1} -hex {2} -r {3} -mol {4} -spe Integrated -o sORF.csv".format(py2, infile, Integrated_Hexamer, Integrated, Integrated_model)
    os.system(cmd)
    nar = pd.read_csv("sORF.csv", sep = "\t")
    nar_df = nar.iloc[:,1:-2]
    Labels = []
    records = SeqIO.parse(infile, "fasta")
    zcurves = []
    codons = []
    kmers = []
    Ts = []
    for i in records:
        Labels.append(label)
        i_seq = str(i.seq)
        zvariable = zcurve(i_seq)
        zcurves.append(zvariable)
        codon_freq = codonFreq(i_seq)
        codons.append(codon_freq)
        kmer_freq = kmer(i_seq, 4)
        kmers.append(kmer_freq)
        T_ratio = Tfreq(i_seq)
        Ts.append(T_ratio)

    # Z curve
    zcurves_df = pd.DataFrame(zcurves, index=Series(Labels))
    zcurves_df.index.name = "Labels"
    # Codon frequency
    codons_df = pd.DataFrame(codons, index=Series(Labels))
    codons_df.index.name = "Labels"
    # kmers
    kmers_df = pd.DataFrame(kmers, index=Series(Labels))
    kmers_df.index.name = "Labels"
    # NAR
    nar_df["Ts"] = Ts 
    nar_df.index = Series(Labels)
    nar_df.index.name = "Labels"
    
    return zcurves_df, codons_df, kmers_df, nar_df


if __name__ == "__main__":
    #prefix = ["Ara", "Hum", "Mou","Pro"]
    #prefix = ["Ara"]
    prefix = ["Bac"]
    test_path = "/home/chuand/small_orf/data/trainAndTest/test"
    feature_out = "/home/chuand/small_orf/data/test_feature/"
    feature_type = ["zcurves", "codons", "kmers", "nar"]
    suffix = {"zcurves":"z",  "codons":"c", "kmers":"k", "nar":"n"}

    for i in prefix:
        codingf = os.path.join(test_path, i+"_coding.fasta")
        noncodingf = os.path.join(test_path, i+"_nocoding.fasta")
        zcurves_c, codons_c, kmers_c, nar_c = GetFeature(codingf, "coding")
        zcurves_nc, codons_nc, kmers_nc, nar_nc = GetFeature(noncodingf, "noncoding")
        # all 4 type features
        positives = pd.concat([zcurves_c, codons_c, kmers_c, nar_c], axis=1)
        negatives = pd.concat([zcurves_nc, codons_nc, kmers_nc, nar_nc], axis=1)
        features = pd.concat([positives, negatives], axis=0)
        outpath = os.path.join(feature_out, "{0}.zckn".format(i))
        features.to_csv(outpath)

        # each of feature type
        for j in feature_type:
            feature_type_c = eval(j+"_c")
            feature_type_nc = eval(j+"_nc")
            features = pd.concat([feature_type_c, feature_type_nc], axis=0)
            outpath = os.path.join(feature_out, i+"."+suffix[j])
            features.to_csv(outpath)

        # two combinations among the four types
        for comNum in range(2,4):
            twoFeatures = combinations(feature_type, comNum)
            for j in twoFeatures:
                j_list = list(j)
                feature_type_c = [eval(k+"_c") for k in j_list]
                feature_type_nc = [eval(k+"_nc") for k in j_list]
                positives = pd.concat(feature_type_c, axis=1)
                negatives = pd.concat(feature_type_nc, axis=1)
                features = pd.concat([positives, negatives], axis=0)
                tmpsuffix = [suffix[k] for k in j_list]
                outpath = os.path.join(feature_out, "{0}.".format(i) + "".join(tmpsuffix))
                features.to_csv(outpath)





