from Bio import SeqIO
import numpy as np
import pandas as pd
import sys

def codonFreq(seq):
    phase_1 = range(0, len(seq), 3)
    codons = []
    freqList = []
    for i in phase_1:
        codon = seq[i:i+3]
        codons.append(codon)
    nums = len(codons)
    for i in codon64:
        freq = codons.count(i)/nums
        freqList.append(freq)
    return freqList

def singleZ(subSeq):
    fa = subSeq.count("A")/len(subSeq)
    fg = subSeq.count("G")/len(subSeq)
    fc = subSeq.count("C")/len(subSeq)
    ft = subSeq.count("T")/len(subSeq)
    x = (fa + fg) - (fc + ft)
    y = (fa + fc) - (fg + ft)
    z = (fa + ft) - (fg + fc)
    return x, y, z

def dinucleotideZ(subSeq1, subSeq2):
    bases = ["A", "C", "G", "T"]
    dibase_list = []
    dinucleotideZ_feature = []
    for i in range(0, len(subSeq1)):
        dibase_list.append(subSeq1[i]+subSeq2[i])
    dibase_list_len = len(dibase_list)
    for i in bases:
        XA = dibase_list.count(i + "A")/dibase_list_len
        XG = dibase_list.count(i + "G")/dibase_list_len
        XC = dibase_list.count(i + "C")/dibase_list_len
        XT = dibase_list.count(i + "T")/dibase_list_len
        XX = (XA + XG) - (XC + XT)
        YX = (XA + XC) - (XG + XT)
        ZX = (XA + XT) - (XG + XC)
        dinucleotideZ_feature.append(XX)
        dinucleotideZ_feature.append(YX)
        dinucleotideZ_feature.append(ZX)
    return dinucleotideZ_feature

def main(infile, label):
    records = SeqIO.parse(infile, "fasta")
    zcurve_variable = []
    row_name = []
    codon_feature = []
    zcurve_codon = []
    zcurve_codon_nar = []
    zcurve_nar = []
    codon_nar = []
    nar=open("D:\\small_orf\\data\\TR_CPPred-sORF.csv").read().split("\n")
    nar = nar[1:-1]
    nar2feature = {}
    for i in nar:
        info = i.split(",")
        geneid = info[0]
        nar_string = info[1:-1]
        nar_feature = [float(j) for j in nar_string]
        if geneid not in nar2feature:
            nar2feature[geneid] = [nar_feature]
        else:
            nar2feature[geneid].append(nar_feature)

    for i in records:
        row_name.append(label)
        i_id = str(i.id)
        i_sequence = str(i.seq)
        phase_1 = range(0, len(i_sequence), 3)
        phase_2 = range(1, len(i_sequence), 3)
        phase_3 = range(2, len(i_sequence), 3)
        phase_1_seq = "".join([i_sequence[j] for j in phase_1])
        phase_2_seq = "".join([i_sequence[j] for j in phase_2])
        phase_3_seq = "".join([i_sequence[j] for j in phase_3])
        x1, y1, z1 = singleZ(phase_1_seq)
        x2, y2, z2 = singleZ(phase_2_seq)
        x3, y3, z3 = singleZ(phase_3_seq)
        dinucleotideZ_in12 = dinucleotideZ(phase_1_seq, phase_2_seq)
        dinucleotideZ_in23 = dinucleotideZ(phase_2_seq, phase_3_seq)
        dinucleotideZ_in13 = dinucleotideZ(phase_1_seq, phase_3_seq)
        tmp = [x1, y1, z1, x2, y2, z2, x3, y3, z3]
        tmp.extend(dinucleotideZ_in12)
        tmp.extend(dinucleotideZ_in23)
        tmp.extend(dinucleotideZ_in13)
        tmp1 = []
        tmp1.extend(tmp)
        zcurve_variable.append(tmp)
        freqList = codonFreq(i_sequence)
        codon_feature.append(freqList)
        tmp5 = []
        tmp5.extend(tmp)
        tmp5.extend(freqList)
        zcurve_codon.append(tmp5)
        if label == "coding":
            narf = nar2feature[i_id][0]
            tmp3 = []
            tmp3.extend(tmp)
            tmp3.extend(freqList)
            tmp3.extend(narf)
            zcurve_codon_nar.append(tmp3)
        if label == "noncoding":
            narf = nar2feature[i_id][0]
            tmp3 = []
            tmp3.extend(tmp)
            tmp3.extend(freqList)
            tmp3.extend(narf)
            zcurve_codon_nar.append(tmp3)
        tmp1.extend(narf)
        zcurve_nar.append(tmp1)

        tmp2 = []
        tmp2.extend(freqList)
        tmp2.extend(narf)
        codon_nar.append(tmp2)
    return row_name, np.array(zcurve_variable), np.array(codon_feature), np.array(zcurve_codon), np.array(zcurve_codon_nar), np.array(zcurve_nar), np.array(codon_nar)

if __name__ == "__main__":
    bases = ["A", "C", "G", "T"]
    codon64=[]

    for base1 in bases:
        for base2 in bases:
            for base3 in bases:
                codon123 = base1 + base2 + base3
                codon64.append(codon123)
    label, z, c, zc, zcn, zn, cn= main("D:\\small_orf\\data\\trainAndTest\\train\\coding_sORF_seq_uniq_v1.fa", "coding")
    label_1, z_1, c_1, zc_1, zcn_1, zn_1, cn_1= main("D:\\small_orf\\data\\trainAndTest\\train\\noncoding_sORF_seq_uniq_v1.fa", "noncoding")
    zcurve_feature =  np.vstack((z, z_1))
    codon_feature =  np.vstack((c, c_1))
    zcurve_codon = np.vstack((zc, zc_1))
    zcurve_codon_nar = np.vstack((zcn, zcn_1))
    zcurve_nar = np.vstack((zn, zn_1))
    codon_nar = np.vstack((cn, cn_1))
    label.extend(label_1)

    df_z = pd.DataFrame(zcurve_feature, index=label)
    df_z.index.name = "Labels"
    df_c = pd.DataFrame(codon_feature, index=label)
    df_c.index.name = "Labels"
    df_zc = pd.DataFrame(zcurve_codon, index=label)
    df_zc.index.name = "Labels"
    df_zcn = pd.DataFrame(zcurve_codon_nar, index=label)
    df_zcn.index.name = "Labels"
    df_zn = pd.DataFrame(zcurve_nar, index=label)
    df_zn.index.name = "Labels"
    df_cn = pd.DataFrame(codon_nar, index=label)
    df_cn.index.name = "Labels"

    df_z.to_csv("D:\\small_orf\\data\\train_feature\\zcurve.csv")
    df_c.to_csv("D:\\small_orf\\data\\train_feature\\codon.csv")
    df_zc.to_csv("D:\\small_orf\\data\\train_feature\\zcurveCodon.csv")
    df_zn.to_csv("D:\\small_orf\\data\\train_feature\\zcurveNAR.csv")
    df_cn.to_csv("D:\\small_orf\\data\\train_feature\\codonNAR.csv")
    df_zcn.to_csv("D:\\small_orf\\data\\train_feature\\zcurveCodonNAR.csv")
    

