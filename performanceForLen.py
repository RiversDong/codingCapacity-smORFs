import os
from sklearn.metrics import make_scorer, accuracy_score, recall_score, precision_score, f1_score, matthews_corrcoef, roc_auc_score, confusion_matrix

def getResult(infile):
    f = open(infile).read().split("\n")
    f = f[1:-1]
    results = []
    for i in f:
        info = i.split("\t")
        res = int(info[1])
        results.append(res)
    return results

path = "/home/chuand/small_orf/data/splitFasta/"
#path = "/home/chuand/small_orf/data/fly_difflen/"
# remember to change the model name
# when changing the file name
#files = ["split_eu_"]
files = ["split_pro_"]
#files = ["Fru_"]
for ii in files:
    #for j in ["1-100", "100-200", "200-300"]:
    for j in ["100", "200", "300"]:
        sample_labels = []
        # coding = os.path.join(path, ii+"coding_"+j+"_cdhit.fasta")
        coding = os.path.join(path, ii+"coding"+j)
        # noncoding = os.path.join(path, ii+"noncoding_"+j+"_cdhit.fasta")
        noncoding = os.path.join(path, ii+"noncoding"+j)
        # pro or eu 
        cmd = "python /home/chuand/small_orf/sORFPredictor/sORFPredictor.py\
                -i {0} -o {1} -m pro ".format(coding, "resultP")
        os.system(cmd)
        predictions = getResult("resultP.prediction")
        positives = [1 for i in range(0,len(predictions))]
        sample_labels.extend(positives)

        cmd = "python /home/chuand/small_orf/sORFPredictor/sORFPredictor.py\
                -i {0} -o {1} -m pro ".format(noncoding, "resultN")
        os.system(cmd)
        predict_negative = getResult("resultN.prediction")
        predictions.extend(predict_negative)
        negatives = [0 for i in range(0,len(predict_negative))]
        sample_labels.extend(negatives)
        ACC=round(accuracy_score(sample_labels, predictions),4)
        Recall=round(recall_score(sample_labels, predictions),4)
        Precision=round(precision_score(sample_labels, predictions),4)
        F1=round(f1_score(sample_labels, predictions),4)
        MCC = round(matthews_corrcoef(sample_labels, predictions),4)
        cnf_matrix = confusion_matrix(sample_labels, predictions)
        TP = cnf_matrix[0,0]
        FP = cnf_matrix[0,1]
        FN = cnf_matrix[1,0]
        TN = cnf_matrix[1,1]
        print(TP, TN, FP,FN)
        Sp = float(TN)/(TN+FP)
        print(j,TP, TN, FP,FN, Recall, Sp, ACC, MCC, Precision, Recall)