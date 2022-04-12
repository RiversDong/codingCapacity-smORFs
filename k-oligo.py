from Bio import SeqIO
import pandas as pd
import sys
from functools import reduce
import numpy as np

k = int(sys.argv[1])
infp = sys.argv[2]
labelp = sys.argv[3]
infn = sys.argv[4]
labeln = sys.argv[5]
outfeature = sys.argv[6]

kmer = reduce(lambda x, y: [i + j for i in x for j in y], [['A', 'T', 'C', 'G']] * k)
def countf(kstrings, ks = kmer):
    strings_len = len(kstrings)
    res = []
    for i in kmer:
        i_frequency = round(kstrings.count(i)/strings_len, 4)
        res.append(i_frequency)
    return res

def getfeature(fasta_file, seq_type, k):    
    records = SeqIO.parse(fasta_file, "fasta")
    features = []
    labels = []
    for i in records:
        i_id = str(i.id)
        i_seq = str(i.seq)
        i_seq_len = len(i_seq) - k + 1
        kstrings = []
        for j in range(0, i_seq_len):
            kstrings.append(i_seq[j:j+k])
        labels.append(seq_type)
        res = countf(kstrings)
        features.append(res)     
    out_feature = np.array(features)
    return out_feature, labels

if __name__ == "__main__":
    codingFeatures, seq_labels= getfeature(infp, labelp, k)
    ncodingFeatures, ncodingLables = getfeature(infn, labeln, k)
    seq_features = np.vstack((codingFeatures, ncodingFeatures))
    seq_labels.extend(ncodingLables)
    df = pd.DataFrame(seq_features, index=seq_labels)
    df.index.name = "Labels"
    df.to_csv(outfeature)
    
    
    