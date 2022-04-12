import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
path = "D:\\small_orf\\result\\"
files = ["5fold_zcurve.csv", "5fold_codon.csv", "5fold_zcurveCodon.csv", "5fold_zcurveNAR.csv", "5fold_codonNAR.csv", "5fold_zcurveCodonNAR.csv"]
#files = ["5fold_zcurve.csv"]
fig, axs = plt.subplots(6, 1)
i_index = 0
for i in files:
    result = os.path.join(path, i)
    data = pd.read_csv(result, sep=",", header=None)
    parameter = data[0]
    accuracy = data[1]
    recall = data[2]
    precision = data[3]
    f1 = data[4]
    mcc = data[5]
    auc = data[6]
    max_index = mcc.idxmax()
    max_mcc = mcc[max_index]
    print(i,max_index, parameter[max_index], max_mcc)
    x = range(0, len(parameter))
    axs[i_index].plot(x, accuracy, '-', color='lightgrey', label='ACC', linewidth=0.5)
    axs[i_index].plot(x, mcc, '-', label='MCC', linewidth=0.5)
    axs[i_index].scatter(max_index, max_mcc, color="red", label='MCC='+str(max_mcc))   
    axs[i_index].set_xticklabels([])
    axs[i_index].set_xticks([])
    axs[i_index].set_yticklabels([])
    axs[i_index].set_yticks([])
    i_index = i_index + 1
plt.tight_layout()
fig.savefig("D:\\small_orf\\Figure\\grid.jpeg", dpi=300)