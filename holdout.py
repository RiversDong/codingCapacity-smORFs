import pandas as pd
import random
from sklearn.metrics import make_scorer, accuracy_score, recall_score, precision_score, f1_score, matthews_corrcoef, roc_auc_score
from sklearn.ensemble import RandomForestClassifier

f = "/home/chuand/small_orf/data/train_feature/train.zckn"
data = pd.read_csv(f, sep=",")

labels = data["Labels"]
p_labels = labels[labels=="coding"]
n_labels = labels[labels=="noncoding"]
n_labels.reset_index(drop=True, inplace=True)

features = data.iloc[:,1:]
p_features = features[labels=="coding"]
n_features = features[labels=="noncoding"]
n_features.reset_index(drop=True, inplace=True)

p_row = p_features.shape[0]
p_row_index = range(0, p_row)
n_row = n_features.shape[0]
n_row_index = range(0, n_row)

# random index for each of the two data types
for i in [0.1, 0.2, 0.3]:
    OUT = open("/home/chuand/small_orf/result/holdout_"+str(i), "w")
    for j in range(0,100):
        p_index = random.sample(p_row_index, int(len(p_row_index)*i))
        n_index = random.sample(n_row_index, int(len(n_row_index)*i))
        
        print(i, p_index)
        print(i, n_index)
        print("\n")
        
        remaining_p = list(set(p_row_index).difference(set(p_index)))
        remaining_n = list(set(n_row_index).difference(set(n_index)))
        
        # select random samples among the two data types
        predict_label_1 = p_labels[p_index]
        predict_label_2 = n_labels[n_index]
        predict_label = pd.concat([predict_label_1,predict_label_2],axis=0)
        predict_label = (predict_label=="coding").astype(int)
        predict_feature1 = p_features.iloc[p_index,:]
        predict_feature2 = n_features.iloc[n_index,:]
        predict_feature = pd.concat([predict_feature1,predict_feature2], axis=0)

        training_label_1 = p_labels[remaining_p]
        training_label_2 = n_labels[remaining_n]
        training_label = pd.concat([training_label_1, training_label_2], axis=0)
        training_label = (training_label=="coding").astype(int)
        train_feature1 = p_features.iloc[remaining_p,:]
        train_feature2 = n_features.iloc[remaining_n,:]
        train_feature = pd.concat([train_feature1, train_feature2], axis=0)

        # train momdel and predict a random samples
        rf=RandomForestClassifier(random_state=1, n_estimators=200, max_depth=15, min_samples_split=10)
        rf.fit(train_feature, training_label)
        y_prediction = rf.predict(predict_feature)
        acc=accuracy_score(predict_label, y_prediction)
        recall=recall_score(predict_label, y_prediction)
        precision=precision_score(predict_label, y_prediction)
        f1=f1_score(predict_label, y_prediction)
        mcc = matthews_corrcoef(predict_label, y_prediction)
        info = [str(i), str(round(acc, 4)), str(round(recall, 4)), str(round(precision, 4)), str(round(f1, 4)), str(round(mcc,4))]
        OUT.write("\t".join(info) + "\n")
        #break
    OUT.close()


