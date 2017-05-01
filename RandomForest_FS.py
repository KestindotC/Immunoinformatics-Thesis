## Script used to build the model for TCGA Colorectal Cancer 
## Immunogenomic and Tumor genomic analysis. (RandomForest version)
## Used in 2017, Academic Research Use.
## Author: Kestin Chang

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.model_selection import train_test_split,KFold
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score,confusion_matrix,roc_curve,auc

########################
#  Data Preprocessing  #
########################
print('Data Pre-process...')
df = pd.DataFrame.from_csv('EmCD8_Features_train.txt',header=0, sep='\t', index_col=False)
df_l = pd.DataFrame.from_csv('EmCD8_Labels_train.txt',header=0, sep='\t', index_col=False)
feature_name = df.columns.values

X = df.as_matrix()
y = df_l.as_matrix()
y = y[:, 0]
n_samples, n_features = X.shape
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.33)


########################
#  Feature Selection   #
########################
print('Feature Selection Start...')
fs_clf = ExtraTreesClassifier(n_estimators = 16000, max_features = "sqrt")
fs_clf.fit(X_train, y_train)
ranking = fs_clf.feature_importances_

model = SelectFromModel(fs_clf, prefit=True, threshold=0.001)
X_important_train = model.transform(X_train)
X_important_test = model.transform(X_test)
X_important_all = model.transform(X)

print("Extract %d important features!!"%X_important_train.shape[1])

global important_features
important_features = feature_name[model.get_support(indices=True)]
# Build default dictionary for final importance report
global final_features
from collections import defaultdict
final_features = defaultdict(lambda: 0.0)
def feature_counter(tmp_imp):
    for f_ind in range(len(tmp_imp)):
        final_features[important_features[f_ind]] += tmp_imp[f_ind]


#############################
#  Random Forest Classifier #
#############################
from scipy import interp
print("Classify Start...")

global line_col
line_col = ['crimson','green','yellow','navy','grey']
acc = []
cv = KFold(n_splits=5)
clf = RandomForestClassifier(n_estimators = 100,max_features=0.5)
fold = 1
mean_tpr = 0.0
mean_fpr = np.linspace(0, 1, 100)

# Start Cross Validation using our important features
for train_index, test_index in cv.split(X_important_test):
    #print("TRAIN:", train_index, "TEST:", test_index)
    X_cvtrain, X_cvtest = X_important_test[train_index], X_important_test[test_index]
    y_cvtrain, y_cvtest = y_test[train_index], y_test[test_index]
    clf.fit(X_cvtrain, y_cvtrain)
    prediction_score = clf.predict_proba(X_cvtest)  # Get predict possibility of classes
    predictions = clf.predict(X_cvtest)             # Get predict classes
    tmp_imp = clf.feature_importances_
    feature_counter(tmp_imp)
    fp, tp, thresholds = roc_curve(y_cvtest, prediction_score[:,1]) # For class 1
    mean_tpr += interp(mean_fpr, fp, tp)
    mean_tpr[0] = 0.0
    roc_auc = auc(fp, tp)
    plt.plot(fp, tp, lw=1, linestyle='--',
             color=line_col.pop(),
             label='ROC fold %d (area = %0.2f)' %(fold, roc_auc))
    #acc = cross_val_score(clf, X_cvtrain, y_cvtrain, cv=10, scoring='accuracy') 
    #Alternative expression for simply calculate cv accuracy
    acc.append(accuracy_score(y_cvtest,predictions))
    fold += 1

######################
#  Model Evaluation  #
######################

print('Mean Accuracy:%f'%np.mean(acc))

import csv
df_final = pd.DataFrame(X_important_all,columns=important_features)
df_final.to_csv('./Out_ImportantFeaturesMatrix_EmCD8_train.txt',sep="\t",index=False)
w = csv.writer(open("Out_RankFeatures_EmCD8_train.csv", "w"))
for key, val in final_features.items():
    w.writerow([key,val])


#######################
# Visualization Model #
#######################
from matplotlib.backends.backend_pdf import PdfPages
# For ROC Curve
pp = PdfPages('ROC_curve.pdf')
mean_tpr /= cv.get_n_splits(X_important_test, y_test)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
plt.plot(mean_fpr, mean_tpr, color='skyblue', linestyle='solid',label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)
plt.title('Receiver Operating Characteristic Curve')
plt.xlabel('1-Specificity')
plt.ylabel('Sensitivity')
plt.legend(loc="lower right")
#plt.show()
plt.savefig(pp, format='pdf')
pp.savefig()
pp.close()

# For Feature Importance 
# @@@@@ Transfer to R ggplot
#import seaborn as sns
#feature_table = [(v, k) for k, v in final_features.iteritems()]
#data_imp = pd.DataFrame(feature_table,columns=["importance","features"])
#data = data_imp[0:100]
#sns.set(style="whitegrid")
#fig = sns.stripplot(x=data['importance'],y=data['features'],size=4, orient="h",palette="Reds_r", edgecolor="gray")
#fig.set(xlabel="Feature Importance", ylabel="")
#plt.show()

