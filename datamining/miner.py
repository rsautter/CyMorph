import pandas as pd
import sys
import numpy as np
from sklearn.model_selection import StratifiedKFold as StratifiedKFold
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, roc_auc_score
from sklearn import svm as svm
from sklearn.neural_network import MLPClassifier
from imblearn.combine import SMOTETomek
from imblearn.over_sampling import SMOTE

#Valores default
gal = pd.read_csv("merged.csv")
kfolds = 10
tech = ["MLP1","MLP2","SVM1","SVM2"]
variables = ["sA3","sA2","sS2","sS3","sH","sGa"]
smoteIt = True
smoteTomekIt = True
k_neighbors = 5
aim = "Zoo1"
savetxt = "ConfusionMetrics.txt"
metrics = ["roc","acc","conf", "prec"]
ell,sp = "E","S"

##Falta interface com usuario


#Seleciona o set de treinamento e teste
sel = StratifiedKFold(n_splits=kfolds)
data = gal[variables+[aim]].dropna()

if smoteIt and not(smoteTomekIt):
	oversamp = SMOTE('minority',k_neighbors=k_neighbors)
	a,b = oversamp.fit_sample(data.loc[:,variables],data.loc[:,aim])
	print("Smoted","Sp:",len(np.where(b==sp)[0]),"Ell:",len(np.where(b==ell)[0]))
	print(a)
	a2 = np.column_stack((a,b))
	data = pd.DataFrame(data=a2,columns=data.columns.values)
elif smoteTomekIt:
	oversamp = SMOTETomek('auto')
	a,b = oversamp.fit_sample(data.loc[:,variables],data.loc[:,aim])
	print("Smoted","Sp:",len(np.where(b==sp)[0]),"Ell:",len(np.where(b==ell)[0]))
	print(a)
	a2 = np.column_stack((a,b))
	data = pd.DataFrame(data=a2,columns=data.columns.values)

precEllList, precSpList, accList, classifierList = [],[],[],[]
for technique in tech:
	print(technique)
	for train, test in sel.split(data[variables],data[aim]):
		dtrain = data.loc[train,variables].dropna()
		dknown = data.loc[train,aim].dropna()	
		dtest = data.loc[test,variables].dropna()
		dtknown = data.loc[test,aim].dropna()
		
		if technique == "SVM1":		
			classifier = svm.SVC()
		if technique == "SVM2":
		if technique == "MLP1":		
			print("Layer Size",data.shape[1])
			classifier = MLPClassifier(activation='logistic',hidden_layer_sizes=(data.shape[1],))
		if technique == "MLP2":		
			print("Base Layer Size",data.shape[1])
			classifier = MLPClassifier(activation='logistic',hidden_layer_sizes=(data.shape[1],2*data.shape[1],data.shape[1]/2))
		classifierList.append(technique)
		classifier.fit(dtrain,dknown)
		res = classifier.predict(dtest)
		dtknown = dtknown.tolist()
		if "prec" in metrics:
			precEllList.append(precision_score(dtknown,res, pos_label=ell))
			precSpList.append(precision_score(dtknown,res, pos_label=sp))
			print("Precision",sp, precision_score(dtknown,res, pos_label=sp))
			print("Precision",ell, precision_score(dtknown,res, pos_label=ell))
		if "acc" in metrics:
			accList.append(accuracy_score(dtknown,res))
			print("Accuracy", accuracy_score(dtknown,res))
		if "conf" in metrics:
			print("Confusion Matrix",confusion_matrix(dtknown,res))
	print("\n")

total = np.transpose(np.array([precEllList,precSpList,accList,classifierList]))
np.savetxt(savetxt,total,header="precision_Ell,precision_Sp,accuracy,technique", comments='', fmt="%s")

