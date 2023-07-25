# I cannot write the data as txt, csv, etc. because of memory size.


# Necessary packages

import time
import random
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score

###################### Data preparation ####################################


# L: lncRNA (240*1)
L = pd.read_excel("./01-lncRNAs-240.xlsx", sheet_name=0, header=None)
# D: diseases (412*1)
D = pd.read_excel("./02-diseases-412.xlsx", sheet_name=0, header=None)
# M: miRNA (495*1)
M = pd.read_excel("./03-miRNAs-495.xlsx", sheet_name=0, header=None)
# LL: lncRNA-lncRNA functional similarities (240*240)
LL = pd.read_excel("./04-lncRNA-lncRNA.xlsx", sheet_name=0, header=None)
# LD: lncRNA-disease associations (240*412)
LD = pd.read_excel("./05-lncRNA-disease.xlsx", sheet_name=0, header=None)
# MD: miRNA-disease associations (495*412)
MD = pd.read_excel("./06-miRNA-disease.xlsx", sheet_name=0, header=None)
# DD: disease-disease semantic similarities (412*412)
DD = pd.read_excel("./07-disease-disease.xlsx", sheet_name=0, header=None)
# LM: lncRNA-miRNA interactions (240*495)
LM = pd.read_excel("./08-lncRNA-miRNA.xlsx", sheet_name=0, header=None)




# Constructing sample dataset

# Represent lncRNA by 1147 features, L1147: (240*1147)
L1147 = pd.concat([LL, LM, LD], axis=1)
# Adding lncRNA name column for L1147
L1148 = pd.concat([L.iloc[:, 0], L1147], axis=1)

# Represent disease by 1147 features, D1147: (412*1147)
# DL: the transfer matrix of LD
DL = LD.transpose()
# DM: the transfer matrix of MD
DM = MD.transpose()
# D1147: (412*1147)
D1147 = pd.concat([DL, DM, DD], axis=1)
# Adding disease name column for D1147
D1148 = pd.concat([D.iloc[:, 0], D1147], axis=1)

# Merge L1148 and D1148 to LDALL (98880*(2+1147+1147=2296))
LDALL = L1148.merge(D1148, how='cross')





### Adjust column position of LDALL (98880*(2+1147+1147=2296)) 
LDALL.columns = [i for i in range(1,2297)]
d1 = LDALL.iloc[:, 0]
d2 = LDALL.iloc[:, 1148]
d3 = LDALL.drop(columns=[LDALL.columns[0], LDALL.columns[1148]])

LDALL = pd.concat([d1, d2, d3], axis=1)
LDALL.to_csv("./lncRNA-disease-ALL.csv", index=False, header=False)

# Read the adjusted dataframe from the Excel file
LDALL = pd.read_csv("./lncRNA-disease-ALL.csv", header=None)





### exclude columns with full zero values,(98880*(3+1952=1955)) 
# Select columns 1 and 2
d1 = LDALL.iloc[:, [0, 1]]

# Exclude columns 1 and 2
d2 = LDALL.iloc[:, 2:]
d2 = d2.loc[:, d2.sum() > 0]

# Combine d1 and d2
LDExcl0 = pd.concat([d1, d2], axis=1)

# Write to Excel
LDExcl0.to_csv("./lncRNA-disease-Excl0.csv", index=False, header=False)

# Read from Excel
LDExcl0 = pd.read_csv("./lncRNA-disease-Excl0.csv", header=None)





# Construct known lncRNA-disease association pairs
xy = pd.DataFrame(np.where(LD == 1))
xy = xy.transpose()
xy.columns = ["rows", "cols"]
xy = xy.sort_values(["cols", "rows"])
cL = L.iloc[xy.iloc[:,0],0].reset_index().drop("index", axis=1)
cD = D.iloc[xy.iloc[:,1],0].reset_index().drop("index", axis=1)
LDA = pd.DataFrame({'lncRNA': cL, 'disease': cD})
LDA = pd.DataFrame({'lncRNA': cL.iloc[:,0], 'disease': cD.iloc[:,0]})
LDA.to_csv("./lncRNA-disease-associations-2697.csv", header=False, index=False)

# Read the lncRNA-disease associations from the Excel file
LDA = pd.read_csv("./lncRNA-disease-associations-2697.csv", header=None)





### set label for each sample, 1 for known LDA, 0 for unknown LDA (98880*1955)
LDExcl0.insert(0, 'label', 0)
for i in range(len(LDExcl0)):
    for j in range(len(LDA)):
        if LDExcl0.iloc[i,1] == LDA.iloc[j, 0] and LDExcl0.iloc[i,2] == LDA.iloc[j,1]:
            LDExcl0.iloc[i,0] = 1




### regularize samples using maximum-minimum regularization method (98880*1955)
B1 = LDExcl0.iloc[:, 3:]
scaler = MinMaxScaler()
B1s = scaler.fit_transform(B1)
B1s = pd.DataFrame(B1s, columns=B1.columns)
sLDExcl0 = pd.concat([LDExcl0.iloc[:,0:3], B1s], axis=1)
sLDExcl0.to_csv("./lncRNA-disease-Excl0-regulation-scaled.csv", header=False, index=False)
sLDExcl0 = pd.read_csv("./lncRNA-disease-Excl0-regulation-scaled.csv", header=None)





### consturct positive samples (2697*1955)
PositiveSample = sLDExcl0[sLDExcl0["label"]==1]
PositiveSample.to_csv("./PositiveSample.csv", header=False, index=False)
PositiveSample = pd.read_csv("./PositiveSample.csv", header=None)




### consturct unlabeled samples ((98880-2697=96183)*1955)
UnlabeledSample = sLDExcl0[sLDExcl0["label"]==0]
UnlabeledSample.to_csv("./UnlabeledSample.csv", header=False, index=False)
UnlabeledSample = pd.read_csv("./UnlabeledSample.csv", header=None)





### consturct negative samples (2697*1955)
random.seed(1234)
sp = random.sample(range(UnlabeledSample.shape[0]), PositiveSample.shape[0])
sps = sorted(sp)
NegativeSample = UnlabeledSample.iloc[sps,:].copy()
NegativeSample.to_csv("./NegativeSample.csv", header=False, index=False)
NegativeSample = pd.read_csv("./NegativeSample.csv", header=None)





### construct training sample set by combining positive and negative samples (5394*1955)
TrainingSample = pd.concat([PositiveSample, NegativeSample], axis=0)
TrainingSample.to_csv("./TrainingSample.csv", index=False, header=False)
TrainingSample = pd.read_csv("./TrainingSample.csv", header=None)
pd.DataFrame(TrainingSample.columns).to_csv("./TrainingSample_colnames.csv")





########################### compute variable importance score using RandonForest #############################
B1 = pd.concat([TrainingSample.iloc[:,0], TrainingSample.iloc[:,3:]], axis=1)

### compute variable importance score

def compute_feature_importance(data):
    # Split the data into features (X) and the target variable (y)
    X = data.iloc[:, 1:]
    y = data.iloc[:, 0]

    # Initialize a Random Forest Classifier with mtry=651 and other parameters
    # 651 is the one-third of the total number of features. (1952/3)
    rf = RandomForestClassifier(n_estimators=500, max_features=651, random_state=123)

    # Train the Random Forest model
    rf.fit(X, y)

    # Get variable importance scores
    im = rf.feature_importances_

    # Repeat 9 times
    for i in range(2, 11):
        print(i)
        # Train the Random Forest model again
        rf.fit(X, y)

        # Accumulate variable importance scores
        im += rf.feature_importances_

    # Compute average variable importance score
    im = im / 10

    # Create a DataFrame with feature names and corresponding variable importance scores
    fs = pd.DataFrame({'Feature': X.columns, 'Importance': im})

    # Sort features by their variable importance scores in descending order
    fs = fs.sort_values(by='Importance', ascending=False)

    # Write the DataFrame to an Excel file
    fs.to_csv("./FeatureScore.csv", index=False, header=False)

# Assuming 'B1' is a DataFrame in Python
data = B1.copy()

# Measure the time taken to compute feature importance
start_time = time.time()
compute_feature_importance(data)
end_time = time.time()
execution_time = end_time - start_time
print("Execution Time: ", execution_time, " seconds")




############# training RandomForest model with top 50, 100, ..., 1950 features ###############################
# Read the training sample set consisting of 5394 lncRNA-disease pairs
B = pd.read_csv("./TrainingSample.csv", header=None)
# Read the variable importance score of each feature
fs = pd.read_csv("./FeatureScore.csv", header=None)

def compute_classification_accuracy(B):
    # Initialize matrix to record classification accuracy
    gaccuracy = np.zeros((39, 2))

    tt = 50
    while tt < fs.shape[0]:
        print("Number of features:", tt)
        # Initialize array to store accuracy in each fold of cross-validation
        laccuracy = np.zeros(10)
        # Initialize mean accuracy in cross-validation
        lmeanaccuracy = 0

        # Construct training sample subset consisting of X1(label) + top tt features
        ttt = fs.iloc[:tt, 0].tolist()
        B1 = B[ttt]
        B4 = pd.concat([B['label'], B1], axis=1)

        # Random resampling in 10-fold cross-validation
        ind = np.random.choice(10, size=B4.shape[0], replace=True, p=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
        print("ind:", ind)
        
        for i in range(10):
            # Train Random Forest model
            rf = RandomForestClassifier(n_estimators=500, max_features=int(tt/3))
            rf.fit(B1[ind != i], B4['label'][ind != i])

            # Predict using Random Forest model
            pred = rf.predict(B1[ind == i])

            # Threshold predictions at 0.5
            pred = np.where(pred >= 0.5, 1, 0)
            # print("pred:", pred)

            # Compute classification accuracy
            laccuracy[i] = accuracy_score(B4['label'][ind == i], pred)
            print("laccuracy:", laccuracy)
            lmeanaccuracy += laccuracy[i] / 10
            

        gaccuracy[int(tt/50)-1, 0] = tt
        gaccuracy[int(tt/50)-1, 1] = lmeanaccuracy
        print("lmeanaccuracy:", lmeanaccuracy)
        print("gaccuracy:\n", gaccuracy)

        tt += 50
        

    gaccuracy = pd.DataFrame(gaccuracy)
    gaccuracy.to_csv("./TrainingSample-gaccuracy.csv", index=False, header=False)

start_time = time.time()
compute_classification_accuracy(TrainingSample)
end_time = time.time()
execution_time = end_time - start_time
print("Execution Time: ", execution_time, " seconds")





########################## 5-fold crossing validation on training sample set with 300 features ###############
### read variable importance score of each feature
fs = pd.read_csv("./FeatureScore.csv", header=None)

# read training sample set consisting of 5394 lncRNA-disease pairs (5394*(3+1952))
B = pd.read_csv("./TrainingSample.csv", header=None)
# extract subset consisting of top 300 featues
tt = 300
ttt = fs.iloc[0:tt, 0]
B1 = B[ttt]
B2 = B.iloc[:, 0]
# TB is training sample set without column X2（lncRNA name）and X3（disease name）
TB = pd.concat([B2, B1], axis=1)

# Read unlabeled sample set consisting of 96183 lncRNA-disease pairs
BB = pd.read_csv("./UnlabeledSample.csv", header=None)
# Extract subset consisting of top 300 features
tt = 300
ttt = fs.iloc[0:tt, 0]
B1 = BB[ttt]
B2 = BB.iloc[:, 0]
# NB is unlabeled sample set without column X2 (lncRNA name) and X3 (disease name)
NB = pd.concat([B2, B1], axis=1)

def FiveFoldCrossingValidation(TB, NB):
    sumauc = 0
    sumap = 0
    PTB = TB.iloc[0:2696,:]
    NTB = TB.iloc[2697:5393,:]
    
    PTB.index = list(range(0,2696))
    NTB.index = list(range(0,2696))

    for i in range(1, 5):
        # Training sample set
        PTB1 = PTB.drop(list(range((540 * (i - 1)), 540 * i + 1)))
        NTB1 = NTB.drop(list(range((540 * (i - 1)), 540 * i + 1)))
        TrainB = pd.concat([PTB1, NTB1])
        
        # Test sample set
        PTB2 = PTB.iloc[(540 * (i - 1)):(540 * i) + 1]
        TestB = pd.concat([PTB2, NB])
        
        # Train RandomForest Model with parameters
        rf = RandomForestClassifier(n_estimators=500, max_features=100)
        rf.fit(TrainB.iloc[:,1:], TrainB.iloc[:,0])
        
        # Predict using RandomForest Model
        pred = rf.predict(TestB.iloc[:,1:])
        
        # Compute AUC value
        auc = roc_auc_score(TestB.iloc[:,0], pred)
        print(auc)
        sumauc += auc

        # Compute AUPR value
        ap = average_precision_score(TestB.iloc[:,0], pred)
        print(ap)
        sumap += ap

    # Last fold (i = 5)
    i = 5
    PTB1 = PTB.drop(list(range((540 * (i - 1)), 2696)))
    NTB1 = NTB.drop(list(range((540 * (i - 1)), 2696)))
    TrainB = pd.concat([PTB1, NTB1])
    
    PTB2 = PTB.iloc[(540 * (i - 1)):2698]
    TestB = pd.concat([PTB2, NB])
    
    rf = RandomForestClassifier(n_estimators=500, max_features=100)
    rf.fit(TrainB.iloc[:,1:], TrainB.iloc[:,0])
    
    pred = rf.predict(TestB.iloc[:,1:])
    
    auc = roc_auc_score(TestB.iloc[:,0], pred)
    print(auc)
    sumauc += auc
    sumauc /= 5
    print(sumauc)
    
    ap = average_precision_score(TestB.iloc[:,0], pred)
    print(ap)
    sumap += ap
    sumap /= 5
    print(sumap)
    
# Measure execution time
start_time = time.time()
FiveFoldCrossingValidation(TB, NB)
print("Execution time:", time.time() - start_time)




################################### predict all lncRNA-disease samples with 300 fearues#######################
### read variable importance score of each feature
fs = pd.read_csv("./FeatureScore.csv", header=None)

# read training sample set consisting of 5394 lncRNA-disease pairs (5394*(3+1952))
B = pd.read_csv("./TrainingSample.csv", header=None)
# extract subset consisting of top 300 featues
tt = 300
ttt = fs.iloc[0:tt, 0]
B1 = B[ttt]
B2 = B.iloc[:, 0]
# TB is training sample set without column X2（lncRNA name）and X3（disease name）
TB = pd.concat([B2, B1], axis=1)

# Read unlabeled sample set consisting of 96183 lncRNA-disease pairs
BB = pd.read_csv("./UnlabeledSample.csv", header=None)
# Extract subset consisting of top 300 features
tt = 300
ttt = fs.iloc[0:tt, 0]
B1 = BB[ttt]
B2 = BB.iloc[:, 0]
# NB is unlabeled sample set without column X2 (lncRNA name) and X3 (disease name)
NB = pd.concat([B2, B1], axis=1)


### train RandomForest Model with parameter，try=the number of features（300）/3
rf = RandomForestClassifier(n_estimators=500, max_features=100)
rf.fit(TB.iloc[:,1:], TB.iloc[:,0])

### predict all unknwon samples using RandomForest Model
pred = rf.predict(NB.iloc[:,1:])
NB1 = BB.iloc[:,0:3]
UnlabeledSampleScore = pd.concat([NB1, pd.DataFrame(pred)], axis = 1)
UnlabeledSampleScore.to_csv("./UnlabeledSampleScore-300-features.csv", index=False, header=False)
