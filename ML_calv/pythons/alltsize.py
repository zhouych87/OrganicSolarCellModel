import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pprint
from sklearn import ensemble
from sklearn import linear_model
from sklearn import svm
from sklearn.model_selection import train_test_split,GridSearchCV,cross_val_score,cross_val_predict
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.metrics import mean_squared_error, r2_score,mean_absolute_error,make_scorer
from sklearn.preprocessing import StandardScaler
#from lightgbm import LGBMRegressor
from sklearn import decomposition, datasets
from sklearn import tree
from sklearn.pipeline import Pipeline
from sklearn.kernel_ridge import KernelRidge

import sys 
import pickle

std_slc = StandardScaler()
pca = decomposition.PCA()

v_col=["m","n","VLL","VHH","VLH","VHL"]
vdata = pd.read_csv("v.out", header=None,names=v_col, delim_whitespace=True)
#fr=0.95
j=sys.argv[1]
a=float(sys.argv[2])
g=float(sys.argv[3])

nd=j+'.dat'
data = pd.read_csv(nd, header=None,delim_whitespace=True,dtype=np.float64)
dim=len(data.iloc[1])
dft = data.iloc[:,2:dim] # skip m n
dft = dft.loc[:,(dft != dft.iloc[0]).any()]
rdn=0
label="VLL"
df=pd.concat([dft,vdata.VLL.abs()],axis=1)

def dstestmain(df,fr,rdn,a,g):
    train_dataset = tdf.sample(frac=fr, random_state=rdn)  
    test_dataset = tdf.drop(train_dataset.index)
    train_labels = train_dataset.pop(label) # kick out label column from train_stats and give to train_labels
    x_train=train_dataset
    y_train=train_labels
    test_labels = test_dataset.pop(label)
    x_test=test_dataset
    y_test =test_labels
    reg=KernelRidge(kernel='rbf',alpha=a,gamma=g)
    reg.fit(x_train, y_train)
    y_train_predicted = reg.predict(x_train)
    y_test_predicted = reg.predict(x_test)
    r2_train = r2_score(y_train, y_train_predicted)
    r2_test = r2_score(y_test, y_test_predicted)
    mae=mean_absolute_error(y_test, y_test_predicted)
    return r2_train, r2_test, mae,reg

for fr in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]:
    smae=np.empty(1);sr2=np.empty(1)
    for rdn in range(20):
        tdf=df.copy()
        r2_train, r2_test, mae,reg=dstestmain(tdf,fr,rdn,a,g)
        smae=np.append(smae.copy(),mae)
        sr2=np.append(sr2.copy(),r2_test)
        print("Fr {} Data {}  Random {} Train_R2 {} Test_R2 {} \
            MAE {}".format(fr,j, rdn, r2_train, r2_test,mae))
    smae=np.delete(smae,0,axis=0)
    sr2=np.delete(sr2,0,axis=0)
    print("Fr {} Data {} Method {} average R2 {} average MAE {}".format(fr, nd, j,sr2.mean(),smae.mean()))

'''
for i in {1..128}
do
t=0.000000001*$i*5
sed "s/xxxxxxxxx/$t/" v0.py > $i.py
nohup python $i.py > $i.log &
done
'''
    
