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

def train(f,j,a,g):
    nd=f+'/'+j+'.dat'
    data = pd.read_csv(nd, header=None,delim_whitespace=True,dtype=np.float64)
    dim=len(data.iloc[1])
    dft = data.iloc[:,2:dim] # skip m n
    dft = dft.loc[:,(dft != dft.iloc[0]).any()]
    v_col=["m","n","VLL","VHH","VLH","VHL"]
    vdata = pd.read_csv(f+'/''v.out', header=None,names=v_col, delim_whitespace=True)
    rdn=0
    label="VLL"
    df=pd.concat([dft,vdata.VLL.abs()],axis=1)
    
    train_dataset = df.sample(frac=1.0,random_state=0)  
    train_labels = train_dataset.pop(label) # kick out label column from train_stats and give to train_labels
    x_train=train_dataset
    y_train=train_labels
    reg=KernelRidge(kernel='rbf',alpha=a,gamma=g)
    reg.fit(x_train, y_train)
    y_train_predicted = reg.predict(x_train)
    r2 = r2_score(y_train, y_train_predicted)
    mae=mean_absolute_error(y_train, y_train_predicted)
    return r2, mae,reg

def test(f,j,reg):
    nd=f+'/'+j+'.dat'
    data = pd.read_csv(nd, header=None,delim_whitespace=True,dtype=np.float64)
    dim=len(data.iloc[1])
    dft = data.iloc[:,2:dim] # skip m n
    dft = dft.loc[:,(dft != dft.iloc[0]).any()]
    v_col=["m","n","VLL","VHH","VLH","VHL"]
    vdata = pd.read_csv(f+'/''v.out', header=None,names=v_col, delim_whitespace=True)
    rdn=0
    label="VLL"
    df=pd.concat([dft,vdata.VLL.abs()],axis=1)
    test_dataset = df.sample(frac=1.0,random_state=0)  
    test_labels = test_dataset.pop(label)
    x_test=test_dataset
    y_test=test_labels
    y_test_predicted = reg.predict(x_test)
    r2= r2_score(y_test, y_test_predicted)
    mae=mean_absolute_error(y_test, y_test_predicted)
    return r2, mae,reg


f=sys.argv[1]
j=sys.argv[2]
a=float(sys.argv[3])
g=float(sys.argv[4])
r2, mae,reg=train(f,j,a,g)
print("Train_set {} Descriptor {}".format(f,j))

for sets in ['00','1','2','3','4','5','6','7','8','9','10','50ns','film']:
    r2, mae,reg=test(sets,j,reg)
    print("Test_set {} R2 {}  MAE {}".format(j,r2,mae))

'''
for i in {1..128}
do
t=0.000000001*$i*5
sed "s/xxxxxxxxx/$t/" v0.py > $i.py
nohup python $i.py > $i.log &
done
'''
    
