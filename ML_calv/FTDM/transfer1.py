import numpy as np
import pandas as pd
import cudf, cuml
from cuml.ensemble import RandomForestRegressor as cuRF
from cuml.kernel_ridge import KernelRidge

#import matplotlib.pyplot as plt
#import pprint
from sklearn.metrics import mean_squared_error, r2_score,mean_absolute_error,make_scorer
import sys 


files=['00','1','2','3','4','5','6','7','8','9','10','50ns','film']
#files=['film','tet']
label="VLL"
dim=37
alpha=0.0509
gamma=1.12e-4

#after hyperparameter optimiztion, but not good
#alpha=0.06138
#gamma=0.0001126

for j in files:
    nd=j+'/du3dpi.dat'
    v_col=["z0","m","n","VLL","VHH","VLH","VHL"]
    vdata = cudf.read_csv(j+"/v.out", header=None,names=v_col, delim_whitespace=True,dtype=np.float32)
    datainitial = cudf.read_csv(nd, header=None,delim_whitespace=True,dtype=np.float32)
    data = datainitial.iloc[:,3:dim] # skip m n
    df=cudf.concat([data,vdata.VLL.abs()],axis=1)
    data=None
    df=df.to_pandas()
    df=df.astype(np.float32)
    df = df.loc[:,(df != df.iloc[0]).any()]
    df = cudf.DataFrame.from_pandas(df)
    train_dataset = df.copy()
    train_labels = train_dataset.pop(label) # kick out label column from train_stats and give to train_labels
    x_train=train_dataset
    y_train=train_labels
    y_train=cudf.DataFrame(y_train)
    reg=KernelRidge(kernel='rbf',alpha=alpha,gamma=gamma)
    reg.fit(x_train, y_train)
    print("{} trained".format(j))
    for i in files:
        vdata = cudf.read_csv(i+"/v.out", header=None,names=v_col, delim_whitespace=True,dtype=np.float32)
        datainitial = cudf.read_csv(i+'/du3dpi.dat', header=None,delim_whitespace=True,dtype=np.float32)
        data = datainitial.iloc[:,4:dim] # skip m n
        df=cudf.concat([data,vdata.VLL.abs()],axis=1)
        data=None
        df=df.to_pandas()
        df=df.astype(np.float32)
        df = df.loc[:,(df != df.iloc[0]).any()]
        df = cudf.DataFrame.from_pandas(df)
        test_dataset = df.copy()
        test_labels = test_dataset.pop(label)
        x_test=test_dataset
        y_test =test_labels
        y_test=cudf.DataFrame(y_test)
        y_test_predicted = reg.predict(x_test)
        y_test=cudf.DataFrame(y_test).to_pandas()
        y_test_predicted=cudf.DataFrame(y_test_predicted).to_pandas()
        r2_test = r2_score(y_test, y_test_predicted)
        mae=mean_absolute_error(y_test, y_test_predicted)
        print("Trained on {},test set {:} Test_R2 {:.3f} MAE {:.6f}".format(j,i, r2_test,mae))
