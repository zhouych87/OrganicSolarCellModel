import numpy as np
import pandas as pd
import cudf, cuml
from cuml.ensemble import RandomForestRegressor as cuRF
from cuml.kernel_ridge import KernelRidge

#import matplotlib.pyplot as plt
#import pprint
from sklearn.metrics import mean_squared_error, r2_score,mean_absolute_error,make_scorer
import sys 

v_col=["z0","m","n","VLL","VHH","VLH","VHL"]
#vdata = pd.read_csv("v.out", header=None,names=v_col, delim_whitespace=True)
vdata = cudf.read_csv("v.out", header=None,names=v_col, delim_whitespace=True,dtype=np.float64)
fr=0.95

ndlist=['cip','um0.1/u3dm','dpi0.1/du3dpi','ftdm/du3dpi']
lnd=['CIP','UM','UIM-PI','FTDM']
#j=sys.argv[1]
# 0.01 348; 0.03-118; 0.1-37

method='KRR'
files=['cip','u3dm','d3dpi/du3dpi']

#j='cip'
#j='um0.1/u3dm'
#j='um1.0/u3dm'
#j='dpi0.1/du3dpi'
#j='dpi1.0/du3dpi'



label="VLL"
dim=37

for j in range(1,4):
    tnd=ndlist[j]
    nd=tnd+'.dat'
    if (tnd=='cip'):
        datainitial = cudf.read_csv(nd, header=None,delim_whitespace=True,dtype=np.float64)
        data = datainitial.iloc[:,4:] 
        alpha=0.01
        gamma=1600
    elif (tnd=='um0.1/u3dm' ):
        datainitial = cudf.read_csv(nd, header=None,delim_whitespace=True,dtype=np.float64)
        data = datainitial.iloc[:,4:34] 
        alpha=0.00004625
        gamma=0.00011875
    elif (tnd=='dpi0.1/du3dpi'):
        datainitial = cudf.read_csv(nd, header=None,delim_whitespace=True,dtype=np.float64)
        data = datainitial.iloc[:,4:34] 
        alpha=0.0004
        gamma=0.0021
    elif (tnd=='ftdm/du3dpi'):
        datainitial = cudf.read_csv(nd, header=None,delim_whitespace=True,dtype=np.float64)
        data = datainitial.iloc[:,4:dim] # skip m n
        alpha=0.0525
        gamma=0.00027
    df = cudf.concat([data,vdata.VLL.abs()],axis=1)
    df = df.to_pandas()
    df = df.astype(np.float64)
    df = df.loc[:,(df != df.iloc[0]).any()]
    df = cudf.DataFrame.from_pandas(df)
    for i in range(1,10):
        fr=0.1*i
        mmae=[]
        for rdn in range(20):
            tdf=df.copy()
            train_dataset = tdf.sample(frac=fr, random_state=rdn)
            train_labels = train_dataset.pop(label) # kick out label column from train_stats and give to train_labels
            x_train=train_dataset
            y_train=train_labels
            y_train=cudf.DataFrame(y_train)
            test_dataset = tdf.drop(train_dataset.index)
            test_labels = test_dataset.pop(label)
            x_test=test_dataset
            y_test =test_labels
            y_test=cudf.DataFrame(y_test)
            reg=KernelRidge(kernel='rbf',alpha=alpha,gamma=gamma)
            #reg=KernelRidge(kernel='laplacian',alpha=alpha,gamma=gamma)
            reg.fit(x_train,y_train)
            y_train_predicted = reg.predict(x_train)
            y_train=cudf.DataFrame(y_train).to_pandas()
            y_train_predicted=cudf.DataFrame(y_train_predicted).to_pandas()
            r2_train = r2_score(y_train, y_train_predicted)
            y_test_predicted = reg.predict(x_test)
            y_test=cudf.DataFrame(y_test).to_pandas()
            y_test_predicted=cudf.DataFrame(y_test_predicted).to_pandas()
            r2_test = r2_score(y_test, y_test_predicted)
            mae=mean_absolute_error(y_test, y_test_predicted)
            print("Train portion {:.2f} Random {} Train_R2 {:.3f} Test_R2 {:.3f} MAE {:.6f}".format(fr,rdn, r2_train, r2_test,mae))
            
            mmae.append(mae)
            if (mae>20.0):
                print("MAE is too large(>20), skip")
                break
        nmae=np.array(mmae)
        print("Descriptors {} Train_portion {:.2f}  average_MAE {:.6f}".format(lnd[j],fr,nmae.mean()))