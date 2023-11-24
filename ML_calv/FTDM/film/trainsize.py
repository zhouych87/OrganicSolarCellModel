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
vdata = cudf.read_csv("v.out", header=None,names=v_col, delim_whitespace=True,dtype=np.float32)
fr=0.95


#j=sys.argv[1]
# 0.01 348; 0.03-118; 0.1-37


#j='0.03m';dim=118
#j='0.01m';dim=348
nd='du3dpi.dat'

#data = pd.read_csv(nd, header=None,delim_whitespace=True,dtype=np.float64)
datainitial = cudf.read_csv(nd, header=None,delim_whitespace=True,dtype=np.float32)
#dim=len(data.iloc[1])
#dim=37 # 0.01 348; 0.03-118; 0.1-37
label="VLL"
dim=37

data = datainitial.iloc[:,4:dim] # skip m n
df=cudf.concat([data,vdata.VLL.abs()],axis=1)
data=None
df=df.to_pandas()
df=df.astype(np.float32)
df = df.loc[:,(df != df.iloc[0]).any()]
df = cudf.DataFrame.from_pandas(df)

alpha=0.0509
gamma=1.12e-4

for i in range(1,10):
    fr=0.1*i
    mmae=[]
    for rdn in range(50):
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
    print("Train portion {:.2f}  average MAE {:.6f}".format(fr,nmae.mean()))


'''
for rdn in range(10):
    label="VLL"
    tdf=df.copy()
    train_dataset = tdf.sample(frac=fr, random_state=rdn)  
    test_dataset = tdf.drop(train_dataset.index)
    train_labels = train_dataset.pop(label) # kick out label column from train_stats and give to train_labels
    x_train=train_dataset
    y_train=train_labels
    y_train=cudf.DataFrame(y_train)
    test_labels = test_dataset.pop(label)
    x_test=test_dataset
    y_test =test_labels
    y_test=cudf.DataFrame(y_test)
    tdf=None;train_dataset=None;test_dataset=None
    for alpha in alphac:
        for gamma in gc:
            reg=KernelRidge(kernel='rbf',alpha=alpha,gamma=gamma)
            reg.fit(x_train, y_train)
            y_train_predicted = reg.predict(x_train)
            y_train=cudf.DataFrame(y_train).to_pandas()
            y_train_predicted=cudf.DataFrame(y_train_predicted).to_pandas()
            r2_train = r2_score(y_train, y_train_predicted)
            y_test_predicted = reg.predict(x_test)
            y_test=cudf.DataFrame(y_test).to_pandas()
            y_test_predicted=cudf.DataFrame(y_test_predicted).to_pandas()
            r2_test = r2_score(y_test, y_test_predicted)
            mae=mean_absolute_error(y_test, y_test_predicted)
            print("alpha {} gamma {} Random {} Train_R2 {} Test_R2 {} MAE {}".format(alpha,gamma,rdn, r2_train, r2_test,mae))
            if (mae>20.0):
                print("MAE is too large(>20), skip")
                break
'''
#sets=np.delete(sets,0,axis=0)
#print("Fr {} Data {} Method {} average MAE {}".format(fr, nd, j,sets.mean()))

'''
for i in {1..128}
do
t=0.000000001*$i*5
sed "s/xxxxxxxxx/$t/" v0.py > $i.py
nohup python $i.py > $i.log &
done
'''

