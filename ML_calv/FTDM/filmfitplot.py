import numpy as np
import pandas as pd
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import mean_squared_error, r2_score,mean_absolute_error,make_scorer
import sys 
import matplotlib.pyplot as plt

cuda=1

files=['00','1','2','3','4','5','6','7','8','9','10','50ns','film']
files=['film','tet']
label="VLL"
dim=37
#alpha=0.0509
#gamma=1.12e-4

#4.0 2.5
#0.125893	0.316228

alpha=0.125893
gamma=0.316228

if (cuda==1):
    import cudf, cuml
    from cuml.kernel_ridge import KernelRidge
    v_col=["z0","m","n","VLL","VHH","VLH","VHL"]
else:
    from sklearn.kernel_ridge import KernelRidge
    v_col=["m","n","VLL","VHH","VLH","VHL"]



def gpudata(j):
    nd=j+'/du3dpi.dat'
    dim=37
    vdata = cudf.read_csv(j+"/v.out", header=None,names=v_col, delim_whitespace=True,dtype=np.float64)
    datainitial = cudf.read_csv(nd, header=None,delim_whitespace=True,dtype=np.float64)
    data = datainitial.iloc[:,4:dim]/datainitial.iloc[1,3] # skip m n
    df=cudf.concat([data,vdata.VLL.abs()],axis=1)
    data=None
    vdata=None
    df=df.to_pandas()
    df=df.astype(np.float64)
    df = df.loc[:,(df != df.iloc[0]).any()]
    df = cudf.DataFrame.from_pandas(df)
    return df

def cpudata(j):
    nd=j+'/du3dpi.dat'
    vdata = pd.read_csv(j+"/v.out", header=None,names=v_col, delim_whitespace=True,dtype=np.float64)
    datainitial = pd.read_csv(nd, header=None,delim_whitespace=True,dtype=np.float64)
    data = datainitial.iloc[:,3:dim]/datainitial.iloc[1,2] # skip m n
    df=pd.concat([data,vdata.VLL.abs()],axis=1)
    df=df.astype(np.float64)
    df = df.loc[:,(df != df.iloc[0]).any()]
    return df

plts = 50
plta = 1.0
pltb = 0.8 
x = np.linspace(0,260,260)
y=x
plt.close()
plt.ion()
plt.figure(figsize=(16,7),constrained_layout=True)
plt.rcParams['font.size'] = '18'
tmp=120
for j in files:
    if (cuda==1):
        df=gpudata(j)
    else:
        df=cpudata(j)
    
    train_dataset = df.copy()
    train_labels = train_dataset.pop(label) # kick out label column from train_stats and give to train_labels
    x_train=train_dataset
    y_train=train_labels
    if (cuda==1):
        y_train=cudf.DataFrame(y_train)
    
    reg=KernelRidge(kernel='rbf',alpha=alpha,gamma=gamma)
    reg.fit(x_train, y_train)
    tmp=tmp+1
    plt.subplot(tmp)
    plt.plot(x,y, 'black')
    test_dataset = df.copy()
    test_labels = test_dataset.pop(label)
    x_test=test_dataset
    y_test =test_labels
    y_test_predicted = reg.predict(x_test)
    if (cuda==1):
        y_test=cudf.DataFrame(y_test).to_pandas()
        y_test_predicted=cudf.DataFrame(y_test_predicted).to_pandas()
    
    r2_test = r2_score(y_test, y_test_predicted)
    mae=mean_absolute_error(y_test, y_test_predicted)
    print("Trained on {},test set {:} Test_R2 {:.3f} MAE {:.6f}".format(j,j, r2_test,mae))
    plt.scatter(y_test, np.abs(y_test_predicted),c="c", s=50, marker=".", alpha=pltb)
    for i in files:
        if (j!=i):
            if (cuda==1):
                df=gpudata(i)
            else:
                df=cpudata(i)
            
            test_dataset = df.copy()
            test_labels = test_dataset.pop(label)
            x_test=test_dataset
            y_test =test_labels
            y_test_predicted = reg.predict(x_test)
            if (cuda==1):
                y_test=cudf.DataFrame(y_test).to_pandas()
                y_test_predicted=cudf.DataFrame(y_test_predicted).to_pandas()
            
            r2_test = r2_score(y_test, y_test_predicted)
            mae=mean_absolute_error(y_test, y_test_predicted)
            print("Trained on {},test set {:} Test_R2 {:.3f} MAE {:.6f}".format(j,i, r2_test,mae))
            plt.scatter(y_test, np.abs(y_test_predicted),
            c="red", s=plts, marker="+", alpha=plta,label="MAE =%.3f" % mae)
    plt.xlabel('DFT calculated EC (meV)')
    plt.ylabel('ML predicted EC (meV)')
    plt.grid(True)

figname='transferfilmv3.png'
plt.savefig(figname,format='png')
