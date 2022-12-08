# abs 
# Author: Yecheng Zhou <zhouych87@gmail.com>
#  https://zhuanlan.zhihu.com/p/554079213
#https://cloud.tencent.com/developer/article/1827841
# License: BSD 3 clause

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import ensemble
from sklearn import linear_model
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.metrics import mean_squared_error, r2_score,mean_absolute_error
from sklearn.model_selection import GridSearchCV 
#from lightgbm import LGBMRegressor
#import autosklearn
#import autosklearn.regression
#from autogluon.tabular import TabularDataset, TabularPredictor
#import autogluon.core as ag
import sys 


def slctd(path,sd):
    df=pd.DataFrame()
    sd=sd.upper()
    if ("REF" in sd ):
        fn=path+"/refpen.dat" 
        refdata = pd.read_csv(fn, header=0, delim_whitespace=True)
        reffeature=['xshort','ynoraml','zshort','anglelong','angleshort','anglenormal']
        df[reffeature]=refdata[reffeature].copy() 
    if ("PI" in sd ):
        fn=path+"/pigeomty.dat"
        pidata = pd.read_csv(fn, header=None,  delim_whitespace=True)
        dim=len(pidata.iloc[1])
        pidatam=pidata.iloc[:,dim-9:dim]
        pifeature=['Dcen','Datom', 'Dcenp1', 'Dcenp2', 'CosL', 'CosS', 'CosN', 'CosO', 'AreaPJ']
        pidatam.columns=pifeature
        df[pifeature]=pidatam.copy()
    if ("MOL" in sd ):
        fn=path+"/mgeomty.dat" 
        mdata = pd.read_csv(fn, header=0, delim_whitespace=True)
        mfeature=['Dcen','Datom', 'Dcenp1', 'Dcenp2', 'CosL', 'CosS', 'CosN', 'CosO', 'AreaPJ']
        df[mfeature]=mdata[mfeature].copy() 
    if ("U3D" in sd ):
        fn=path+"/u3dmrs.dat" # u,c,p,m,v,e #u0.24,c0.178,p0.23,m0.27,v0.30,e0.20
        ums3dfeature=[]
        for i in range(32):
            tmp='S'+str(i)
            ums3dfeature.append(tmp)
        ums3ddata = pd.read_csv(fn, header=0, delim_whitespace=True)
        df[ums3dfeature]=ums3ddata[ums3dfeature].copy()   
    if ("C3D" in sd ):
        fn=path+"/c3dmrs.dat" # u,c,p,m,v,e #u0.24,c0.178,p0.23,m0.27,v0.30,e0.20
        cms3dfeature=[]
        for i in range(32):
            tmp='S'+str(i)
            cms3dfeature.append(tmp)
        cms3ddata = pd.read_csv(fn, header=0, delim_whitespace=True)
        df[cms3dfeature]=cms3ddata[cms3dfeature].copy()
    if ("P3D" in sd ):
        fn=path+"/p3dmrs.dat" # u,c,p,m,v,e #u0.24,c0.178,p0.23,m0.27,v0.30,e0.20
        pms3dfeature=[]
        for i in range(32):
            tmp='S'+str(i)
            pms3dfeature.append(tmp)
        pms3ddata = pd.read_csv(fn, header=0, delim_whitespace=True)
        df[pms3dfeature]=pms3ddata[pms3dfeature].copy()
    if ("M3D" in sd ):
        fn=path+"/m3dmrs.dat" # u,c,p,m,v,e #u0.24,c0.178,p0.23,m0.27,v0.30,e0.20
        mms3dfeature=[]
        for i in range(32):
            tmp='S'+str(i)
            mms3dfeature.append(tmp)
        mms3ddata = pd.read_csv(fn, header=0, delim_whitespace=True)
        df[mms3dfeature]=mms3ddata[mms3dfeature].copy() # u,c,p,m,v,e 
    if ("V3D" in sd ):
        fn=path+"/v3dmrs.dat" # u,c,p,m,v,e #u0.24,c0.178,p0.23,m0.27,v0.30,e0.20
        vms3dfeature=[]
        for i in range(32):
            tmp='S'+str(i)
            vms3dfeature.append(tmp)
        vms3ddata = pd.read_csv(fn, header=0, delim_whitespace=True)
        df[vms3dfeature]=vms3ddata[vms3dfeature].copy()
    if ("E3D" in sd ):
        fn=path+"/e3dmrs.dat" # u,c,p,m,v,e #u0.24,c0.178,p0.23,m0.27,v0.30,e0.20
        ems3dfeature=[]
        for i in range(32):
            tmp='S'+str(i)
            ems3dfeature.append(tmp)
        ems3ddata = pd.read_csv(fn, header=0, delim_whitespace=True)
        df[ems3dfeature]=ems3ddata[ems3dfeature].copy()
    df= df.loc[:,~df.T.duplicated(keep='first')] # remove repeat columns ,float become object
    v_col=["m","n","VLL","VHH","VLH","VHL"]
    fn=path+"/v.out"
    vdata = pd.read_csv(fn, header=None,names=v_col, delim_whitespace=True)
    df["v"] = np.abs(vdata["VLL"])
    return df

def molmethod(mm,rdn):
    reg=ensemble.RandomForestRegressor(random_state=rdn)
    mm=mm.upper() 
    if (mm == "LGBM") :
        reg = LGBMRegressor(random_state=rdn)
    if (mm == "RF"):
        reg = ensemble.RandomForestRegressor(random_state=rdn)
    if (mm == "LR") :
        reg = linear_model.LinearRegression()
    if (mm == "DECISIONTREE") :
        reg = DecisionTreeRegressor(random_state=rdn)
    if (mm == "SGD") :
        reg = linear_model.SGDRegressor(random_state=rdn,max_iter=10000)
    if (mm == "LASSO") :
        reg =linear_model.Lasso(random_state=rdn,max_iter=10000)
    if (mm == "ENETCV") :
        reg = linear_model.ElasticNetCV(random_state=rdn,max_iter=10000)
    if (mm == "ENET") :
        reg = linear_model.ElasticNet(random_state=rdn,max_iter=10000)
    if (mm == "RIDGE") :
        reg = linear_model.Ridge(random_state=rdn,max_iter=20000)
    if (mm == "SVR") :
        reg=  svm.SVR()
    if (mm == "NUSVR") :
        reg = svm.NuSVR()
    if (mm == "KNN") :
        reg = KNeighborsRegressor()
    if (mm == "GB") :
        reg = ensemble.GradientBoostingRegressor(random_state=rdn)
    if (mm == "AUTOSK") :
        reg = autosklearn.regression.AutoSklearnRegressor(
        time_left_for_this_task=300,
        memory_limit=None,
        resampling_strategy='cv',resampling_strategy_arguments={'folds': 5})
    return reg 

def dstestmain(nd,method,fr,rdn):
    df=slctd(".",nd)
    label="v"
    #print(df)
    train_dataset = df.sample(frac=fr, random_state=rdn)  
    test_dataset = df.drop(train_dataset.index)
    train_labels = train_dataset.pop(label) # kick out label column from train_stats and give to train_labels
    x_train=train_dataset
    y_train=train_labels
    test_labels = test_dataset.pop(label)
    x_test=test_dataset
    y_test =test_labels
    reg = molmethod(method,rdn)
    reg.fit(x_train, y_train)
    y_train_predicted = reg.predict(x_train)
    r2_train = r2_score(y_train, y_train_predicted)
    #print("Data: {} Method: {} Train R2: {} ".format(nd, method, r2_train))
    y_test_predicted = reg.predict(x_test)
    r2_test = r2_score(y_test, y_test_predicted)
    mae=mean_absolute_error(y_test, y_test_predicted)
    print("Data: {} Method: {} Random: {} Train_R2: {} Test_R2: {} MAE: {}".format(nd, method,rdn, r2_train, r2_test,mae))
    #return r2_train, r2_test, mae
    return y_train, y_train_predicted,y_test, y_test_predicted,r2_train, r2_test, mae


def t3plot(i,j,fr):
    import numpy as np
    import matplotlib.pyplot as plt
    import t3 
    
    s = 50
    a = 1.0
    b = 0.8 
    x = np.linspace(0,120,240)
    y=x
    
    plt.close()
    plt.ion()
    plt.figure(figsize=(12.8,9.6))
    
    for rdn in range(1,10):
        y_train, y_train_predicted,y_test, y_test_predicted,r2_train, r2_test, mae=t3.dstestmain(i,j,fr,rdn)
        tmp=330+rdn
        plt.subplot(tmp)
        plt.plot(x,y, 'black')
        plt.scatter( y_train, y_train_predicted, 
                    c="c", s=50, marker="+", alpha=b)
        
        plt.scatter( y_test, y_test_predicted,edgecolor='r',
                    c="orange", s=s, marker="^", alpha=a,
                    label="MAE =%.3f" % mae)
        plt.legend()
    
    figname='./fig/'+j+'_'+i+'_'+str(int(fr*100))+'_py1.png'
    plt.savefig(figname,format='png')



#cfr=[0.95,0.90,0.8,0.7,0.6,0.5,0.4]
if __name__ == '__main__':
    #mth=["rf","lr","decisiontree","sgd","lasso","enetcv","enet","ridge","knn","gb","svr","nusvr","lgbm"]
    #named=["pi","mol","u3d","c3d","p3d","m3d","v3d","e3d"]
    #print("sys 1,",sys.argv[1])
    #print("sys 2,",sys.argv[2])
    mth=['knn','gb','decisiontree','rf']
    named=["ref","mol","u3d","e3d"]
    named=["pi"]
    fr=float(sys.argv[1])
    for i in named: 
        for j in mth:
            for rdn in range(100):
                dstestmain(i,j,fr,rdn)
    '''
    mth=["rf","lr","decisiontree","sgd","lasso","enetcv","enet","ridge","knn","gb","svr","nusvr","lgbm"]
    named=["pi","mol","u3d","c3d","p3d","m3d","v3d","e3d"]
    
    mth=['knn','gb','decisiontree','rf']
    named=["ref","mol","u3d","e3d"]
    tfr=[0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]
    
    for j in mth:
        for i in named:
            for fr in tfr:
                t3plot(i,j,fr)
'''




