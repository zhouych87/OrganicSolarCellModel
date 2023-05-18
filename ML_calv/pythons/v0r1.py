import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pprint
from sklearn import ensemble
from sklearn import linear_model
from sklearn import svm
from sklearn.model_selection import train_test_split,GridSearchCV 
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

#import autosklearn
#import autosklearn.regression
#from autogluon.tabular import TabularDataset, TabularPredictor
#import autogluon.core as ag
import sys 

std_slc = StandardScaler()
pca = decomposition.PCA()


def molmethod(X,mm,rdn):
    def MSE(y_true,y_pred):
        mse = mean_squared_error(y_true, y_pred)
        return mse
    def R2(y_true,y_pred):    
         r2 = r2_score(y_true, y_pred)
         return r2
    def two_score(y_true,y_pred):    
        MSE(y_true,y_pred) #set score here and not below if using MSE in GridCV
        score = R2(y_true,y_pred)
        print('R2: %2.3f   MSE: %2.3f' %(score,MSE(y_true,y_pred)))
        return score
    def two_scorer():
        return make_scorer(two_score, greater_is_better=True) # change for false if using MSE

    scorer = {'MAE': 'neg_mean_absolute_error', 'r2': 'r2'}  
    score='neg_mean_absolute_error'
    scoremk = {'MAE': 'neg_mean_absolute_error', 'r2': make_scorer(r2_score)}

    reg=ensemble.RandomForestRegressor(random_state=rdn)
    mm=mm.upper() 
    if (mm == "LGBM") :
        parameters = {'learning_rate': [0.01,0.001,0.0001,0.00001],
        'max_depth': [4,6,8],'num_leaves': [20,30,40]}
        regt = LGBMRegressor(random_state=rdn)
        reg=GridSearchCV(regt, param_grid=parameters, scoring='neg_mean_absolute_error', cv=3)
        
    if (mm == "RF"):
        parameters = {'min_samples_leaf': [1,2,3]} #'max_depth': [4,6,8],'min_samples_leaf': [1,2,3],
        regt = RandomForestRegressor()
        reg=GridSearchCV(regt, param_grid=parameters, scoring='neg_mean_absolute_error', cv=3)
    if (mm == "RFR"):
        parameters = {'min_samples_leaf': [1,2,3]} #'max_depth': [4,6,8],'min_samples_leaf': [1,2,3],
        regt = RandomForestRegressor()
        reg=GridSearchCV(regt, param_grid=parameters, scoring='r2', cv=3)
    if (mm == "RFS"):
        parameters = {'min_samples_leaf': [1,2,3]} #'max_depth': [4,6,8],'min_samples_leaf': [1,2,3],
        regt = RandomForestRegressor()
        reg=GridSearchCV(regt, param_grid=parameters, scoring='neg_mean_absolute_error', cv=3)
    if (mm == "RFT"):
        parameters = {'min_samples_leaf': [1,2,3]} #'max_depth': [4,6,8],'min_samples_leaf': [1,2,3],
        regt = RandomForestRegressor()
        reg=GridSearchCV(regt, param_grid=parameters,scoring=two_scorer(),cv=3, refit='r2') 
    if (mm == "RFMK"):
        parameters = {'min_samples_leaf': [1,2,3]} #'max_depth': [4,6,8],'min_samples_leaf': [1,2,3],
        regt = RandomForestRegressor()
        reg=GridSearchCV(regt, param_grid=parameters,scoring=scoremk,cv=3, refit='r2') 
    if (mm == "RFRF"):
        parameters = {'min_samples_leaf': [1,2,3]} #'max_depth': [4,6,8],'min_samples_leaf': [1,2,3],
        #'min_samples_split': [2,4,6]}  #{'max_depth': 4, 'min_samples_leaf': 1, 'min_samples_split': 2}
        regt = RandomForestRegressor()
        reg=GridSearchCV(regt, param_grid=parameters, scoring=scorer, cv=3,refit='r2')

    if (mm == "LR") :
        parameters = {"fit_intercept": [True, False],
        "solver": ['svd', 'cholesky', 'lsqr', 'sparse_cg', 'sag', 'saga']}
        regt = linear_model.LinearRegression()
        reg=GridSearchCV(regt, param_grid=parameters, scoring='neg_mean_absolute_error', cv=3)
        
    if (mm == "DECISIONTREE") :
        n_components = list(range(1,X.shape[1]+1,1))
        criterion = ['mse','mae','friedman_mse','poisson']
        max_depth = [2,4,6,8,10,12]
        regt = tree.DecisionTreeRegressor()
        parameters = dict(pca__n_components=n_components,
                      regt__criterion=criterion,
                      regt__max_depth=max_depth)
        pipe = Pipeline(steps=[('std_slc', std_slc),('pca', pca),('regt',regt)])
        reg=GridSearchCV(pipe, parameters)
    if (mm == "SGD") :
        regt = linear_model.SGDRegressor(random_state=rdn,max_iter=10000)
        reg=GridSearchCV(regt, param_grid=parameters, scoring='neg_mean_absolute_error', cv=3)
    if (mm == "LASSO") :
        regt =linear_model.Lasso(random_state=rdn,max_iter=10000)
        reg=GridSearchCV(regt, param_grid=parameters, scoring='neg_mean_absolute_error', cv=3)
    if (mm == "ENETCV") :
        regt = linear_model.ElasticNetCV(random_state=rdn,max_iter=10000)
        reg=GridSearchCV(regt, param_grid=parameters, scoring='neg_mean_absolute_error', cv=3)
    if (mm == "ENET") :
        regt = linear_model.ElasticNet(random_state=rdn,max_iter=10000)
        reg=GridSearchCV(regt, param_grid=parameters, scoring='neg_mean_absolute_error', cv=3)
    if (mm == "RIDGE") :
        regt = linear_model.Ridge(random_state=rdn,max_iter=20000)
        reg=GridSearchCV(regt, param_grid=parameters, scoring='neg_mean_absolute_error', cv=3)
    if (mm=="KRR"):
        parameters = {'alpha':np.logspace(-10, 3, 53),'gamma':np.logspace(-8, 3, 45)}
                      # #default 1.0 , #default none
        regt=KernelRidge(kernel=#'poly','sigmoid , 
            'rbf')
        reg=GridSearchCV(regt, param_grid=parameters, scoring=scorer, cv=3,refit='r2')
    if (mm == "SVR") :
        parameters = {'C':[1,5,10,15,25,50],
        'gamma': [0.0001,0.0005,0.001,0.005,0.01,0.05,0.1]}
        regt=svm.SVR()
        reg=GridSearchCV(regt, param_grid=parameters, scoring='neg_mean_absolute_error', cv=5)
    if (mm == "NUSVR") :
        regt = svm.NuSVR()
        reg=GridSearchCV(regt, param_grid=parameters, scoring='neg_mean_absolute_error', cv=5)
    if (mm == "KNN") :
        regt = KNeighborsRegressor()
        reg=GridSearchCV(regt, param_grid=parameters, scoring='neg_mean_absolute_error', cv=5)
    if (mm == "GB") :
        regt = ensemble.GradientBoostingRegressor(random_state=rdn)
        reg=GridSearchCV(regt, param_grid=parameters, scoring='neg_mean_absolute_error', cv=5)
    if (mm == "AUTOSK") :
        reg = autosklearn.regression.AutoSklearnRegressor(
        time_left_for_this_task=300,
        memory_limit=None,
        resampling_strategy='cv',resampling_strategy_arguments={'folds': 5})
    return reg 

def dstestmain(x_train,y_train,x_test,y_test,rdn,method):  #
    reg = molmethod(x_train,method,rdn)
    reg.fit(x_train, y_train)
    print(f"{reg.best_params_} {reg.best_score_:0.03f}")
    print(f"{reg.cv_results_}")
    best_params = reg.best_params_
    model = reg.best_estimator_
    score = reg.best_score_
    resultsr2=reg.cv_results_['mean_test_r2']
    resultsmae=reg.cv_results_['mean_test_MAE']
    params=reg.cv_results_['params']
    print(method,nd)
    for pm,rr2,mmae in zip(params,resultsr2,resultsmae):
        print(f"{pm} {rr2:0.03f} {mmae:0.03f}")


    y_train_predicted = reg.predict(x_train)
    r2_train = r2_score(y_train, y_train_predicted)
    #print("Data: {} Method: {} Train R2: {} ".format(nd, method, r2_train))
    y_test_predicted = reg.predict(x_test)
    r2_test = r2_score(y_test, y_test_predicted)
    mae=mean_absolute_error(y_test, y_test_predicted)
    #dd=np.vstack((y_test, y_test_predicted))
    #da=pd.DataFrame(dd)
    #da=da.T
    #ds=da[da[0]>5.0]
    #mae=mean_absolute_error(ds[0],ds[1])
    #return r2_train, r2_test, mae
    return y_train, y_train_predicted,y_test, y_test_predicted,r2_train, r2_test, mae

v_col=["m","n","VLL","VHH","VLH","VHL"]
vdata = pd.read_csv("v.out", header=None,names=v_col, delim_whitespace=True)
fr=0.95

#tnd=['c3dmrs.dat','cip.dat','dmrc3d.dat','dmre3d.dat','dmrm3d.dat','dmrp3d.dat','dmru3d.dat','dmrv3d.dat','dpic3d.dat','dpie3d.dat','dpim3d.dat','dpip3d.dat','dpiu3d.dat','dpiv3d.dat',
#'dsall.dat' has problem
tnd=['c3dmrs.dat','cip.dat','dmrc3d.dat','dmre3d.dat','dmrm3d.dat','dmrp3d.dat','dmru3d.dat',\
'dmrv3d.dat','dpic3d.dat','dpie3d.dat','dpim3d.dat','dpip3d.dat','dpiu3d.dat','dpiv3d.dat',\
'e3dmrs.dat','frame4.dat','m3dmrs.dat','mgeomty.dat','neardisfft.dat',\
'p3dmrs.dat','pic3d.dat','picip.dat','pie3d.dat','pigeomty.dat','pim3d.dat','pip3d.dat',\
'piu3d.dat','piv3d.dat','refpen.dat','rmol.dat','rpi.dat','u3dmrs.dat','v3dmrs.dat']

j=sys.argv[1]
nd=j+'.dat'
j='KRR'

data = pd.read_csv(nd, header=None,delim_whitespace=True,dtype=np.float64)
dim=len(data.iloc[1])
dft = data.iloc[:,2:dim] # skip m n
df=pd.concat([dft,vdata.VLL.abs()],axis=1)
df = df.loc[:,(df != df.iloc[0]).any()]
sets=np.empty(1)
for rdn in range(1):
    label="VLL"
    tdf=df.copy()
    train_dataset = tdf.sample(frac=fr, random_state=rdn)  
    test_dataset = tdf.drop(train_dataset.index)
    train_labels = train_dataset.pop(label) # kick out label column from train_stats and give to train_labels
    x_train=train_dataset
    y_train=train_labels
    test_labels = test_dataset.pop(label)
    x_test=test_dataset
    y_test =test_labels
    y_train, y_train_predicted,y_test, y_test_predicted,r2_train, \
    r2_test, mae=dstestmain(x_train,y_train,x_test,y_test,rdn,j)
    sets=np.append(sets.copy(),mae)
    print("Fr {} Data {} Method {} Random {} Train_R2 {} Test_R2 {} MAE {}".format(fr, nd, j,rdn, r2_train, r2_test,mae))

sets=np.delete(sets,0,axis=0)
print("Fr {} Data {} Method {} average MAE {}".format(fr, nd, j,sets.mean()))

'''
for i in {1..128}
do
t=0.000000001*$i*5
sed "s/xxxxxxxxx/$t/" v0.py > $i.py
nohup python $i.py > $i.log &
done
'''
	
