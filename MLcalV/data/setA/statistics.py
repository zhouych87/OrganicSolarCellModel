import numpy as np
import pandas as pd
import os
import sys
#from scipy import stats
 
#os.system('sed -i "s/Train R2/trainr2/g" test.log')

mth=['knn','gb','decisiontree','rf']
named=["ref","mol","u3d","e3d"]
cfr=[95,90,80,70,60,50,40,30,20,15,10,'05','04','03','02','01']

for fr in cfr:
    logn=str(fr)+'t3.log'
    logd= pd.read_csv(logn, header=None, delim_whitespace=True)
    logd.insert(loc=12, column=12, value=int(fr))  
    if (fr == 95):
        df=logd.copy()
    else:
        df=pd.concat([df,logd],axis=0,ignore_index=True)


#df = df[df[LABEL] != 0]  # Remove zero small 
#df = df[abs(df[LABEL]) > 1] #0.001 0.11; 0.01-0.29,0.1-0.17,1.0-0.11; 0-0.29

#          0    1        2      3        4   5          6         7         8         9     10         11  12
#0      Data:  ref  Method:     rf  Random:   0  Train_R2:  0.999970  Test_R2:  0.999839  MAE:   0.265755  95

#data_range = data.loc[data['x3'] >= 2]    


# 1 3 5 7 9 11 12 
cfr=[95,90,80,70,60,50,40,30,20,15,10,5,4,3,2,1]
trainsize=pd.DataFrame.empty()

for i in mth:
    for j in named:
        print(i,j,end=' ') 
        tmean=[]
        for fr in cfr:
            tdf=df.loc[(df[1]==j) & (df[3]==i) & (df[12]==fr)]
            tmae=tdf[11]
            #tmean.append(tmae.mean())
            print(" %.3f"%(tmae.mean()),end=' ')
        print(' ')


cfr=[95,90,80,70,60,50,40,30,20,15,10,5,4,3,2,1]
sfr=['method','Descriptors','95','90','80','70','60','50','40','30','20','15','10','5','4','3','2','1']
dfr=['95','90','80','70','60','50','40','30','20','15','10','5','4','3','2','1']
trainsize=pd.DataFrame(columns=sfr)

for i in mth:
    for j in named:
        sets=[]
        sets.append(i)
        sets.append(j)
        print(i,j,end=' ') 
        for fr in cfr:
            tdf=df.loc[(df[1]==j) & (df[3]==i) & (df[12]==fr)]
            tmae=tdf[11]
            #tmean.append(tmae.mean())
            sets.append(tmae.mean())
            print(" %.3f"%(tmae.mean()),end=' ')
        dset=pd.DataFrame(sets)
        dset=dset.T
        dset.columns=sfr
        trainsize=pd.concat((trainsize,dset))
        print(' ')


train=trainsize.copy()
import matplotlib.pyplot as plt
plt.ion()
plt.figure(figsize=(12.8,9.6))

clr=pd.DataFrame(['black','orange','cyan','red'])
clr=clr.T
clr.columns=['knn','gb','decisiontree','rf']

#mrk=pd.DataFrame(['+','x','1','2'])
mrk=pd.DataFrame(['s','o','>','<'])
mrk=mrk.T
mrk.columns=["ref","mol","u3d","e3d"]

plt.close()
fig, ax = plt.subplots()
for i in mth:
    xx=list(cfr)
    #yyy=pd.DataFrame()
    yyyy=[]
    for j in named:
        y=train[(train['method']==i) & (train['Descriptors']==j)]
        yy=y[dfr]
        #yyy=pd.concat((yyy,yy))
        yyy=list(yy.iloc[0])
        yyyy=[yyyy,yyy]
        #p1=plt.plot(xx,yyy, color=clr[i][0], marker=mrk[j][0],fillstyle = 'none', linestyle='dashed', linewidth=1, markersize=6)
        p1=plt.loglog(xx,yyy, color=clr[i][0], marker=mrk[j][0],fillstyle = 'none', linestyle='dashed', linewidth=1, markersize=6)
    #yyyy=np.array(yyy.T)


plt.xlabel('Training Fraction (%)')
plt.ylabel('Average MAE (meV)')
ax.legend((xx,yyy),['Leberer','This work','UMor','EMor'])

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], linestyle='dashed',color='black', lw=1),
                Line2D([0], [0], linestyle='dashed',color='orange',lw=1),
                Line2D([0], [0], linestyle='dashed',color='cyan', lw=1),
                Line2D([0], [0], linestyle='dashed',color='red', lw=1)]

#lines = ax.plot(data)
leg1=ax.legend(custom_lines, ['KNN','GB','DT','RF'], loc='lower center')
'''
from matplotlib.lines import Line2D
custom_lines1 = [Line2D([0], [0],marker='+', lw=0),
                Line2D([0], [0], marker='x',lw=0),
                Line2D([0], [0], marker='1', lw=0),
                Line2D([0], [0], marker='2', lw=0)]

#lines = ax.plot(data)
leg2=ax.legend(custom_lines1,['Leberer','This work','UMor','EMor'], loc='upper right')
'''


from matplotlib.lines import Line2D
custom_lines1 = [Line2D([0],[0],marker='s',fillstyle = 'none', lw=0),
                Line2D([0], [0], marker='o',fillstyle = 'none',lw=0),
                Line2D([0], [0], marker='>',fillstyle = 'none', lw=0),
                Line2D([0], [0], marker='<',fillstyle = 'none', lw=0)]

#lines = ax.plot(data)
leg2=ax.legend(custom_lines1,['Leberer','This work','UMor','EMor'], loc='lower left')

ax.add_artist(leg1)
