import numpy as np
import pandas as pd
import os
import sys
from time import sleep 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import mpl_toolkits.axisartist as axisartist


custom_lines = [Line2D([0], [0], linestyle='None',marker='s',color='gray', lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='blue', lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='cyan',lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='teal',lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='olive', lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='olivedrab', lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='orange', lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='sienna', lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='red', lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='purple', lw=1)]

clr=['gray','blue','cyan','teal','olive','olivedrab','orange','sienna','red','purple']

nd=['refpen','cip','v3dmrs','p3dmrs','c3dmrs','u3dmrs','e3dmrs','m3dmrs']
lnd=['GEO','CIP','VMor','PMor','CMor','UMor','EMor','MMor']

nd=['refpen','cip','u3dmrs','c3dmrs','dmru3d','dmrc3d','piu3d','pic3d','dpiu3d','dpic3d']
lnd=['GEO','CIP','UM','CM','UIM','CIM','UM-PI','CM-PI','UIM-PI','CIM-PI']


#nd=['cip','dpiu3d','dpic3d','piu3d','u3dmrs','dmru3d']
logn='sumsize.log'
logd= pd.read_csv(logn, header=None, delim_whitespace=True)
df=logd.copy()


plt.close()
plt.ion()

fig=plt.figure(figsize=(12,4.0),constrained_layout=True)
plt.rcParams['font.size'] = '18'

ax0=fig.add_subplot(111)

widthfr=0.1
sets=np.empty(1)
setm=np.array([0,0])
for i in range(len(nd)):
    tdf=df.loc[(df[3]==nd[i]) & (df[1]==0.95)]
    tmae=tdf[11] ##8
    sets=np.append(sets.copy(),tmae.mean())
    setm=np.vstack((setm.copy(),np.array([tmae.mean()-tmae.min(),tmae.max()-tmae.mean()])))

setm=np.delete(setm,0,axis=0)
sets=np.delete(sets,0,axis=0)
#setyerr=setm.reshape(2,len(cfr))
#plt.bar(tcfr,sets,yerr=setm.T,width=widthfr,color=clr[i])
plt.bar(lnd,sets,yerr=setm.T,color=clr,label=lnd) #width=10*widthfr,
#plt.errorbar(tcfr,sets,yerr=setm.T,lw=4,linestyle='None',ecolor=clr[i],marker='*',markeredgecolor='black')
#sns.barplot(data=df, x="island", y="body_mass_g", hue="sex")

plt.xlabel('Descriptor')
plt.ylabel('MAE (meV)')
plt.title("KRR")

plt.grid(True)
#plt.legend(custom_lines, lnd, loc='upper right')
figname='dtest_MAE_bar.png'
plt.savefig(figname,format='png')

