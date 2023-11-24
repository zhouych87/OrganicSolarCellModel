import numpy as np
import pandas as pd
import os
import sys
import _tkinter
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import mpl_toolkits.axisartist as axisartist


lnd=['CIP','UM','UIM-PI','FTDM']
loglist=['trainCIP.log','trainUM.log','trainUIM.log','trainFTDM.log']
clr=['cyan','green','orange','r'] #'cyan','sienna','b','green','orange','r'

custom_lines = [Line2D([0], [0], linestyle='None',marker='s',color=clr[0], lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color=clr[1],lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color=clr[2], lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color=clr[3], lw=1)]

plt.close()
plt.ion()

fig=plt.figure(figsize=(16,8.0),constrained_layout=True)
plt.rcParams['font.size'] = '18'

ax0=fig.add_subplot(111)
cfr=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
widthfr=2


tcfr=np.linspace(10,90,9)
tcfr=tcfr-1.5*widthfr

for i in range(len(loglist)):
    logd= pd.read_csv(loglist[i], header=None, delim_whitespace=True)
    df=logd.dropna()
    sets=np.empty(1)
    setm=np.array([0,0])
    for fr in cfr:
        tdf=df.loc[(np.float64(df[2])==fr)]
        tmae=tdf[10] ##8
        sets=np.append(sets.copy(),tmae.mean())
        setm=np.vstack((setm.copy(),np.array([tmae.mean()-tmae.min(),tmae.max()-tmae.mean()])))
    
    setm=np.delete(setm,0,axis=0)
    sets=np.delete(sets,0,axis=0)
    plt.bar(tcfr,sets,yerr=setm.T,color=clr[i],width = 2.0,align='center', alpha=1.0, ecolor='black', capsize=4)
    tcfr=tcfr+widthfr


ax0.xaxis.set_major_locator(MultipleLocator(10)) 
ax0.yaxis.set_major_locator(MultipleLocator(1)) 
#plt.xscale('log')
#plt.yscale('log')
plt.ylim(0,10)
plt.grid(True)
plt.xlabel('Training portion (%)')
plt.ylabel('MAE (meV)')
plt.legend(custom_lines, lnd, loc='upper right')

figname='size_tet_MAE.png'
plt.savefig(figname,format='png',dpi=300)

