import numpy as np
import pandas as pd
import os
import sys
#from scipy import stats
 
####################################
############ 3D REVISE #############
####################################

named=["u3d","c3d","p3d","m3d","v3d","e3d"]
sets=['snsr','cnsr','snr','cnr','snsqrtr','cnsqrtr','sn','cn','abssnsr','abscnsr','abssnr','abscnr','abssnsqrtr','abscnsqrtr','abssn','abscn']
mth=['knn','gb','decisiontree','rf']

for i in named:  
    for set in sets:
        logn=set+'/test.log'
        logd= pd.read_csv(logn, header=None, delim_whitespace=True)
        logd.insert(loc=12, column=12, value=set)  
        if (set == 'snsr'):
            df=logd.copy()
        else:
            df=pd.concat([df,logd],axis=0,ignore_index=True)


################ revise 3D  large  ####################
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.lines import Line2D
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


custom_lines = [Line2D([0], [0], linestyle='None',marker='s',color='b', lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='cyan',lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='green',lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='orange', lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='red', lw=1),
                Line2D([0], [0], linestyle='None',marker='s',color='sienna', lw=1)]
widthfr=0.4
clr= matplotlib.colors.Colormap('rainbow', N=16)
clr = plt.cm.rainbow(np.linspace(0, 1, 6))
clr=['b','cyan','green','orange','r','sienna']

#tcfr=np.array({i for i in range(16)}).tolist
#cfr=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
#tcfr=np.array(cfr)


for j in mth:
    plt.close()
    plt.ion()
    #plt.figure(figsize=(12.8,9.6))
    fig, ax0 = plt.subplots(figsize=(12.8,9.6))
    
    plt.rcParams['font.size'] = '24'
    for label in (ax0.get_xticklabels() + ax0.get_yticklabels()):
        label.set_fontsize(24)
    tcfr=np.array([i for i in range(16)])
    tcfr=tcfr+0.8
    for ii in range(6):  
        i=named[ii]
        seta=np.empty(1)
        setm=np.empty((1,2))
        for set in sets:
            tdf=df.loc[(df[1]==i) & (df[12]==set)& (df[3]==j)]
            tmae=tdf[11]
            seta=np.append(seta.copy(),tmae.mean())
            setm=np.append(setm.copy(),[[tmae.min(),tmae.max()]],axis=0)
            
        setm=np.delete(setm,0,axis=0)
        seta=np.delete(seta,0,axis=0)
        setyerr=setm.reshape(2,len(sets))
        #plt.bar(tcfr,sets,yerr=setyerr,width=widthfr,color=clr[jj],alpha=0.7)
        plt.errorbar(tcfr,seta,yerr=setyerr,lw=4,linestyle='None',ecolor=clr[ii],marker='*',markeredgecolor='black', markerfacecolor='white') #ecolor=clr[jj]
        #plt.errorbar(tcfr,sets,yerr=setyerr,lw=5,linestyle='None',ecolor=clr[jj]) 
        tcfr=tcfr+0.1
        ax0.xaxis.set_major_locator(MultipleLocator(1)) #
        ax0.yaxis.set_major_locator(MultipleLocator(2))
        plt.ylim(0,30)
        plt.xlim(0.5,16.5)
        plt.tight_layout()
        plt.title(j)
        plt.grid(True)
        plt.xlabel(sets)
        plt.ylabel('MAE (meV)')
    plt.legend(custom_lines, ["UMor","CMor","PMor","MMor","VMor","EMor"], loc='upper right')
    figname=j+'_rev3d.png'  # 10 20 30 70 80 , 40 90 
    plt.savefig(figname,format='png')



################ small 3D   ####################


for j in mth:
    plt.close()
    plt.ion()
    #plt.figure(figsize=(12.8,9.6))
    fig, ax0 = plt.subplots()
    
    plt.rcParams['font.size'] = '12'
    for label in (ax0.get_xticklabels() + ax0.get_yticklabels()):
        label.set_fontsize(12)
    tcfr=np.array([i for i in range(16)])
    tcfr=tcfr+0.8
    for ii in range(6):  
        i=named[ii]
        seta=np.empty(1)
        setm=np.empty((1,2))
        for set in sets:
            tdf=df.loc[(df[1]==i) & (df[12]==set)& (df[3]==j)]
            tmae=tdf[11]
            seta=np.append(seta.copy(),tmae.mean())
            setm=np.append(setm.copy(),[[tmae.min(),tmae.max()]],axis=0)
            
        setm=np.delete(setm,0,axis=0)
        seta=np.delete(seta,0,axis=0)
        setyerr=setm.reshape(2,len(sets))
        #plt.bar(tcfr,sets,yerr=setyerr,width=widthfr,color=clr[jj],alpha=0.7)
        plt.errorbar(tcfr,seta,yerr=setyerr,lw=2,linestyle='None',ecolor=clr[ii],marker='*',markeredgecolor='black', markerfacecolor='white') #ecolor=clr[jj]
        #plt.errorbar(tcfr,sets,yerr=setyerr,lw=5,linestyle='None',ecolor=clr[jj]) 
        tcfr=tcfr+0.1
        ax0.xaxis.set_major_locator(MultipleLocator(1)) #
        ax0.yaxis.set_major_locator(MultipleLocator(2))
        plt.ylim(0,30)
        plt.xlim(0.5,16.5)
        plt.tight_layout()
        plt.title(j)
        plt.grid(True)
        plt.xlabel(sets)
        plt.ylabel('MAE (meV)')
    plt.legend(custom_lines, ["UMor","CMor","PMor","MMor","VMor","EMor"], loc='upper right')
    figname='s'+j+'_rev3d.png'  # 10 20 30 70 80 , 40 90 
    plt.savefig(figname,format='png')


