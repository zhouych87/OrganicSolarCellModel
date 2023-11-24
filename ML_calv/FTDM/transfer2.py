import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


data=pd.read_csv('transfer.csv')
dat=data.iloc[:,1:]
dat.index=data.iloc[:,0]
cname=['A', '1@B', '2@B', '3@B', '4@B', '5@B', '6@B', '7@B', '8@B', '9@B', '10@B', 'C', 'Film']
dat.columns=cname
dat.index=cname
plt.ion()

f, ax = plt.subplots(figsize=(12, 9))
sns.heatmap(dat, annot=True, annot_kws={"fontsize":12},fmt=".2f", linewidths=.5, ax=ax,cmap="PiYG_r")
#sns.heatmap(dat)

plt.xlabel('Train set',fontsize=15, color='k') #x轴label的文本和字体大小
plt.ylabel('Test set',fontsize=15, color='k') #y轴label的文本和字体大小

cax = plt.gcf().axes[-1]
cax.tick_params(labelsize=15)
cbar = ax.collections[0].colorbar
cbar.set_label('MAE (meV)',size=15)

plt.xticks(size=12)
plt.yticks(size=12)


plt.savefig('transfer.png', dpi=300)

#plt.show()
plt.ioff()


