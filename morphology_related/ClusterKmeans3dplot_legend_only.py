# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 23:30:18 2024

@author: 123
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 13:10:08 2024
输入文件包含：
@author: 123
"""

import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

Cluster_Scaled = 'AaAa_Cluster_Scaled.xlsx'
ClusterCenterInformation = 'ClusterCenterInformation.xlsx'
Molecule = 'ICB'

# 读取数据
cluster_data = pd.read_excel(Cluster_Scaled, usecols=[4, 5, 6, 7], names=['Average', 'Farthest', 'Angle', 'Label'])
center_data = pd.read_excel(ClusterCenterInformation, usecols=[0, 1, 2, 3], names=['Average', 'Farthest', 'Angle', 'Label'])

# 设置颜色
n_clusters = len(cluster_data['Label'].unique())  # 获取聚类的数量
colors = plt.cm.viridis(np.linspace(0, 1, n_clusters))  # 生成对应数量的颜色

# 创建空白画布
fig, ax = plt.subplots(figsize=(1.77, 1.77))

# 关闭坐标轴
ax.axis('off')

# 创建自定义图例句柄和标签
handles = []
labels = []

for i in range(n_clusters):
    handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[i], markersize=3, label=f'S{i}'))
    labels.append(f'S{i}')

# 添加类中心点
handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=3, label='Cen'))
labels.append('Cen')

# 设置图例的字体和字号
legend_font = FontProperties(family='Times New Roman', size=6.5, weight='bold')

# 创建图例，并设置其属性
legend = ax.legend(
    handles=handles,
    labels=labels,
    frameon=False,  # 移除图例边框
    prop=legend_font,  # 设置图例文字的字体属性
    loc='upper right',  # 图例位置
    handlelength=0.5,  # 控制图例图标长度
    handletextpad=0.2,  # 控制图例图标与文本之间的距离
    columnspacing=0.5,  # 控制图例列之间的间距
    borderaxespad=0.5,  # 控制图例与坐标轴之间的间距
    labelspacing=0.1  # 控制图例条目之间的垂直间距
)

# 保存图像
output_filename = f'{Molecule}_{Cluster_Scaled[:4]}_cluster_legend_only2.png'
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"Image saved as {output_filename}")

# 显示图像
plt.show()