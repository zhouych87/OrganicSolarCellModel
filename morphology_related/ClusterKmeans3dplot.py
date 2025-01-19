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

Cluster_Scaled = 'eC9AD_Cluster_Scaled.xlsx'
ClusterCenterInformation = 'ClusterCenterInformation.xlsx'
Molecule = 'eC9'

# 定义图像尺寸
initial_fig_size = 1.77  # 初始图像尺寸（英寸）
fig_width = initial_fig_size   # 增加图像宽度，以便更好地展示细节
fig_height = initial_fig_size   # 增加图像高度，以便更好地展示细节

# 读取数据
cluster_data = pd.read_excel(Cluster_Scaled, usecols=[4, 5, 6, 7], names=['Average', 'Farthest', 'Angle', 'Label'])
center_data = pd.read_excel(ClusterCenterInformation, usecols=[0, 1, 2, 3], names=['Average', 'Farthest', 'Angle', 'Label'])

# 检查数据是否正确读取
print("Cluster Data:")
print(cluster_data.head())
print("\nCenter Data:")
print(center_data.head())

# 设置颜色
n_clusters = len(cluster_data['Label'].unique())  # 获取聚类的数量
colors = plt.cm.viridis(np.linspace(0, 1, n_clusters))  # 生成对应数量的颜色

# 创建3D图
fig = plt.figure(figsize=(fig_width, fig_height), facecolor='White')  # 设置背景透明
ax = fig.add_subplot(111, projection='3d', facecolor='White')  # 设置背景透明

# 绘制散点图
for i in range(n_clusters):
    subset = cluster_data[cluster_data['Label'] == i]
    ax.scatter(subset['Average'], subset['Farthest'], subset['Angle'], c=[colors[i]], label=f'S{i}', s=1)  # 调整点的大小

print('plot scatter is finished...')

# 绘制类中心点
for idx, row in center_data.iterrows():
    ax.scatter(row['Average'], row['Farthest'], row['Angle'], c='red', s=4, marker='o', label='Cen' if idx == 0 else "")
    #ax.scatter(row['Average'], row['Farthest'], row['Angle'], c='red', s=5, marker='o', edgecolors='k', label='Cen' if idx == 0 else "")

# 设置坐标轴范围和刻度
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_zlim(0, 1)
ax.set_xticks(np.arange(0, 1.1, 0.3))
ax.set_yticks(np.arange(0, 1.1, 0.3))
ax.set_zticks(np.arange(0, 1.1, 0.2))

# # 添加半刻度线
# for tick in np.arange(0.1, 1, 0.2):
#     ax.plot([tick, tick], [0, 1], [0, 0], color='gray', linestyle='--', linewidth=0.5)
#     ax.plot([0, 1], [tick, tick], [0, 0], color='gray', linestyle='--', linewidth=0.5)
#     ax.plot([0, 1], [0, 0], [tick, tick], color='gray', linestyle='--', linewidth=0.5)

# 定义字体属性
font = {
    'family': 'Times New Roman',  # 字体系列，例如 'serif', 'sans-serif', 'monospace' 等
    'color':  'black',  # 字体颜色
    #'weight': 'bold',  # 字体粗细
    'size': 6.5,  # 字体大小
}

# 设置坐标轴标签
ax.set_xlabel('Average', fontdict=font, labelpad=-13)  # 调整标签与刻度的距离
ax.set_ylabel('Farthest', fontdict=font, labelpad=-13)  # 调整标签与刻度的距离
ax.set_zlabel('Angle', fontdict=font, labelpad=-10)  # 调整标签与刻度的距离


# 设置刻度标签的字体和字号
for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
    for label in axis.get_majorticklabels():
        label.set_fontsize(6.5)
        label.set_color('black')
        #label.set_weight('bold')
        label.set_family('Times New Roman')
        
# 调整刻度标签与刻度线之间的距离
ax.tick_params(axis='x', pad=-6)
ax.tick_params(axis='y', pad=-6)
ax.tick_params(axis='z', pad=-2)

# # 调整坐标轴线条粗细
# for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
#     axis._axinfo['axisline']['linewidth'] = 5  # 设置线条粗细

# # 设置图例的字体和字号
# legend_font = FontProperties(family='Times New Roman', size=6.5, weight='bold')
# # legend = ax.legend(frameon=False, prop=legend_font, loc='upper right')  # 移除图例边框
# # #legend = ax.legend(frameon=False, prop=legend_font, loc=(0.8, 0.8))  # 移除图例边框

# # 创建自定义图例句柄和标签
# handles = []
# labels = []

# for i in range(n_clusters):
#     handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[i], markersize=3, label=f'S{i}'))
#     labels.append(f'S{i}')

# # 添加一个空行
# handles.append(Line2D([0], [0], marker='', color='w', label=''))
# labels.append('')

# # 添加类中心点
# handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=3, label='Cen'))
# labels.append('Cen')

# # 创建图例，并设置其属性
# legend = ax.legend(
#     handles=handles,
#     labels=labels,
#     frameon=False,  # 移除图例边框
#     prop=legend_font,  # 设置图例文字的字体属性
#     loc='upper right',  # 图例位置
#     handlelength=0.5,  # 控制图例图标长度
#     handletextpad=0.2,  # 控制图例图标与文本之间的距离
#     columnspacing=0.5,  # 控制图例列之间的间距
#     borderaxespad=0.5,  # 控制图例与坐标轴之间的间距
#     labelspacing=0.1  # 控制图例条目之间的垂直间距
# )


# # 设置图形的整体左移和上移
# left_offset = 1  # 左移英寸
# top_offset = 1   # 上移英寸

# # 计算左移和上移的比例
# left_ratio = left_offset / fig_width
# top_ratio = top_offset / fig_height

# # 调整子图的位置
# fig.subplots_adjust(left=left_ratio, bottom=0, right=1, top=1 - top_ratio)


# 保存图像
output_filename = f'{Molecule}_{Cluster_Scaled[:4]}_cluster_3dplot.png'
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"Image saved as {output_filename}")

# 显示图像
plt.show()