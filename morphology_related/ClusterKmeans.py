# -*- coding: utf-8 -*-
"""
Spyder Editor
只用Kmeans算法聚类，并且只支持一个聚类数计算，只为了后期画图提供数据
This is a temporary script file.
"""

import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
import numpy as np
from scipy.spatial.distance import cdist

# 读取原始数据文件
file_path = 'eC9AD.xlsx'
data = pd.read_excel(file_path, header=None)

# 提取需要处理的列
columns_to_process = [4, 5, 6]  # 从第5列到第7列

# 确保列索引在有效范围内
if max(columns_to_process) >= data.shape[1]:
    raise IndexError("列索引超出数据框的列数范围")

X = data.iloc[:, columns_to_process].values

# 对第10列数据进行弧度转换和余弦变换
X_radians_9 = np.radians(X[:, 2])  # 注意这里使用的是第10列，即索引为2
X_cos_9 = np.cos(X_radians_9).reshape(-1, 1)

# 将两组距离特征分别做归一化处理
X_part1 = X[:, 0]  # 修改这里以正确引用第一列
X_part2 = X[:, 1]  # 修改这里以正确引用第二列
# 归一化
scaler = MinMaxScaler()
X_scaled_part1 = scaler.fit_transform(X_part1.reshape(-1, 1))  # 第5列
X_scaled_part2 = scaler.fit_transform(X_part2.reshape(-1, 1))  # 第6列
# 合并三个维度的数据
X_pca_scaled = np.hstack((X_scaled_part1, X_scaled_part2, X_cos_9))

# 提取原始数据表的前4列
original_columns = data.iloc[:, :4].values
# 将原始数据表的前4列与处理后的特征合并
X_pca_scaled_with_original = np.hstack((original_columns, X_pca_scaled))
# 创建DataFrame
X_pca_scaled_df = pd.DataFrame(X_pca_scaled_with_original, 
                               columns=['Original_1', 'Original_2', 'Original_3', 'Original_4', 'Average', 'Farthest', 'Angle'])
# 保存到Excel文件
X_pca_scaled_df.to_excel(f'{file_path[:-5]}_Feature_Scaled.xlsx', index=False)

# 聚类设置
n_clusters = 5  # 聚类数
kmeans = KMeans(n_clusters=n_clusters, random_state=0, n_init=10)
kmeans.fit(X_pca_scaled)
labels = kmeans.labels_
centers = kmeans.cluster_centers_

# 将聚类标签添加到原始数据表的最后一列
data['Cluster'] = labels

# 保存带有聚类标签的原始数据表
data.columns = [f'Feature_{i+1}' for i in range(data.shape[1] - 1)] + ['Cluster']
data.to_excel(f'{file_path[:-5]}_Cluster.xlsx', index=False)

data_Cluster_and_Scaled_Part1 = pd.read_excel(f'{file_path[:-5]}_Feature_Scaled.xlsx', usecols=range(7))
data_Cluster_and_Scaled_Part2 = pd.read_excel(f'{file_path[:-5]}_Cluster.xlsx', usecols=[7])
# 合并两个部分
data_Cluster_and_Scaled = pd.concat([data_Cluster_and_Scaled_Part1, data_Cluster_and_Scaled_Part2], axis=1)
# 保存到Excel文件
data_Cluster_and_Scaled.to_excel(f'{file_path[:-5]}_Cluster_Scaled.xlsx', index=False)

print("聚类完成，结果已保存。")

# 计算每个聚类中心到轴0值的距离
def distance_to_axes(center):
    shortest_dist = center[0]
    farest_dist = center[1]
    angle_dist = center[2]  # 因为要找离Angle轴0值最远的，所以用1减去值
    return shortest_dist, farest_dist, angle_dist

# 找出符合条件的类
best_cluster_idx = None
min_shortest_farest_sum = float('inf')
max_angle_dist = float('-inf')

for idx, center in enumerate(centers):
    shortest_dist, farest_dist, angle_dist = distance_to_axes(center)
    if (shortest_dist + farest_dist < min_shortest_farest_sum) and (angle_dist > max_angle_dist):
        min_shortest_farest_sum = shortest_dist + farest_dist
        max_angle_dist = angle_dist
        best_cluster_idx = idx

# 获取目标类的数据
target_cluster_data = data[data['Cluster'] == best_cluster_idx]

# 保存目标类的数据
with pd.ExcelWriter(f'{file_path[:-5]}_Cluster_fit.xlsx') as writer:
    target_cluster_data.to_excel(writer, sheet_name=f'Cluster_{best_cluster_idx}', index=False)

print(f"符合条件的类为：Cluster {best_cluster_idx}，已保存目标类的数据。")


# 计算每类中离聚类中心最近和最远的5个点位置
closest_points = {}
farthest_points = {}

for cluster_idx in range(n_clusters):
    cluster_data = data[data['Cluster'] == cluster_idx]
    cluster_features = cluster_data.iloc[:, 4:7].values  # 提取Shortest, Farest, Angle列
    distances = cdist(cluster_features, [centers[cluster_idx]], 'euclidean').flatten()
    closest_indices = np.argsort(distances)[:5]
    farthest_indices = np.argsort(distances)[-5:]
    closest_points[cluster_idx] = cluster_data.iloc[closest_indices]
    farthest_points[cluster_idx] = cluster_data.iloc[farthest_indices]

# 保存每类中离聚类中心最近和最远的5个点位置
with pd.ExcelWriter(f'{file_path[:-5]}_Cluster_Center.xlsx') as writer:
    for cluster_idx in range(n_clusters):
        # 添加标签列
        closest_points[cluster_idx]['Label'] = 'Average'
        farthest_points[cluster_idx]['Label'] = 'Farthest'
        
        # 合并最近和最远的点
        combined_points = pd.concat([closest_points[cluster_idx], farthest_points[cluster_idx]])
        
        # 保存到同一个子表
        combined_points.to_excel(writer, sheet_name=f'Cluster_{cluster_idx}', index=False)

print("每类中离聚类中心最近和最远的5个点已保存。")


# 保存每个类的聚类中心信息
cluster_centers_df = pd.DataFrame(centers, columns=['Average', 'Farthest', 'Angle'])
cluster_centers_df['Cluster'] = range(n_clusters)
cluster_centers_df.to_excel('ClusterCenterInformation.xlsx', index=False)

print("每个类的聚类中心信息已保存至ClusterCenterInformation.xlsx。")