# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 20:58:17 2024

@author: 123
"""

import pandas as pd
import numpy as np

# 读取两个Excel文件
random_out = pd.read_excel('Randomout.xlsx', header=None)  # 没有表标题
aa_feature_scaled = pd.read_excel('eC9AD_Feature_Scaled.xlsx')

# 重命名Random_out的列，使其与AaAa_Feature_Scaled的前4列名称一致
random_out.columns = aa_feature_scaled.columns[:4]

# 使用merge方法找到匹配的行
scaled_feature_select = pd.merge(aa_feature_scaled, random_out, on=aa_feature_scaled.columns[:4].tolist(), how='inner')

# 保存匹配结果到新的Excel文件
scaled_feature_select.to_excel('Scaled_Feature_Select.xlsx', index=False)

# 复制文件并进行额外的计算
scaled_feature_select_distance = scaled_feature_select.copy()

# 对第5列、第6列和第7列分别做平方运算
scaled_feature_select_distance['Shortest^2'] = scaled_feature_select_distance.iloc[:, 4] ** 2
scaled_feature_select_distance['Farest^2'] = scaled_feature_select_distance.iloc[:, 5] ** 2
scaled_feature_select_distance['Angle^2'] = scaled_feature_select_distance.iloc[:, 6] ** 2
scaled_feature_select_distance['(1-Angle)^2'] = (1 - scaled_feature_select_distance.iloc[:, 6]) ** 2

# 对第10列数据实行“第10列数据减1再取绝对值”的运算
scaled_feature_select_distance['1-Angle^2'] = abs(scaled_feature_select_distance['Angle^2'] - 1)
#scaled_feature_select_distance['(1-Angle)^2'] = abs(scaled_feature_select_distance['Angle^2'] - 1)

# 计算第8列、第9列和第11列数据的平方和后开根号
scaled_feature_select_distance['Distance1'] = np.sqrt(scaled_feature_select_distance[['Shortest^2', 'Farest^2', '1-Angle^2']].sum(axis=1))
scaled_feature_select_distance['Distance_Couple2'] = np.sqrt(scaled_feature_select_distance[['Shortest^2', 'Farest^2', '(1-Angle)^2']].sum(axis=1))
scaled_feature_select_distance['Distance_Nominus'] = np.sqrt(scaled_feature_select_distance[['Shortest^2', 'Farest^2', 'Angle^2']].sum(axis=1))

# 保存最终结果到新的Excel文件
scaled_feature_select_distance.to_excel('Scaled_Feature_Select_Distance.xlsx', index=False)