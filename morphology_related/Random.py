# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 12:27:24 2024
在组数据中随机取出10组数
用于聚类之后的数据中随机取出5个分子对做量化计算，计算交叠积分
@author: 123
"""

import pandas as pd
import random


# 设置要抽取的行数
num_samples = 9  # 你可以根据需要更改这个值
# 定义随机种子
random_seed = 8

# 读取Excel文件
file_path = 'Random.xlsx'
df = pd.read_excel(file_path, header=None)

# 检查DataFrame是否为空或行数少于所需样本数
if len(df) < num_samples:
    print(f"警告: 文件中的行数少于{num_samples}行，将返回所有数据行。")
    sample_df = df
else:
    # 随机选择指定数量的行
    sample_df = df.sample(n=num_samples, random_state=random_seed)

# 打印结果
print(sample_df)

# 保存到新的Excel文件
output_file_path = f'{file_path[:-5]}_Randomout.xlsx'
sample_df.to_excel(output_file_path, index=False, header=None)
print(f"已将{num_samples}个随机样本保存到 {output_file_path}")

# 从输出文件中读取数据并提取第1列和第3列
grep_df = pd.read_excel(output_file_path, header=None, usecols=[0, 2])

# 保存到Grep.xlsx
grep_output_file_path = 'Grep.xlsx'
grep_df.to_excel(grep_output_file_path, index=False, header=None)
print(f"已将第1列和第3列的数据保存到 {grep_output_file_path}")