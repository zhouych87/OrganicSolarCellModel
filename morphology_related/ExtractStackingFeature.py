# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 00:54:34 2024
1、提取结构中的特征，特征输出表中信息顺序为：
1_残基1编号    2_残基1内片段编号    3_残基2编号     4_残基2内片段编号     5_片间平均距离（2和4间）     
6_片间最小（2和4间）   7_片间最大（2和4间）    8_片间质心（2和4间）   9_分子间距离（1和3间）   10_角度 （2和4间）
11_残基1片内距离平均    12_残基1片内距离最小   13_残基1片内距离最大	   14_残基1片内中心距离最小（质心到原子）	  15_残基1片内中心距离最大（质心到原子）	   
16_残基2片内距离平均    17_残基2片内距离最小	18_残基2片内距离最大	   19_残基2片内中心距离最小（质心到原子）	  20_残基2片内中心距离最大（质心到原子）

2、自动获取残基名称和残基计数

3、 入口参数：.tpr文件和.trr文件。还需要要在c和d中输入对残基的切片信息，作为计算堆积的依据
    输出结果：自动获取20维度的特征，用于聚类和后期分析处理
    同时将输出的.txt堆积信息自动转换为.xlsx文件，便于查看
    
    能确保分子和片段的完整性，在角度映射时得到正确的结果。
    
创建于11月9日，修改于11月25日，修正了角度值计算的错误，加入了unwrap()在确保选择分子后得到完整的分子，再在分子中选择片段，
并进行平面映射得到角度。
@author: 123
"""

import MDAnalysis as mda
from MDAnalysis.analysis import distances
import math
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import cdist
import os
import matplotlib.pyplot as plt
from sklearn.cluster import BisectingKMeans, KMeans
from sklearn.preprocessing import MinMaxScaler
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn.utils.extmath import row_norms
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.cluster import BisectingKMeans, KMeans
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans
import math


topology_file = "eC9.tpr"
trajectory_file = "eC9.trr"

print("Loading molecular structure file...")
#u = mda.Universe("3PM6eC9confout.gro","traj.trr")  # 输出的gro和xtc文件
u = mda.Universe(topology_file, trajectory_file)  # 输出的gro和xtc文件

c = [
    "name   C38	C32	C28	C37	C43	N7	C44	N8	C29	C30	C31	O2	C36	C35	Cl1	C34	Cl	C33",
    "name   C11	C10	S1	C6	C5	C9	S3	C12	N3	C8	C2	C1	C39	N1	S5	N2	C40	C7	C4	C3	N4	C13	C14	S2	C15	C16	S4",
    "name   C27	C21	C17	C26	C41	N5	C42	N6	C18	C19	C20	O1	C25	C24	Cl2	C23	Cl3	C22"
]

d = [
    "name   C169	C170	C171	C172 	S23",
    "name   C166	C165	C164	C163	S21",
    "name   C128	C127	C126	C125	S16",
    "name   C98	C97	C96	C95	S13",
    "name   C60	C59	C58	C57	S8",
    "name   C38	C37	C36	C35	S5",
    "name   C176	C175	C174	C173	C168	S22	C167	C188	C187	C186	C185	S24	O5	O6",
    "name   C100	S14	C99	C116	C101	C102	C103	C104	S15	C113	C114	C115	O3	O4",
    "name   C40	S6	C39	C48	C41	C42	C43	C44	S7	C45	C46	C47	O1	O2",
    "name   C148	C147	C146	C133	C134	C135	C136	C137	S18	C132	C131	C150	C149	S19	C151	C152	C153	C154	S20	C130	C129	S17	F5	F6",
    "name   C80	C79	C78	C65	C66	C67	C68	C69	S10	C64	C63	C82	C81	S11	C83	C84	C85	C86	S12	C62	C61	S9	F3	F4",
    "name   C33	C32	C31	C18	C19	C20	C21	C22	S3	C15	C14	C17	C16	S2	C13	C34	S4	C12	C11	C10	C9	S1	F1	F2"
]

boxv = u.dimensions

#自动获取结构中残基名称和残基数量
def count_two_types_of_residues(topology_file, trajectory_file):
    """
    从给定的拓扑文件和轨迹文件中读取两种类型的残基名称及其个数。

    参数:
    topology_file (str): 拓扑文件路径（例如 "3PM6eC9.tpr"）
    trajectory_file (str): 轨迹文件路径（例如 "traj.trr"）

    返回:
    tuple: 一个包含四个值的元组，分别为第一种残基名称、第一种残基计数、第二种残基名称、第二种残基计数
    """
    # 创建 Universe 对象
    u = mda.Universe(topology_file, trajectory_file)

    # 获取所有残基对象
    residues = u.residues

    # 统计每个残基名称出现的次数，并记录每个残基第一次出现的位置
    residue_counts = {}
    first_appearance = {}
    for index, residue in enumerate(residues):
        residue_name = residue.resname
        if residue_name in residue_counts:
            residue_counts[residue_name] += 1
        else:
            residue_counts[residue_name] = 1
            first_appearance[residue_name] = index

    # 检查是否有两种类型的残基
    if len(residue_counts) != 2:
        raise ValueError("The structure file does not contain exactly two types of residues.")

    # 根据第一次出现的位置排序
    sorted_residue_types = sorted(residue_counts.items(), key=lambda x: first_appearance[x[0]])

    # 返回四个值
    residue_type_1_name, residue_type_1_count = sorted_residue_types[0]
    residue_type_2_name, residue_type_2_count = sorted_residue_types[1]

    return residue_type_1_name, residue_type_1_count, residue_type_2_name, residue_type_2_count

# 计算片段内原子间距离的平均值和最小值
def calculate_distance_AtomtoAtom(atoms, box=None):
    positions = np.array([atom.position for atom in atoms])
    if len(positions) < 2:
        return (0, 0, 0)  # 如果只有一个原子，则最大、最小和平均距离都为0    
    if box is None:
        distance_matrix = squareform(pdist(positions, metric='euclidean'))
    else:
        distance_matrix = distances.distance_array(positions, positions, box=box)    
    # 计算距离矩阵的平均值，距离矩阵的对角线元素为 0，因此我们忽略这些元素
    upper_triangle = distance_matrix[np.triu_indices_from(distance_matrix, k=1)]    
    if len(upper_triangle) == 0:
        return (0, 0, 0)  # 如果没有有效的距离（例如只有一个原子），返回0    
    distancesmax_chipin = np.max(upper_triangle)
    distancesmin_chipin = np.min(upper_triangle)
    distancemean_chipin = np.mean(upper_triangle)    
    return (distancesmax_chipin, distancesmin_chipin, distancemean_chipin)

#计算一个分子片段中最大三角形原子的位置索引
def getindx(crd):
    tdis1 = 0.0
    indx = [0, 0, 0]
    for j in range(len(crd)):
        for jj in range(len(crd)):
            if j != jj:
                dis1 = crd[j] - crd[jj]
                tdis = np.dot(dis1, dis1)
                if tdis > tdis1:
                    tdis1 = tdis
                    indx[0] = j
                    indx[1] = jj
    
    laxs = crd[indx[0]] - crd[indx[1]]
    tdis1 = 0.0
    
    for i in range(len(crd)):
        if i != indx[0] and i != indx[1]:
            tv0 = crd[i] - crd[indx[0]]
            theta = math.acos((np.dot(tv0, laxs)) / (np.linalg.norm(laxs) * np.linalg.norm(tv0)))
            tv0_norm = np.linalg.norm(tv0)
            tdis = tv0_norm * np.sin(theta)
            if tdis > tdis1:
                tdis1 = tdis
                indx[2] = i
                
    return indx

# 计算质心坐标的函数
def calculate_centroid(atoms):
    positions = np.array([atom.position for atom in atoms])  # 确保使用正确的属性
    centroid = np.mean(positions, axis=0)
    return centroid

# 计算片段质心到片段内部的最远和最近距离函数
def calculate_distance_CentoAtoms(centroid, atoms, box=None):
    positions = np.array([atom.position for atom in atoms])
    if len(positions) == 0:
        return (0, 0)  # 如果没有原子，返回0    
    centroid = np.array([centroid])    
    if box is None:
        distances_to_atoms = np.linalg.norm(centroid - positions, axis=1)
    else:
        distances_to_atoms = distances.distance_array(centroid, positions, box=box)[0]    
    distancesmax_CFar = np.max(distances_to_atoms)
    distancesmin_CNear = np.min(distances_to_atoms)    
    return (distancesmax_CFar, distancesmin_CNear)

#调用残基自动获取
residue_name_1, residue_count_1, residue_name_2, residue_count_2 = count_two_types_of_residues(topology_file, trajectory_file)
resA_count = residue_count_1
resD_count = residue_count_2
print("residue_name_1:", residue_name_1)
print("residue_name_2:", residue_name_2)
print("residue_count_1:", residue_count_1)
print("residue_count_2:", residue_count_2)

print("Selecting atom groups...")
g1_list = [u.select_atoms('resid ' + str(i) + ' and ' + c[j]) for i in range(1, resA_count + 1) for j in range(3)]
g2_list = [u.select_atoms('resid ' + str(ii) + ' and ' + c[jj]) for ii in range(1, resA_count + 1) for jj in range(3)]
g3_list = [u.select_atoms('resid ' + str(iii) + ' and ' + d[jjj]) for iii in range(resA_count + 1, resA_count + resD_count + 1) for jjj in range(12)]      
g4_list = [u.select_atoms('resid ' + str(iiii) + ' and ' + d[jjjj]) for iiii in range(resA_count + 1, resA_count + resD_count + 1) for jjjj in range(12)] 

satisfying_pairsAA = []
chipinpairs_aa_dis = []
chipinpairs_Ca_dis = []
pairs_AtoA = []
satisfying_pairsDD = []
satisfying_pairsAD = []


print("Calculating distance matrices and finding satisfying pairs AA...")
for i in range(1, resA_count + 1):
    for j in range(3):
        for ii in range(i + 1, resA_count + 1):
            for jj in range(3):
                g1_positions = g1_list[(i - 1) * 3 + j]
                g2_positions = g2_list[(ii - 1) * 3 + jj]
                
                # # 计算两个原子组之间的距离矩阵，考虑周期性边界条件
                dist_matrix = distances.distance_array(g1_positions, g2_positions, box=boxv)
                #计算片段间的平均距离和最大距离
                min_distance = np.min(dist_matrix)
                
                if min_distance < 5.0:
                    print(f"Found satisfying pair: Residue 1: {i}, Fragment: {j}, Residue 2: {ii}, Fragment: {jj}")
                    g3 = u.select_atoms('resid ' + str(i) + ' and all')
                    g4 = u.select_atoms('resid ' + str(ii) + ' and all')
                    g3.unwrap()
                    g4.unwrap()                    
                    g1 = g3.select_atoms(c[j])
                    g2 = g4.select_atoms(c[jj])
                    
                    print("Calculation  averag distancese chip in...")
                    distancesmax_chipin_1, distancesmin_chipin_1, distancemean_chipin_1 = calculate_distance_AtomtoAtom(g1, box=boxv)
                    distancesmax_chipin_2, distancesmin_chipin_2, distancemean_chipin_2 = calculate_distance_AtomtoAtom(g2, box=boxv)
                    print("g1:", g1)
                    print("g2:", g2)
                    print("g3:", g3)
                    print("g4:", g4)

                    
                    # 计算两个分子片段的质心
                    print("Calculation centroid position...")
                    centroid_1 = calculate_centroid(g1)
                    centroid_2 = calculate_centroid(g2)                    
                    # 计算每个质心到其分子片段内其他原子的最远距离
                    print("Calculation the distance from centroid position to the longest atom...")
                    distancesmax_CFar_1, distancesmin_CNear_1 = calculate_distance_CentoAtoms(centroid_1, g1, box=boxv)
                    distancesmax_CFar_2, distancesmin_CNear_2 = calculate_distance_CentoAtoms(centroid_2, g2, box=boxv)
                    
                    #计算两个片段的角度
                    print("Calculation the angle...")
                    va_indices = getindx(g1.positions)
                    vb_indices = getindx(g2.positions)
                    va_positions = [g1.positions[index] for index in va_indices]
                    vb_positions = [g2.positions[index] for index in vb_indices]
                    va1, va2, va3 = va_positions
                    vb1, vb2, vb3 = vb_positions
                    na = np.cross(va1 - va2, va3 - va2)
                    nb = np.cross(vb1 - vb2, vb3 - vb2)
                    cosine_angle = np.dot(na, nb) / (np.linalg.norm(na) * np.linalg.norm(nb))
                    angle_rad = np.arccos(cosine_angle)
                    degree = np.degrees(angle_rad)
                    degree = (degree if degree <= 90 else abs(degree - 180))
                    
                    #计算片段间中心距离和分子中心距离
                    print("Calculation distances fragment 1 and fragment 2, and the distance of the centroid between the two fragments center...")     
                    cmindis = distances.distance_array(g1.center_of_geometry(), g2.center_of_geometry(), box=boxv).min()
                    mcdistance = distances.distance_array(g3.center_of_geometry(), g4.center_of_geometry(), box=boxv).min()
                    
                    print("Calculation distances average between fragment 1 and fragment 2...")                    
                    avg_distance = np.mean(dist_matrix)
                    max_distance = np.max(dist_matrix)
                    
                    satisfying_pairsAA.append(
                        (i, j, ii, jj, avg_distance, min_distance, max_distance, cmindis, mcdistance, degree, distancemean_chipin_1, distancesmin_chipin_1, distancesmax_chipin_1, distancesmin_CNear_1, distancesmax_CFar_1, distancemean_chipin_2, distancesmin_chipin_2, distancesmax_chipin_2, distancesmin_CNear_2, distancesmax_CFar_2, g1_positions, g2_positions)
                    )
              
print("Saving satisfying pairs AA to files...")
with open('satisfying_pairsAA.txt', 'w') as output_file1:
    for i, j, ii, jj, avg_distance, min_distance, max_distance, cmindis, mcdistance, degree, distancemean_chipin_1, distancesmin_chipin_1, distancesmax_chipin_1, distancesmin_CNear_1, distancesmax_CFar_1, distancemean_chipin_2, distancesmin_chipin_2, distancesmax_chipin_2, distancesmin_CNear_2, distancesmax_CFar_2, g1_positions, g2_positions in satisfying_pairsAA:
        output_file1.write(f"{i}\t{j}\t{ii}\t{jj}\t{avg_distance:.6f}\t{min_distance:.6f}\t{max_distance:.6f}\t{cmindis:.6f}\t{mcdistance:.6f}\t{degree:.6f}\t{distancemean_chipin_1:.6f}\t{distancesmin_chipin_1:.6f}\t{distancesmax_chipin_1:.6f}\t{distancesmin_CNear_1:.6f}\t{distancesmax_CFar_1:.6f}\t{distancemean_chipin_2:.6f}\t{distancesmin_chipin_2:.6f}\t{distancesmax_chipin_2:.6f}\t{distancesmin_CNear_2:.6f}\t{distancesmax_CFar_2:.6f}\n")

# 读取 .txt 文件并转换为 DataFrame
df = pd.read_csv('satisfying_pairsAA.txt', sep='\t', header=None)
# 将 DataFrame 保存为 .xlsx 文件，并指定工作表名称
excel_file = f'{topology_file[:-4]}_AA.xlsx'
sheet_name = f'{topology_file[:-4]}_AA'
df.to_excel(excel_file, sheet_name=sheet_name, index=False, header=False)
print("Files have been saved successfully.")



print("Calculating distance matrices and finding satisfying pairs DD...")
for i in range(resA_count + 1, resA_count + resD_count + 1):       
    for j in range(12):         
        for ii in range(i + 1, resA_count + resD_count + 1):
            for jj in range(12):
                g1_positions = g3_list[(i - resA_count - 1) * 12 + j]   
                g2_positions = g4_list[(ii - resA_count - 1) * 12 + jj] 
                
                # # 计算两个原子组之间的距离矩阵，考虑周期性边界条件
                dist_matrix = distances.distance_array(g1_positions, g2_positions, box=boxv)
                #计算片段间的平均距离和最大距离
                min_distance = np.min(dist_matrix)
                
                if min_distance < 5.0:
                    print(f"Found satisfying pair: Residue 1: {i}, Fragment: {j}, Residue 2: {ii}, Fragment: {jj}")
                    g3 = u.select_atoms('resid ' + str(i) + ' and all')
                    g4 = u.select_atoms('resid ' + str(ii) + ' and all')
                    g3.unwrap()
                    g4.unwrap()                    
                    g1 = g3.select_atoms(d[j])
                    g2 = g4.select_atoms(d[jj])
                    
                    print("Calculation  averag distancese chip in...")
                    distancesmax_chipin_1, distancesmin_chipin_1, distancemean_chipin_1 = calculate_distance_AtomtoAtom(g1, box=boxv)
                    distancesmax_chipin_2, distancesmin_chipin_2, distancemean_chipin_2 = calculate_distance_AtomtoAtom(g2, box=boxv)
                    print("g1:", g1)
                    print("g2:", g2)
                    print("g3:", g3)
                    print("g4:", g4)

                    
                    # 计算两个分子片段的质心
                    print("Calculation centroid position...")
                    centroid_1 = calculate_centroid(g1)
                    centroid_2 = calculate_centroid(g2)                    
                    # 计算每个质心到其分子片段内其他原子的最远距离
                    print("Calculation the distance from centroid position to the longest atom...")
                    distancesmax_CFar_1, distancesmin_CNear_1 = calculate_distance_CentoAtoms(centroid_1, g1, box=boxv)
                    distancesmax_CFar_2, distancesmin_CNear_2 = calculate_distance_CentoAtoms(centroid_2, g2, box=boxv)
                    
                    #计算两个片段的角度
                    print("Calculation the angle...")
                    va_indices = getindx(g1.positions)
                    vb_indices = getindx(g2.positions)
                    va_positions = [g1.positions[index] for index in va_indices]
                    vb_positions = [g2.positions[index] for index in vb_indices]
                    va1, va2, va3 = va_positions
                    vb1, vb2, vb3 = vb_positions
                    na = np.cross(va1 - va2, va3 - va2)
                    nb = np.cross(vb1 - vb2, vb3 - vb2)
                    cosine_angle = np.dot(na, nb) / (np.linalg.norm(na) * np.linalg.norm(nb))
                    angle_rad = np.arccos(cosine_angle)
                    degree = np.degrees(angle_rad)
                    degree = (degree if degree <= 90 else abs(degree - 180))
                    
                    #计算片段间中心距离和分子中心距离
                    print("Calculation distances fragment 1 and fragment 2, and the distance of the centroid between the two fragments center...")     
                    cmindis = distances.distance_array(g1.center_of_geometry(), g2.center_of_geometry(), box=boxv).min()
                    mcdistance = distances.distance_array(g3.center_of_geometry(), g4.center_of_geometry(), box=boxv).min()
                    
                    print("Calculation distances average between fragment 1 and fragment 2...")                    
                    avg_distance = np.mean(dist_matrix)
                    max_distance = np.max(dist_matrix)
                    
                    satisfying_pairsDD.append(
                        (i, j, ii, jj, avg_distance, min_distance, max_distance, cmindis, mcdistance, degree, distancemean_chipin_1, distancesmin_chipin_1, distancesmax_chipin_1, distancesmin_CNear_1, distancesmax_CFar_1, distancemean_chipin_2, distancesmin_chipin_2, distancesmax_chipin_2, distancesmin_CNear_2, distancesmax_CFar_2, g1_positions, g2_positions)
                    )
                  
print("Saving satisfying pairs DD to files...")
with open('satisfying_pairsDD.txt', 'w') as output_file2:
    for i, j, ii, jj, avg_distance, min_distance, max_distance, cmindis, mcdistance, degree, distancemean_chipin_1, distancesmin_chipin_1, distancesmax_chipin_1, distancesmin_CNear_1, distancesmax_CFar_1, distancemean_chipin_2, distancesmin_chipin_2, distancesmax_chipin_2, distancesmin_CNear_2, distancesmax_CFar_2, g1_positions, g2_positions in satisfying_pairsDD:
        output_file2.write(f"{i}\t{j}\t{ii}\t{jj}\t{avg_distance:.6f}\t{min_distance:.6f}\t{max_distance:.6f}\t{cmindis:.6f}\t{mcdistance:.6f}\t{degree:.6f}\t{distancemean_chipin_1:.6f}\t{distancesmin_chipin_1:.6f}\t{distancesmax_chipin_1:.6f}\t{distancesmin_CNear_1:.6f}\t{distancesmax_CFar_1:.6f}\t{distancemean_chipin_2:.6f}\t{distancesmin_chipin_2:.6f}\t{distancesmax_chipin_2:.6f}\t{distancesmin_CNear_2:.6f}\t{distancesmax_CFar_2:.6f}\n")

# 读取 .txt 文件并转换为 DataFrame
df = pd.read_csv('satisfying_pairsDD.txt', sep='\t', header=None)
# 将 DataFrame 保存为 .xlsx 文件，并指定工作表名称
excel_file = f'{topology_file[:-4]}_DD.xlsx'
sheet_name = f'{topology_file[:-4]}_DD'
df.to_excel(excel_file, sheet_name=sheet_name, index=False, header=False)
print("Files have been saved successfully.")



print("Calculating distance matrices and finding satisfying pairs AD...")
for i in range(1, resA_count + 1):           
    for j in range(3):            
        for ii in range(resA_count + 1, resA_count + resD_count + 1):
            for jj in range(12):
                g1_positions = g1_list[(i - 1) * 3 + j]     
                g2_positions = g3_list[(ii - resA_count - 1) * 12 + jj]
                
                # # 计算两个原子组之间的距离矩阵，考虑周期性边界条件
                dist_matrix = distances.distance_array(g1_positions, g2_positions, box=boxv)
                #计算片段间的平均距离和最大距离
                min_distance = np.min(dist_matrix)
                
                if min_distance < 5.0:
                    print(f"Found satisfying pair: Residue 1: {i}, Fragment: {j}, Residue 2: {ii}, Fragment: {jj}")
                    g3 = u.select_atoms('resid ' + str(i) + ' and all')
                    g4 = u.select_atoms('resid ' + str(ii) + ' and all')
                    g3.unwrap()
                    g4.unwrap()                    
                    g1 = g3.select_atoms(c[j])
                    g2 = g4.select_atoms(d[jj])
                    
                    print("Calculation  averag distancese chip in...")
                    distancesmax_chipin_1, distancesmin_chipin_1, distancemean_chipin_1 = calculate_distance_AtomtoAtom(g1, box=boxv)
                    distancesmax_chipin_2, distancesmin_chipin_2, distancemean_chipin_2 = calculate_distance_AtomtoAtom(g2, box=boxv)
                    print("g1:", g1)
                    print("g2:", g2)
                    print("g3:", g3)
                    print("g4:", g4)
                    
                    # 计算两个分子片段的质心
                    print("Calculation centroid position...")
                    centroid_1 = calculate_centroid(g1)
                    centroid_2 = calculate_centroid(g2)                    
                    # 计算每个质心到其分子片段内其他原子的最远距离
                    print("Calculation the distance from centroid position to the longest atom...")
                    distancesmax_CFar_1, distancesmin_CNear_1 = calculate_distance_CentoAtoms(centroid_1, g1, box=boxv)
                    distancesmax_CFar_2, distancesmin_CNear_2 = calculate_distance_CentoAtoms(centroid_2, g2, box=boxv)
                    
                    #计算两个片段的角度
                    print("Calculation the angle...")
                    va_indices = getindx(g1.positions)
                    vb_indices = getindx(g2.positions)
                    va_positions = [g1.positions[index] for index in va_indices]
                    vb_positions = [g2.positions[index] for index in vb_indices]
                    va1, va2, va3 = va_positions
                    vb1, vb2, vb3 = vb_positions
                    na = np.cross(va1 - va2, va3 - va2)
                    nb = np.cross(vb1 - vb2, vb3 - vb2)
                    cosine_angle = np.dot(na, nb) / (np.linalg.norm(na) * np.linalg.norm(nb))
                    angle_rad = np.arccos(cosine_angle)
                    degree = np.degrees(angle_rad)
                    degree = (degree if degree <= 90 else abs(degree - 180))
                    
                    #计算片段间中心距离和分子中心距离
                    print("Calculation distances fragment 1 and fragment 2, and the distance of the centroid between the two fragments center...")     
                    cmindis = distances.distance_array(g1.center_of_geometry(), g2.center_of_geometry(), box=boxv).min()
                    mcdistance = distances.distance_array(g3.center_of_geometry(), g4.center_of_geometry(), box=boxv).min()
                    
                    print("Calculation distances average between fragment 1 and fragment 2...")                    
                    avg_distance = np.mean(dist_matrix)
                    max_distance = np.max(dist_matrix)
                    
                    satisfying_pairsAD.append(
                        (i, j, ii, jj, avg_distance, min_distance, max_distance, cmindis, mcdistance, degree, distancemean_chipin_1, distancesmin_chipin_1, distancesmax_chipin_1, distancesmin_CNear_1, distancesmax_CFar_1, distancemean_chipin_2, distancesmin_chipin_2, distancesmax_chipin_2, distancesmin_CNear_2, distancesmax_CFar_2, g1_positions, g2_positions)
                    )
                  
print("Saving satisfying pairs AD to files...")
with open('satisfying_pairsAD.txt', 'w') as output_file3:
    for i, j, ii, jj, avg_distance, min_distance, max_distance, cmindis, mcdistance, degree, distancemean_chipin_1, distancesmin_chipin_1, distancesmax_chipin_1, distancesmin_CNear_1, distancesmax_CFar_1, distancemean_chipin_2, distancesmin_chipin_2, distancesmax_chipin_2, distancesmin_CNear_2, distancesmax_CFar_2, g1_positions, g2_positions in satisfying_pairsAD:
        output_file3.write(f"{i}\t{j}\t{ii}\t{jj}\t{avg_distance:.6f}\t{min_distance:.6f}\t{max_distance:.6f}\t{cmindis:.6f}\t{mcdistance:.6f}\t{degree:.6f}\t{distancemean_chipin_1:.6f}\t{distancesmin_chipin_1:.6f}\t{distancesmax_chipin_1:.6f}\t{distancesmin_CNear_1:.6f}\t{distancesmax_CFar_1:.6f}\t{distancemean_chipin_2:.6f}\t{distancesmin_chipin_2:.6f}\t{distancesmax_chipin_2:.6f}\t{distancesmin_CNear_2:.6f}\t{distancesmax_CFar_2:.6f}\n")

# 读取 .txt 文件并转换为 DataFrame
df = pd.read_csv('satisfying_pairsAD.txt', sep='\t', header=None)
# 将 DataFrame 保存为 .xlsx 文件，并指定工作表名称
excel_file = f'{topology_file[:-4]}_AD.xlsx'
sheet_name = f'{topology_file[:-4]}_AD'
df.to_excel(excel_file, sheet_name=sheet_name, index=False, header=False)
print("Files have been saved successfully.")

   
                    
