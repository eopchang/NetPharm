#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:27:52 2017
@author: Chang-Eop



H_list(herb_ID로 구성)의 본초들에 대해 target degree matrix 생성.
row: targets (TAR00001, TAR00002, ...)
col: herbs. (1,2,...(herb_ID))
value: degree 

"""

import os
os.chdir('/Users/Chang-Eop/Desktop/GitHub/NetPharm')



import numpy as np
import pandas as pd
import networkx as nx


H_name = pd.read_excel('02_Info_Herbs_Name.xlsx')
M_info = pd.read_excel('03_Info_Molecules.xlsx')
T_info = pd.read_excel('04_Info_Targets.xlsx')
D_info = pd.read_excel('05_Info_Diseases.xlsx')
H_M = pd.read_excel('23_Herbs_Molecules_Relationships.xlsx')
M_T = pd.read_excel('34_Molecules_Targets_Relationships.xlsx')   
T_D = pd.read_excel('45_Targets_Diseases_Relationships.xlsx')

T_info_toGene = pd.read_excel('04_Info_Targets_forMatching(curated).xlsx') #target name을 공식적인 gene name으로 변환위해 필요.

H_list = [335, 336]


ob_th = 30
dl_th = .18

#결과 3차원 어레이로     
N_Herb = H_name.shape[0]
N_Tar = T_info.shape[0] 

           
result_T = np.zeros((N_Tar, N_Herb))


for H in H_list: #[335]:#
    H_i = H - 1#herb_ID는 index+1임. H_i: index, H:herb_ID

    print('Herb_'+ str(H) + ' of '+str(N_Herb), ':')

    
    # ADME filtering
    M_info_th = M_info[M_info.ob > ob_th]
    M_info_th = M_info_th[M_info_th.drug_likeness > dl_th]
        
    #herbe-molecule 목록에서 ADME 필터링 한 moldecule만 남기기
    H_M_th = H_M[H_M.Mol_ID.isin(M_info_th.MOL_ID)]  #H_M에 NaN이 몇개 있음 -_-;; 확인 필요
    
    
    #중복 없는 herb 리스트
    herbs_unique = H_M_th.herb_ID.unique()
    
    ##C-T network construction #####'node type'을 attribute로 갖음.
    #해당 본초의 networkx graph 객체(compound-target network)를 value로 갖는 dictionary
    mols = H_M_th.Mol_ID[H_M_th.herb_ID == H]
    
    if len(mols) == 0:
        continue #threshold 넘는 컴파운드 없는 경우엔 다음 루프로.
    
    CT_network = nx.Graph()
    CT_network.add_nodes_from(np.array(mols))
    for j in mols:
        targets = np.array(M_T.TARGET_ID[M_T.MOL_ID == j])
        CT_network.add_nodes_from(targets)
        for k in targets:
            CT_network.add_edge(j,k)
    
    
    #CT_network의 degree 정보 추출
    # DataFrame 구조에서 ndarray구조 변환후 작업하고 다시 DataFrame으로 comeback(걍 그게 편해서)
    CT_degrees_items = CT_network.degree().items()
    CT_degrees_pd = pd.DataFrame(list(CT_degrees_items))
    CT_degrees_pd['type'] = 0
    CT_degrees_np = np.array(CT_degrees_pd)
    
    for i in CT_degrees_np[:,0]:
        if i[0] == 'M':
            CT_degrees_np[CT_degrees_np[:,0]==i,2] = 'Compound'
            CT_degrees_np[CT_degrees_np[:,0]==i,0] = list(M_info_th[M_info_th.MOL_ID == i].molecule_name)
        else:
            CT_degrees_np[CT_degrees_np[:,0]==i,2] = 'Target'
            CT_degrees_np[CT_degrees_np[:,0]==i,0] = list(T_info[T_info.TAR_ID == i].target_name)
    
    CT_degrees_pd = pd.DataFrame(CT_degrees_np, columns = ['Node', 'Degree','Type'])
    
    #degree descending 순으로 배열
    CT_degrees_pd = CT_degrees_pd.sort_index(by = 'Degree', ascending=False)
    
    #CT_degree 정보  T,M 으로 구분
    CT_degrees_T = CT_degrees_pd[CT_degrees_pd.Type =='Target']
    CT_degrees_M = CT_degrees_pd[CT_degrees_pd.Type =='Compound']
    
    #reuslt_Tmatrix에 degree 값 채워넣기.
    for i,j in enumerate(CT_degrees_T.Node):
        position = T_info.index[T_info.target_name == j][0]
        result_T[position, H_i] = np.array(CT_degrees_T)[i,1] #T_info.target_name 순서
        
        




