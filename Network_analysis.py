# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 17:04:01 2017

@author: Chang-Eop
"""

import os
os.chdir('/Users/Chang-Eop/Desktop/GitHub/NetPharm')

import numpy as np
import pandas as pd
import networkx as nx
import pylab

ob_th = 30
dl_th = 0.18

H_name = pd.read_excel('02_Info_Herbs_Name.xlsx')
M_info = pd.read_excel('03_Info_Molecules.xlsx')
T_info = pd.read_excel('04_Info_Targets.xlsx')
H_M = pd.read_excel('23_Herbs_Molecules_Relationships.xlsx')
M_T = pd.read_excel('34_Molecules_Targets_Relationships.xlsx')   


#총 herb, compound 갯수 확인
print(H_M.herb_ID.unique().size, H_M.Mol_ID.unique().size)

# ADME filtering
M_info = M_info[M_info.ob > ob_th]
M_info = M_info[M_info.drug_likeness > dl_th]

#herbe-molecule 목록에서 ADME 필터링 한 moldecule만 남기기
H_M = H_M[H_M.Mol_ID.isin(M_info.MOL_ID)]  #H_M에 NaN이 몇개 있음 -_-;; 확인 필요


#중복 없는 herb 리스트
herbs_unique = H_M.herb_ID.unique()

#모든 본초를 key로, 각 본초의 networkx graph 객체(compound-target network)를 value로 갖는 dictionary
All_herbs_CT = {} 
for i in herbs_unique:
    mols = H_M.Mol_ID[H_M.herb_ID == i]
    G = nx.Graph()
    G.add_nodes_from(np.array(mols))
    for j in mols:
         targets = np.array(M_T.TARGET_ID[M_T.MOL_ID == j])
         G.add_nodes_from(targets)
         for k in targets:
             G.add_edge(j,k)
    All_herbs_CT['herb_ID_'+ str(i)] = G
    
#관심 처방의 본초 리스트
Formulae_list = ['herb_ID_4', 'herb_ID_23', 'herb_ID_221', 'herb_ID_273']

#관심 처방의 본초별 그래프 merge
Formulae_merged = All_herbs_CT[Formulae_list[0]]
for i in range(len(Formulae_list)):
    Formulae_merged = nx.compose(Formulae_merged, All_herbs_CT[Formulae_list[i]])
    
#node type에 따라 색깔 지정 list 생성
nodes_color = Formulae_merged.nodes()
for i in range(len(nodes_color)):
    if nodes_color[i][0:3] == 'MOL':
        nodes_color[i] = 'r'
    else:
        nodes_color[i] = 'b'
            
nx.draw_networkx(Formulae_merged,with_labels=True, node_color = nodes_color)
    


    