# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 17:04:01 2017

@author: Chang-Eop
"""

import os
os.chdir('/Users/Chang-Eop/Dropbox/가천대 Neural Network & Systems Medicine Lab/프로젝트/이원융_네트워크약리학리뷰/TCMSP_data')

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

#대상 herb, compound, target 갯수 확인
no_H = H_M.herb_ID.unique().size
no_M = H_M.Mol_ID.unique().size
no_T = M_T.TARGET_ID.unique().size
print(no_H, no_M, no_T)


HMT_matrix = np.zeros((no_H, no_M, int(max(M_T.TARGET_ID)[3::]))) 
#m-t 2D matrix가 no_H 개 깔림.
#herb는 ADME filtering으로 사라진 본초 빼고 나머지 갯수만큼만 생성
#molecule은 ADME filtering으로 걸러낸 나머지 갯수만큼만 생성
# target은 target ID의 일련번호 중 최댓값 크기로 생성!!!!!!!!!!!!!!!!

## 앞으로 분석에서 레퍼런스로 삼을 리스트 (인덱싱 0부터)
H_ref = H_name.herb_cn_name[H_M.herb_ID.unique()-1]
M_ref = pd.DataFrame(H_M.Mol_ID.unique(), index = range(no_M), columns = ['MOL_ID'])
##T-ref는 필요 없음. 인덱스 +1해서 TAR 일련번호 찾으면 됨(indexing 0부터이므)



#중복 없는 herb 리스트
herbs_unique = H_M.herb_ID.unique()

All_herbs_CT = {} #모든 본초를 key로, 각 본초의 networkx graph 객체(compound-target network)를 value로 갖는 dictionary
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
    
    
    
    
    


    