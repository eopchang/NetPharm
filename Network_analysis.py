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

#관심 처방의 본초 리스트
Formulae_list = ['herb_ID_4', 'herb_ID_23', 'herb_ID_221', 'herb_ID_273']

H_name = pd.read_excel('02_Info_Herbs_Name.xlsx')
M_info = pd.read_excel('03_Info_Molecules.xlsx')
T_info = pd.read_excel('04_Info_Targets.xlsx')
D_info = pd.read_excel('05_Info_Diseases.xlsx')
H_M = pd.read_excel('23_Herbs_Molecules_Relationships.xlsx')
M_T = pd.read_excel('34_Molecules_Targets_Relationships.xlsx')   
T_D = pd.read_excel('45_Targets_Diseases_Relationships.xlsx')


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
    

    

##C-T network construction & visualize##    
#관심 처방의 본초별 그래프 merge
Formulae_merged = All_herbs_CT[Formulae_list[0]]
for i in range(len(Formulae_list)):
    Formulae_merged = nx.compose(Formulae_merged, All_herbs_CT[Formulae_list[i]])
    

#node type에 따라 색깔 지정 list 생성
nodes_color = Formulae_merged.nodes()
for i in range(len(nodes_color)):
    if nodes_color[i][0] == 'M':
        nodes_color[i] = 'r'
    else:
        nodes_color[i] = 'c'

##visualize        
pylab.figure()  

#position of nodes  
pos=nx.spring_layout(Formulae_merged)

#position of labels
pos_labels = {}
y_off = .01 #상황에 맞게 조정
for i, v in pos.items():
    pos_labels[i] = (v[0], v[1]+y_off)
    
nx.draw_networkx(Formulae_merged, pos, with_labels=False, node_color = nodes_color)   
nx.draw_networkx_labels(Formulae_merged, pos_labels)


##T-D network construction & visualize
#target 인수로 받아 해당 disease(list 구조) 출력하는 함수 정의.
def disease(target):
    t_ID = T_D['target_ID']
    d_ID = T_D['disease_ID']
    result = np.array(d_ID[t_ID == target])
    return result
    
#관심처방의 해당 타겟들 추출
Targets = Formulae_merged.nodes()
for i in Targets:
    if i[0] == 'M':
        Targets = [x for x in Targets if x != i] #remove로 지우면 같은 분자 2개 이상일 경우 첫번째만 지우므로


#########################특정 Targets 미리 구성했다면 여기서부터#######
#T-D network construction
G = nx.Graph()
G.add_nodes_from(Targets) #Targets: list
for i in Targets:
    G.add_nodes_from(disease(i))
    for j in disease(i):
             G.add_edge(i,j)
TD_network = G
    

#node type에 따라 색깔 지정 list 생성
nodes_color2 = TD_network.nodes()
for i in range(len(nodes_color2)):
    if nodes_color2[i][0] == 'D':
        nodes_color2[i] = 'g'
    else:
        nodes_color2[i] = 'c'


##visualize
#시각화 위해 노드 번호 int로 바꾼 새로운 그래프(TD_network_int) 만듦.      
TD_network_int = nx.convert_node_labels_to_integers(TD_network)
pylab.figure()  

#position of nodes  
pos=nx.spring_layout(TD_network_int) #position of nodes

#position of labels
pos_labels = {}
y_off = 0 #상황에 맞게 조정
for i, v in pos.items():
    pos_labels[i] = (v[0], v[1]+y_off)
    
nx.draw_networkx(TD_network_int,pos, with_labels=False, node_color = nodes_color2)   
nx.draw_networkx_labels(TD_network_int,pos_labels, font_size = 8)


#T: int형으로 바뀐 node 들의 실제 이름 
isD = [i[0] == 'D' for i in TD_network.nodes()] #node가 disease인지(target인지) 여부
T = TD_network_int.nodes()
for i in T:
    if isD[i]: 
        T[i] = list(D_info[D_info['DIS_ID'] == TD_network.nodes()[i]]['disease_name'])
    else:
        T[i] = list(T_info[T_info['TAR_ID'] == TD_network.nodes()[i]]['target_name'])
        
T = pd.DataFrame(T)

# int로 라벨링 된 노드들의 실제 이름 엑셀파일로 저장
T.to_excel('TD_index.xlsx')

from collections import defaultdict
DT_table =  defaultdict(list)#Disease와 해당 Target 들 테이블로 정리
for i in TD_network_int.edges():
    if isD[i[0]]:
        DisName = T[0][i[0]] #frame구조에서 series 뽑은 다음(T[0]) index(i[0])
        TarName = T[0][i[1]]
        DT_table[DisName].append(TarName)
    else:
        DisName = T[0][i[1]] #frame구조에서 series 뽑은 다음(T[0]) index(i[0])
        TarName = T[0][i[0]]
        DT_table[DisName].append(TarName)
        
DT_table_frame = pd.DataFrame(list(DT_table.items()))
DT_table_frame.to_csv('DT_table') #value에 list 구조때문에 excel로는 저장이 안됨.