#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 23:37:37 2017

@author: Chang-Eop
"""


#CT_network 그래프 객체는 CT_network
#TD_network 그래프 객체는 TD_network
import os
os.chdir('/Users/Chang-Eop/Desktop/GitHub/NetPharm')

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


#관심 처방의 본초 리스트
Formulae_list = [336] #herb ID 

H_name = pd.read_excel('02_Info_Herbs_Name.xlsx')
M_info = pd.read_excel('03_Info_Molecules.xlsx')
T_info = pd.read_excel('04_Info_Targets.xlsx')
D_info = pd.read_excel('05_Info_Diseases.xlsx')
H_M = pd.read_excel('23_Herbs_Molecules_Relationships.xlsx')
M_T = pd.read_excel('34_Molecules_Targets_Relationships.xlsx')   
T_D = pd.read_excel('45_Targets_Diseases_Relationships.xlsx')

T_info_toGene = pd.read_excel('04_Info_Targets_forMatching(curated).xlsx') #target name을 공식적인 gene name으로 변환위해 필요.

n = 10 #threshold bin size       

N_Tar = T_info.shape[0] 
N_Dis = D_info.shape[0]   

                
result_T = np.zeros((N_Tar, n))
result_D = np.zeros((N_Dis, n))                        
                             
            
#Compound filtering threshold 범위 설정 
max_ob = (max(M_info.ob)-min(M_info.ob))*0.2
max_dl = (max(M_info.drug_likeness)-min(M_info.drug_likeness))*0.2

ob_range = np.linspace(min(M_info.ob),max_ob,n)
dl_range = np.linspace(min(M_info.drug_likeness),max_dl,n)


for t in range(n):
    print(str(t+1) + ' of '+ str(n))
    ob_th = ob_range[t]
    dl_th = dl_range[t]

    # ADME filtering
    M_info = M_info[M_info.ob > ob_th]
    M_info = M_info[M_info.drug_likeness > dl_th]
    
    #herbe-molecule 목록에서 ADME 필터링 한 moldecule만 남기기
    H_M = H_M[H_M.Mol_ID.isin(M_info.MOL_ID)]  #H_M에 NaN이 몇개 있음 -_-;; 확인 필요
    
    
    #중복 없는 herb 리스트
    herbs_unique = H_M.herb_ID.unique()
    
    #모든 본초를 key로, 각 본초의 networkx graph 객체(compound-target network)를 value로 갖는 dictionary
    herbs_CT = {} 
    for i in Formulae_list:
        mols = H_M.Mol_ID[H_M.herb_ID == i]
        G = nx.Graph()
        G.add_nodes_from(np.array(mols))
        for j in mols:
             targets = np.array(M_T.TARGET_ID[M_T.MOL_ID == j])
             G.add_nodes_from(targets)
             for k in targets:
                 G.add_edge(j,k)
        herbs_CT['herb_ID_'+ str(i)] = G
        
    
        
    
    ##C-T network construction ########################### 
    #CT_network: graph 구조의 CT_netwowrk. 'node type'을 attribute로 갖음.
    #관심 처방의 본초별 그래프 merge
    CT_network = herbs_CT['herb_ID_'+ str(Formulae_list[0])]
    for i in range(len(Formulae_list)):
        CT_network = nx.compose(CT_network, herbs_CT['herb_ID_'+ str(Formulae_list[0])])
        
    
    
    
    #CT_network의 degree 정보 추출
    # DataFrame 구조에서 ndarray구조 변환후 작업하고 다시 DataFrame으로 comeback(걍 그게 편해서)
    CT_degrees_items = CT_network.degree().items()
    CT_degrees_pd = pd.DataFrame(list(CT_degrees_items))
    CT_degrees_pd['type'] = 0
    CT_degrees_np = np.array(CT_degrees_pd)
    
    for i in CT_degrees_np[:,0]:
        if i[0] == 'M':
            CT_degrees_np[CT_degrees_np[:,0]==i,2] = 'Compound'
            CT_degrees_np[CT_degrees_np[:,0]==i,0] = list(M_info[M_info.MOL_ID == i].molecule_name)
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
        result_T[position,t] = np.array(CT_degrees_T)[i,1] #T_info.target_name 순서
        
        

    
    
    
    
    ##T-D network construction#####################################
    #target 인수로 받아 해당 disease(list 구조) 출력하는 함수 정의.
    def disease(target):
        t_ID = T_D['target_ID']
        d_ID = T_D['disease_ID']
        result = np.array(d_ID[t_ID == target])
        return result
        
    #관심처방의 해당 타겟들 추출
    Targets = CT_network.nodes()
    for i in Targets:
        if i[0] == 'M':
            Targets = [x for x in Targets if x != i] #remove로 지우면 같은 분자 2개 이상일 경우 첫번째만 지우므로
    
    
    #########################특정 Targets 미리 구성했다면 여기서부터#######(disease 함수는 필요)
    #T-D network construction
    G = nx.Graph()
    G.add_nodes_from(Targets) #Targets: list
    for i in Targets:
        G.add_nodes_from(disease(i))
        for j in disease(i):
                 G.add_edge(i,j)
    TD_network = G
    
    
    
    
    #TD_network의 degree 정보 추출
    # DataFrame 구조에서 ndarray구조 변환후 작업하고 다시 DataFrame으로 comeback(걍 그게 편해서)
    TD_degrees_items = TD_network.degree().items()
    TD_degrees_pd = pd.DataFrame(list(TD_degrees_items))
    TD_degrees_pd['type'] = 0 
    TD_degrees_np = np.array(TD_degrees_pd)
    
    for i in TD_degrees_np[:,0]:
        if i[0] == 'D':
            TD_degrees_np[TD_degrees_np[:,0]==i,2] = 'Disease'
            TD_degrees_np[TD_degrees_np[:,0]==i,0] = list(D_info[D_info.DIS_ID == i].disease_name)
        else:
            TD_degrees_np[TD_degrees_np[:,0]==i,2] = 'Target'
            TD_degrees_np[TD_degrees_np[:,0]==i,0] = list(T_info[T_info.TAR_ID == i].target_name)
    
    TD_degrees_pd = pd.DataFrame(TD_degrees_np, columns = ['Node', 'Degree', 'Type'])
    
    #degree descending 순으로 배열
    TD_degrees_pd = TD_degrees_pd.sort_index(by = 'Degree', ascending=False)
    
    #CT_degree 정보 T, D로 구분
    TD_degrees_T = TD_degrees_pd[TD_degrees_pd.Type =='Target']
    TD_degrees_D = TD_degrees_pd[TD_degrees_pd.Type =='Disease']
    
    
    #reuslt_D matrix에 degree 값 채워넣기.
    for i,j in enumerate(TD_degrees_D.Node):
        position = D_info.index[D_info.disease_name == j][0]
        result_D[position,t] = np.array(TD_degrees_D)[i,1] #D_info.disease_name 순서
        
    
result_T_nonzero = result_T[np.sum(result_T,1) >0,:]
result_D_nonzero = result_D[np.sum(result_D,1) >0,:]

result_T_nonzero_name = T_info.target_name[np.sum(result_T,1) >0]
result_D_nonzero_name = D_info.disease_name[np.sum(result_D,1) >0]


