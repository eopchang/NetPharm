#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 22:47:18 2017

@author: Chang-Eop
"""
#간소화 버전
#Output:
    #1. CT_network 그래프 객체
    #2. TD_network 그래프 객체
    #3. CT_network의 두 종류 노드(compound, target)별 degree 정리된 엑셀파일
    #4. TD_network의 두 종류 노드(Target, Disease)별 degree 정리된 엑셀파일

#CT_network 그래프 객체는 CT_network
#TD_network 그래프 객체는 TD_network
import os
os.chdir('/Users/Chang-Eop/Desktop/GitHub/NetPharm')

import numpy as np
import pandas as pd
import networkx as nx


#Compound filtering threshold 설정 
ob_th = 30 #Oral bioavility(default = 30)
dl_th = 0.18 #drug-likeness(default = 0.18)


title_CT_degrees_T = 'CT_degrees_T_ginseng.xlsx' #compound-target network의 target node degrees
title_CT_degrees_M = 'CT_degrees_M_ginseng.xlsx' #compound-target network의 compound node degrees
title_TD_degrees_T = 'TD_degrees_T_ginseng.xlsx' #target-disease network의 target node degrees
title_TD_degrees_D = 'TD_degrees_D_ginseng.xlsx' #target-disease network의 disease node degrees

#관심 처방의 본초 리스트
Formulae_list = [336]

H_name = pd.read_excel('02_Info_Herbs_Name.xlsx')
M_info = pd.read_excel('03_Info_Molecules.xlsx')
T_info = pd.read_excel('04_Info_Targets.xlsx')
D_info = pd.read_excel('05_Info_Diseases.xlsx')
H_M = pd.read_excel('23_Herbs_Molecules_Relationships.xlsx')
M_T = pd.read_excel('34_Molecules_Targets_Relationships.xlsx')   
T_D = pd.read_excel('45_Targets_Diseases_Relationships.xlsx')

T_info_toGene = pd.read_excel('04_Info_Targets_forMatching(curated).xlsx') #target name을 공식적인 gene name으로 변환위해 필요.



# ADME filtering
M_info = M_info[M_info.ob > ob_th]
M_info = M_info[M_info.drug_likeness > dl_th]

#herbe-molecule 목록에서 ADME 필터링 한 moldecule만 남기기
H_M = H_M[H_M.Mol_ID.isin(M_info.MOL_ID)]  #H_M에 NaN이 몇개 있음 -_-;; 확인 필요



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
    herbs_CT[i] = G
    

    

##C-T network construction ########################### 
#CT_network: graph 구조의 CT_netwowrk. 'node type'을 attribute로 갖음.
#관심 처방의 본초별 그래프 merge
CT_network = herbs_CT[Formulae_list[0]]
for i in range(len(Formulae_list)):
    CT_network = nx.compose(CT_network, herbs_CT[Formulae_list[i]])
    



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

#CT_degree 정보 excel로 저장
CT_degrees_T = CT_degrees_pd[CT_degrees_pd.Type =='Target']
CT_degrees_T.to_excel(title_CT_degrees_T)

CT_degrees_M = CT_degrees_pd[CT_degrees_pd.Type =='Compound']
CT_degrees_M.to_excel(title_CT_degrees_M)




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

#CT_degree 정보 excel로 저장
TD_degrees_T = TD_degrees_pd[TD_degrees_pd.Type =='Target']
TD_degrees_T.to_excel(title_TD_degrees_T)

TD_degrees_D = TD_degrees_pd[TD_degrees_pd.Type =='Disease']
TD_degrees_D.to_excel(title_TD_degrees_D)





#이하 GSEA 분석. 작성중.

import gseapy as gp

glist = [str(list(T_info_toGene[T_info_toGene.TAR_ID == x].gene_name)[0]) for x in T_info_toGene.TAR_ID]
#gene name이 int인 경우 있어 str로 변환.

glist = [x.upper() for x in glist]
#gseapy사용법 따라 대문자로 변환

enrichr_results = gp.enrichr(gene_list=glist, description='test_name', gene_sets= 'Human_Phenotype_Ontology',
                             outdir='enrichr_test', cutoff=0.5, scale=0.8, no_plot=True)

results_sig = enrichr_results[enrichr_results.ix[:,3] < .05].Term
print(results_sig.shape)                             


#'GO_Molecular_Function_2015' 
#'PPI_Hub_Proteins'
#'OMIM_Disease'
#'WikiPathways_2016'
#'Human_Phenotype_Ontology'
#'GO_Biological_Process_2015'



