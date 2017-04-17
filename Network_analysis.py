# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 17:04:01 2017

@author: Chang-Eop
"""
#CT_network 그래프 객체는 CT_network
#TD_network 그래프 객체는 TD_network
import os
os.chdir('/Users/Chang-Eop/Desktop/GitHub/NetPharm')

import numpy as np
import pandas as pd
import networkx as nx
import pylab, copy

#Compound filtering threshold 설정 
ob_th = 30 #Oral bioavility(default = 30)
dl_th = 0.18 #drug-likeness(default = 0.18)

title_CT_index = 'CT_index_ginseng.xlsx'
title_TD_index = 'TD_ginseng.xlsx'
title_DT_table = 'DT_table_ginseng' #주의: CSV 포맷으로 저장할거임
title_CT_degrees_T = 'CT_degrees_T_ginseng.xlsx' #compound-target network의 target node degrees
title_CT_degrees_M = 'CT_degrees_M_ginseng.xlsx' #compound-target network의 compound node degrees
title_TD_degrees_T = 'TD_degrees_T_ginseng.xlsx' #target-disease network의 target node degrees
title_TD_degrees_D = 'TD_degrees_D_ginseng.xlsx' #target-disease network의 disease node degrees

#관심 처방의 본초 리스트
Formulae_list = ['herb_ID_336']

H_name = pd.read_excel('02_Info_Herbs_Name.xlsx')
M_info = pd.read_excel('03_Info_Molecules.xlsx')
T_info = pd.read_excel('04_Info_Targets.xlsx')
D_info = pd.read_excel('05_Info_Diseases.xlsx')
H_M = pd.read_excel('23_Herbs_Molecules_Relationships.xlsx')
M_T = pd.read_excel('34_Molecules_Targets_Relationships.xlsx')   
T_D = pd.read_excel('45_Targets_Diseases_Relationships.xlsx')

T_info_toGene = pd.read_excel('04_Info_Targets_forMatching(curated).xlsx') #target name을 공식적인 gene name으로 변환위해 필요.



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
#Formulae_merged: graph 구조의 CT_netwowrk. 'node type'을 attribute로 갖음.
#관심 처방의 본초별 그래프 merge
Formulae_merged = All_herbs_CT[Formulae_list[0]]
for i in range(len(Formulae_list)):
    Formulae_merged = nx.compose(Formulae_merged, All_herbs_CT[Formulae_list[i]])
    



#node type에 따라 색깔 지정 list 생성 & nodes type에 따라 node attribute 생성
nodes = Formulae_merged.nodes() #for indexing
nodes_color = Formulae_merged.nodes() #for color
a = dict.fromkeys(Formulae_merged.nodes()) #for attribute

for i in range(len(nodes)):
    if nodes[i][0] == 'M':
        a[nodes[i]] = 'Mol'
        nodes_color[i] = 'r'
        
    else:
        a[nodes[i]] = 'Tar'
        nodes_color[i] = 'c'
        
#CT_network에 'node type'attribute 부여
nx.set_node_attributes(Formulae_merged, 'node type', a)

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

CT_nodes = []
for i in Formulae_merged.nodes():
    if i[0] == 'T':
        ts = list(T_info[T_info['TAR_ID'] == i]['target_name'])
        CT_nodes.append(ts)
    else:
        ms = list(M_info[M_info['MOL_ID'] == i]['molecule_name'])
        CT_nodes.append(ms)
    
CT_index = pd.DataFrame(Formulae_merged.nodes())
CT_index[1] = CT_nodes
CT_index.to_excel(title_CT_index)



CT_network = Formulae_merged

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

# int로 바꾸면서 networkx가 node 순서 마음대로 바꿔버리므로 그 순서에 맞게 node color도 재배열해야 함. 
nodes_color_rearanged = np.array(nodes_color2)[np.array(TD_network_int.nodes())]

nx.draw_networkx(TD_network_int,pos, with_labels=False, node_color = nodes_color_rearanged)   
nx.draw_networkx_labels(TD_network_int,pos_labels, font_size = 8)


#T: int형으로 바뀐 node 들의 실제 이름 
isD = [i[0] == 'D' for i in TD_network.nodes()] #node가 disease인지(target인지) 여부
T = TD_network_int.nodes()
T.sort() #int로 바꾸면서 노드 순서 뒤바꼈으므로 다시 배열하여 노드 리스트 작성하고자 함.
for i in T:
    if isD[i]: 
        T[i] = list(D_info[D_info['DIS_ID'] == TD_network.nodes()[i]]['disease_name'])
    else:
        T[i] = list(T_info[T_info['TAR_ID'] == TD_network.nodes()[i]]['target_name'])
        
T = pd.DataFrame(T)

# int로 라벨링 된 노드들의 실제 이름 엑셀파일로 저장
T.to_excel(title_TD_index)

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
DT_table_frame.to_csv(title_DT_table) #value에 list 구조때문에 excel로는 저장이 안됨.



CT_network.nodes()    
list(T_info_toGene[T_info_toGene.TAR_ID == 'TAR00087'].gene_name)                 
[list(T_info_toGene[T_info_toGene.TAR_ID == x].gene_name)[0] for x in T_info_toGene.TAR_ID]




#이하 GSEA 분석. 작성중.

import gseapy as gp

glist = [str(list(T_info_toGene[T_info_toGene.TAR_ID == x].gene_name)[0]) for x in T_info_toGene.TAR_ID]
#gene name이 int인 경우 있어 str로 변환.

glist = [x.upper() for x in glist]
#gseapy사용법 따라 대문자로 변환

enrichr_results = gp.enrichr(gene_list=glist, description='test_name', gene_sets= 'Human_Phenotype_Ontology',
                             outdir='enrichr_kegg', cutoff=0.5, scale=0.8, no_plot=True)

results_sig = enrichr_results[enrichr_results.ix[:,3] < .05].Term
print(results_sig.shape)                             


#'GO_Molecular_Function_2015' 
#'PPI_Hub_Proteins'
#'OMIM_Disease'
#'WikiPathways_2016'
#'Human_Phenotype_Ontology'
#'GO_Biological_Process_2015'

