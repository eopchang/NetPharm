#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#networkx로 생성된 그래프 객체를 cytoscape에 import할 수 있는 엑셀파일로 저장해주는 모듈
"""
Created on Sun Mar 26 20:28:22 2017

@author: Neural Network & Sysyems Medicine Lab 
"""
    


import pandas as pd
import networkx as nx

def interactions(G, title):
    G.edges()
    G_dic_frame = pd.DataFrame(G.edges(), columns = ['node 1', 'node 2']) #단순 interaction 쌍 without attribute
    G_dic_frame.index.name = 'index' #cytoscape import시 에러막기 위해 이름 필요
    #연결 정보 저장
    G_dic_frame.to_excel(title + '.xlsx') 
    
   

def node_attr(G, title):
    #dictionary 구조에서 key순서 바뀌므로 반드시 sorted해야 함. 
    G_dic_frame_node = pd.DataFrame(sorted(G.nodes()), columns = ['node name'])
    G_dic_frame_node.index.name = 'index' #cytoscape import시 에러막기 위해 이름 필요
    attrs_node = list(list(G.node.values())[0])
    for i in attrs_node:
        node_attr = list(sorted(nx.get_node_attributes(G,i).values()))
        G_dic_frame_node[i+' (attr)'] = node_attr
        
    #node attribute 저장                 
    G_dic_frame_node.to_excel(title + '.xlsx')