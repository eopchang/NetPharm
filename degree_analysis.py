# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 19:21:17 2017

@author: Chang-Eop
"""


#03_Info_Molecules.xlsx파일의 drug-likenss 컬럼명이 제대로 인식이 안되는 관계로  엑셀파일상에서 drug_likeness로 수정.
#34_Molecules_Targets_Relationships의 MOL ID 컬럼명 제대로 인식되도록 MOL_ID로 바꿈.
import os
os.chdir('/Users/Chang-Eop/Desktop/GitHub/NetPharm')

import numpy as np
import pandas as pd
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

#HMT matrix를 nbinary adjancy matrix로.
for i in range(no_H):
    print(i, 'of',no_H, '1st round!')
    #해당 herb의 moleucle 추출
    mols = H_M.Mol_ID[H_M.herb_ID == herbs_unique[i]]
    for j in mols:
         targets = np.array(M_T.TARGET_ID[M_T.MOL_ID == j])
         for k in targets:
             HMT_matrix[i,H_M.Mol_ID.unique() == j,int(k[3::])-1] = 1
             

#col= herbs, row = degree of compounds
#각 허브별로 row size 동일. 즉, 해당 본초에 포함되지 않은 compound도 빈칸으로 들어가있음. 
Degrees = np.zeros((no_M, no_H)) 
for i in range(no_H):
    print(i, 'of',no_H, '2nd round!!')
    Degrees[:,i] = np.sum(HMT_matrix[i,:,:],1)
    

#각 본초별로 degree distribution을 해당 폴더에 PNG 파일로 저장
# compound가 1개인 경우 histogram 그리기에 에러가 남 -> 재낌. 
for i in range(no_H):
    print(i, 'of',no_H, 'FIGURES SAVED!!!!!!!!!!!!!!!!!!!!!!!')
    DegDist = Degrees[:,i]
    DegDist = DegDist[DegDist !=0]
    if DegDist.size == 1:
        print('PASS!!!PASS!!!PASS!!!PASS!!!!PASS!!!!PASS!!!! ', i)
    else:
        pylab.hist(DegDist, bins=20)
        pylab.title(np.array(H_name.ix[[i],[0,2,3,4]])[0][1:])
        pylab.xlabel('Degrees'), pylab.ylabel('Frequency')
        pylab.savefig('/Users/Chang-Eop/Dropbox/가천대 Neural Network & Systems Medicine Lab/프로젝트/이원융_네트워크약리학리뷰/TCMSP_data/DegreeDists/'+str(np.array(H_name.ix[[i],[0,2,3,4]])[0][1:][0]))
        #pylab.show()
        pylab.close()
    

