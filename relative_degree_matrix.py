#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 20:19:22 2017

@author: Chang-Eop

WideRange_all_herbs.py로 생성된 result_T_all.npy,
 result_D_all.npy 파일이 폴더내에 존재해야 함.
figure fontsize는 적절히 조절해 사용할 것.
"""
import os
os.chdir('/Users/Chang-Eop/Desktop/GitHub/NetPharm')



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#관심 본초 herb_ID
Herb = [336]

#관심 본초가 타 본초 대비 상위 어느정도일때(0-1) visualize할것인지. (0.7=상위 30%)
R = 0.7



result_T_all = np.load('result_T_all.npy')
result_D_all = np.load('result_D_all.npy')


n = 10 #WideRange_all_herbs와 동일해야 함
d_th = 2#WideRange_all_herbs와 동일해야 함

T_info = pd.read_excel('04_Info_Targets.xlsx')
D_info = pd.read_excel('05_Info_Diseases.xlsx')
H_name = pd.read_excel('02_Info_Herbs_Name.xlsx')

N_Herb = H_name.shape[0]


#전체 본초들의 평균 결과
result_T_mean = np.mean(result_T_all, axis = 2)
result_D_mean = np.mean(result_D_all, axis = 2)

result_T_nonzero_mean = result_T_mean[np.sum(result_T_mean,1) > d_th*n,:]
result_D_nonzero_mean = result_D_mean[np.sum(result_D_mean,1) > d_th*n,:]

result_T_nonzero_name = T_info.target_name[np.sum(result_T_mean,1) > d_th*n]
result_D_nonzero_name = D_info.disease_name[np.sum(result_D_mean,1) > d_th*n]

plt.figure()
plt.title('Targets\n' +"Mean of " + str(N_Herb) + " Herbs")
plt.imshow(result_T_nonzero_mean)
plt.yticks(range(result_T_nonzero_name.size),result_T_nonzero_name, Fontsize = 4)
plt.xticks(range(n), range(n), Fontsize=4)



plt.figure()
plt.title('Diseases\n'+"Mean of " + str(N_Herb) + " Herbs")
plt.imshow(result_D_nonzero_mean)
plt.yticks(range(result_D_nonzero_name.size),result_D_nonzero_name, Fontsize = 4)
plt.xticks(range(n), range(n), Fontsize=4)







#######개별 본초들 확인
x = Herb[0] -1 #관심 본초. 

result_T = result_T_all[:,:,x]
result_D = result_D_all[:,:,x]
    
result_T_nonzero = result_T[np.sum(result_T,1) > d_th*n,:]
result_D_nonzero = result_D[np.sum(result_D,1) > d_th*n,:]

result_T_nonzero_name = T_info.target_name[np.sum(result_T,1) > d_th*n]
result_D_nonzero_name = D_info.disease_name[np.sum(result_D,1) > d_th*n]

plt.figure()
plt.title('Targets')
plt.imshow(result_T_nonzero)
plt.yticks(range(result_T_nonzero_name.size),result_T_nonzero_name, Fontsize = 4)
plt.xticks(range(n), range(n), Fontsize=4)


plt.figure()
plt.title('Diseases')
plt.imshow(result_D_nonzero)
plt.yticks(range(result_D_nonzero_name.size),result_D_nonzero_name, Fontsize = 4)
plt.xticks(range(n), range(n), Fontsize=4)



###########본초별로 오름차순 정리하고 관심본초당 각 행, 열에서의 비율 계산.
result_T_all_sorted = np.sort(result_T_all, axis = 2)
result_D_all_sorted = np.sort(result_D_all, axis = 2)

#먼저 본초별 각 행,열에서의 순위 역순 계산
result_T_rel = np.zeros(result_T.shape)
for i in range(result_T.shape[0]):
    for j in range(result_T.shape[1]):
        result_T_rel[i,j] = np.where(result_T_all_sorted[i,j] == result_T[i,j])[0][0]#동점시에 가장 작은 값

result_D_rel = np.zeros(result_D.shape)
for i in range(result_D.shape[0]):
    for j in range(result_D.shape[1]):
        result_D_rel[i,j] = np.where(result_D_all_sorted[i,j] == result_D[i,j])[0][0]#동점시에 가장 작은 값

#순위 역순을 전체 본초 갯수로 나눠서 0-1 사이의 비율로 변환        
result_T_rel = result_T_rel/(np.ones(np.shape(result_T_rel))*N_Herb)
result_D_rel = result_D_rel/(np.ones(np.shape(result_D_rel))*N_Herb)


####################################



result_T_rel_nonzero = result_T_rel[np.sum(result_T_rel,1) > R*n,:]
result_D_rel_nonzero = result_D_rel[np.sum(result_D_rel,1) > R*n,:]

result_T_nonzero_name = T_info.target_name[np.sum(result_T_rel,1) >R*n]
result_D_nonzero_name = D_info.disease_name[np.sum(result_D_rel,1) >R*n]

plt.figure()
plt.title('Targets\n'+ 'Relative Score')
plt.imshow(result_T_rel_nonzero)
plt.yticks(range(result_T_nonzero_name.size),result_T_nonzero_name, Fontsize = 4)
plt.xticks(range(n), range(n), Fontsize=4)


plt.figure()
plt.title('Diseases\n'+'Relative Score')
plt.imshow(result_D_rel_nonzero)
plt.yticks(range(result_D_nonzero_name.size),result_D_nonzero_name, Fontsize = 6)
plt.xticks(range(n), range(n), Fontsize=6)