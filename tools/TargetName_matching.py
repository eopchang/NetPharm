#04_Info_Targets.xls와 drugbank_protein_identifiers.xlsx 이용. 

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 05:43:48 2017

@author: Chang-Eop
"""

import os
os.chdir('/Users/Chang-Eop/Desktop/GitHub/NetPharm')

import numpy as np
import pandas as pd
import re

pro_tcmsp = pd.read_excel('04_Info_Targets.xlsx')
drpi = pd.read_excel('drugbank_protein_identifiers.xlsx')

p = re.compile('\w', re.I) #대소문자 구별없이 문자열, 혹은 숫자 패턴

#protein name에서 대소문자, 숫자만 뽑음으로써 공백, '-'등의 표기방식 차이로 인한 문제를 제거
for i in range(pro_tcmsp.shape[0]):
    new_name = p.findall(pro_tcmsp.ix[i,'target_name'])
    new_name = "".join(new_name)
    pro_tcmsp.ix[i,'target_name'] = new_name
    
for i in range(drpi.shape[0]):
    new_name = p.findall(drpi.ix[i,'Name'])
    new_name = "".join(new_name)
    drpi.ix[i,'Name'] = new_name
                   
                   
        
pro_tcmsp['NewGeneName'] = np.arange(pro_tcmsp.shape[0])

count = 0
for i in range(pro_tcmsp.shape[0]):
    matched = drpi[drpi['Name'] == pro_tcmsp.ix[i,'target_name']]
    if matched.size == 0:
        pass
    else:
        if matched.ix[:,2].shape[0] > 1:
            pro_tcmsp.ix[i,'NewGeneName'] = str(list(matched.ix[:,2]))
            
            count = count +1
            #print(matched)
            
        else:
            names_string = list(matched.ix[:,2])[0]#스트링
            names_string = str(names_string) #nan이 있는 경우 에러 막기 위해
            pro_tcmsp.ix[i,'NewGeneName'] = names_string.partition(" ")[0]
            count = count + 1
        

pro_tcmsp.to_excel('04_Info_Targets_forMatching.xlsx') 
    
    
        

              
              