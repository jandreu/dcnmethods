#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 11:42:54 2022

Copyright (C) University of Essex

"""

## Extract dafa from SPM files
import scipy.io
import os
import numpy as np
from sklearn.svm import SVC
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import KFold
from itertools import combinations
from matplotlib import pyplot as plt
import pandas as pd
#from sklearn.model_selection import permutation_test_score
from sklearn.metrics import make_scorer, f1_score
from sklearn.inspection import permutation_importance
from sklearn_genetic import GASearchCV
from sklearn_genetic.space import Continuous, Categorical, Integer
from sklearn_genetic.callbacks import DeltaThreshold
from sklearn.feature_selection import RFECV
from itertools import compress
import pickle as pkl
import os
import sys
sys.path.append('./BPCA/')    
from pca_all_impute import PCAImputer  
from tools import permutation_test_score_alt, f_importances_coeff


permutations=100
population=10
generations=20
processes=-1

base = os.path.basename(sys.argv[0])
filename = os.path.splitext(base)
dirname = filename[0]

sep = os.sep
path =  '.' + sep + 'results_' + dirname
isExist = os.path.exists(path)
if not isExist:
    os.makedirs(path)
patho =  '.' + sep + 'outputs_' + dirname
isExist = os.path.exists(patho)
if not isExist:
    os.makedirs(patho)

dirfolder = '.' + sep + 'Input for MVPA'
subjdata = []
for file in os.listdir(dirfolder):
    d = os.path.join(dirfolder, file)
    if os.path.isdir(d):
        data = scipy.io.loadmat(d + os.sep +  'processed.mat')
        subjdata.append(data['d'])

##Find hemodynamic species included in data
hspecies = []
for subjd in subjdata:
    vnames = subjd.dtype.names
    hspecies = [x.replace('_betas','') for x in vnames if '_betas' in x]
    hspecies_vname =  [x for x in vnames if '_betas' in x]
    

## Find conditions for each subject and hemodynamic species
subj_cond = []        
for subjd in subjdata:
    cond_hsp = []     
    for i in range(0,len(hspecies)):
        design = subjd['design']
        cond_hsp.append([x[0] for x in design[0][0][0][i][0]['name'][0]])
    subj_cond.append(cond_hsp)  

## create analysis datasets
X = np.empty(len(hspecies_vname), dtype=object)
for i, _ in enumerate(X):
    X[i] = []
Y = np.empty(len(hspecies_vname), dtype=object)
for i, _ in enumerate(Y):
    Y[i] = []
    


channels =[]
for i in range(0,len(subjdata)):
    subjd = subjdata[i]
    ch = subjd['chlab']
    channels.append(list(ch[0][0][:,2]))
    hemo_conditions = subj_cond[i]
    for ii, hspecie in enumerate(hspecies_vname):
        beta = subjd[hspecie][0][0]
        conditions = hemo_conditions[ii]
        for iii, condition in enumerate(conditions):
            X[ii].append(beta[iii])
            Y[ii].append(condition)
            

subjdata = []  #clear subject data        
## convert to numpy array
for i in range(0,len(Y)):
    Y[i] = np.array(Y[i])
    X[i] = np.array(X[i])
 
##impute NaN by mean for each condition
imp = PCAImputer(method='bpca')

for idx, x in enumerate(X):
    y = Y[idx]
    for condition in conditions:
        xaux =  x[np.where(y==condition)]
        x_imp = imp.fit_transform(xaux)[0]
        x[np.where(y==condition)] = x_imp
    scaler = MinMaxScaler([0, 10])
    scaled_x = scaler.fit_transform(x)
    X[idx] = scaled_x

filename = path + sep + 'data_X.pkl'
with open(filename,'wb') as f:
    pkl.dump(X, f)
filename = path + sep + 'data_Y.pkl'
with open(filename,'wb') as f:
    pkl.dump(Y, f)

### define parameters search
param_grid = {'C': Continuous(1e-6, 1e+6, distribution='log-uniform'),
              'shrinking': Categorical([True,False])}

### tune optimizer
cv = KFold(n_splits=5, shuffle=True)

base_estimator = SVC(kernel='linear');

clf = GASearchCV(estimator=base_estimator,
                               cv=cv,
                               scoring='accuracy',
                               population_size=population,
                               generations=generations,
                               param_grid=param_grid,
                               n_jobs=processes,
                               error_score='raise',
                               keep_top_k=4)
callback = DeltaThreshold(threshold=0.01, generations=10, metric='fitness_max')

comb = set(list(combinations(conditions, 2)))
n_comb = sum(1 for _ in comb)
models = np.ndarray(shape=(len(X), n_comb), dtype=object)
dataset = np.ndarray(shape=(len(X), n_comb), dtype=object)
selected_names = np.ndarray(shape=(len(X), n_comb), dtype=object)

names = channels[0]
names.append('Score')
names.append('P-value')


filename = patho + sep + 'names.pkl'
with open(filename,'wb') as f:
    pkl.dump(names, f)
filename = patho + sep + 'hspecies.pkl'
with open(filename,'wb') as f:
    pkl.dump(hspecies, f)

## Run channel selection and optimization for each condition pair and hemo species
for i, _ in  enumerate(X):
    x = X[i]
    y = Y[i]
    comb = set(list(combinations(conditions, 2)))
    index_names = []
    for ii, comp in enumerate(comb):
        xt = x[np.where((y==comp[0]) | ( y==comp[1]))]
        yt = y[np.where((y==comp[0]) | ( y==comp[1]))]
        
        ## select channels
        selector = RFECV(base_estimator,cv=5,step=1, n_jobs=-1)
        xt = selector.fit_transform(xt,yt)
        support = selector.support_
        s_names = list(compress(channels[0],support))
        selected_names[i,ii] = s_names
        
        # optimize parameters
        clf.fit(xt, yt, callbacks=callback)
        models[i,ii] = [clf.best_estimator_,clf.best_params_]
        dataset[i,ii] = [xt,yt]
        index_names.append(comp[0] + '  vs ' + comp[1] )
        print("****** MODEL : " + str(i) + " AND " + str(ii) + " FINISHED ********")
        ## Save iter
        filename = patho + sep + 'models.pkl'
        with open(filename,'wb') as f:
            pkl.dump(models, f)
        filename = patho + sep + 'dataset.pkl'
        with open(filename,'wb') as f:
            pkl.dump(dataset, f)
        filename = patho + sep + 'index_names.pkl'
        with open(filename,'wb') as f:
            pkl.dump(index_names, f)
            

## permutation test to obtain signficances
for i, model_hemo in enumerate(models):
    df_ = []
    for iii in range(0,2):
        df_aux = pd.DataFrame(index=index_names, columns=names)
        df_.append(df_aux.fillna(0))

    for ii, model in enumerate(model_hemo):
         params= model[1]
         score_stats, coeff_stats = permutation_test_score_alt(base_estimator, dataset[i,ii][0], dataset[i,ii][1], fit_params=params, n_jobs=processes, n_permutations=permutations)
         filename = patho + sep + 'perm_scores_' + str(i) + '_' + str(ii) + '.pkl'
         with open(filename,'wb') as f:
             pkl.dump(score_stats, f)
         filename = patho + sep + 'perm_coeff_' + str(i) + '_' + str(ii) + '.pkl'
         with open(filename,'wb') as f:
             pkl.dump(coeff_stats, f)
         score = score_stats[0]
         pvalue = score_stats[2]
         stats = pd.Series([score,pvalue], index=['Score','P-value'])
         best_channels, all_channels = f_importances_coeff(coeff_stats, selected_names[i,ii])
         row = []
         row.append(pd.concat([best_channels,stats],axis=0))
         row.append(pd.concat([all_channels,stats],axis=0))
         for iii in range(0,len(row)):
             for index, value in row[iii].items():
                 df_[iii].loc[index_names[ii],index]= value
         print("****** PERMUTATION : " + str(i) + " AND " + str(ii) + " FINISHED ********")
         
    for iii in range(0,len(df_)):
        if iii==0:
            df_[iii].to_excel(path + sep  + hspecies[i] + "_significance.xlsx")
            df_[iii].to_csv(path + sep  + hspecies[i] + "_significance.csv")
        else:
            df_[iii].to_excel(path + sep  + hspecies[i] + "_value.xlsx")
            df_[iii].to_csv(path + sep  + hspecies[i] + "_value.csv")
 
     
     


## Build dataset for each hemoglobin species
#for hsp in hspecies:
      
         # sbjhb = subjd['%s_betas' % (hsp)][0][0]

