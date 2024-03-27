#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 18:59:30 2023

@author: roshnishetty
"""

## Add this file to the notebooks directory in the original experimental patch-seq data repo: https://github.com/jcamunas/patchseq/tree/master
# Line 125 saves beta_glucose180.csv based on filter_condition2

import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import pickle
import numpy as np
%matplotlib inline
from singlet.dataset import Dataset
from singlet import SampleSheet, CountsTable
import matplotlib.patches as mpatches


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#%%


runfile('helper_functions.py')
runfile('predictions.py')

#%%

dict_phenotype_norm = {'CellSize_pF': 'Cell size',
                   'NormalizedTotalCapacitance_fF/pF': 'Total Exocitosis',
                   'NormalizedFirstDepolarizationCapacitance_fF/pF': 'Early exocytosis',
                   'NormalizedLateDepolarizationCapacitance': 'Late exocytosis',
                 'CalciumIntegralNormalizedtoCellSize_pC/pF': 'Ca2+ entry',
                   'CapacitanceNormalizedtoCalcium_fF/pC': 'Exocytosis norm Ca2+',
                   'NormalizedEarlyPeakCalciumCurrentAmplitude_pA/pF': 'Early Ca2+ current',
                 'NormalizedLateCalciumCurrentAmplitude_pA/pF': 'Late Ca2+ current',
                   'NormalizedLateCalciumChannelConductance_pS/pF' : 'Late Ca2+ Conductance',
                   'ReversalPotentialbyramp_mV': 'Reversal potential',
                   'NormalizedPeakSodiumCurrentAmplitude_pA/pF': 'Peak Na+ current',
                  'NormalizedSodiumChannelConductance_pS/pF': 'Na+ conductance'}

phenotype_norm = ['Cell size', 'Total Exocitosis','Early exocytosis','Late exocytosis',
                          'Ca2+ entry','Exocytosis norm Ca2+', 'Early Ca2+ current','Late Ca2+ current',
                          'Late Ca2+ Conductance','Reversal potential','Peak Na+ current','Na+ conductance']

#%%

root_folder = './../data/'
analysis_folder = './../analysis/'
resources_folder= './../resources/'
fig_folder = './../figures/fig2_suppfig3/'
correlations_folder = analysis_folder +'correlations_betacell/beta_correlations_all/beta_ND'
folder_gsea_sum = analysis_folder +'/GSEA_betacell/grouped'

#%%

#Load and reformat patch-seq dataset
filename =  root_folder + 'patchclamp_wcryo_human.counts.tab'
metadata = root_folder + 'patchclamp_wcryo_human.metadata.tab'

test = pd.read_csv(filename,sep='\t')
annotation = pd.read_csv(metadata,'\t')
annotation = annotation.drop(['cell_type'],axis=1)

annotation['DiabetesStatus'].replace({'heathy': 'healthy'}, inplace=True)
annotation['DiabetesStatus'].replace({'elevated HbA1c': 'T2D'}, inplace=True)
#add cell type information already computed
ct_all = pd.read_csv(analysis_folder + 'cell_typing_merged.csv', sep='\t', index_col=0, names=['cell_type'])
ct_FACS = pd.read_csv(analysis_folder +'cell_typing_FACS_endocrine.csv', sep='\t', index_col=0, names=['cell_type'])
ct_all = pd.concat([ct_all,ct_FACS])
annotation = annotation.join(ct_all['cell_type'], on='name')

annotation.rename(columns=dict_phenotype_norm, inplace=True)
#remove commas from data ephys
annotation[phenotype_norm] = annotation[phenotype_norm].apply(lambda x: pd.to_numeric(x.astype(str)
                                                   .str.replace(',',''), errors='coerce'))
annotation[phenotype_norm] = annotation[phenotype_norm].astype(float)

#remove genes not seen in 5 cells or having 10 counts total
test= filter_genes_pp(test, min_cells=5, min_counts=10)
#create dataset
ds = Dataset(counts_table=CountsTable(test),samplesheet=SampleSheet(annotation))
#remove nans in glucse
ds.samplesheet = ds.samplesheet[~ds.samplesheet[['Glucose_mM']].isnull().values]
#save unnormalized table
ds_pclamp_raw = ds.copy()
# Get data and normalize gene expression to combine cell size plot and marker genes
ds_norm = ds.copy()
ds_norm.counts = ds_norm.counts.normalize()
ds_norm.counts.pseudocount = 1
ds_norm.counts = ds_norm.counts.log(base=2)
ds_norm.counts.pseudocount = 1
ds_norm_pclamp = ds_norm.copy()

filter_condition = {'Cryopreserved': ['Yes']}
ds_t1d = filter_samplesheet(ds_norm, filter_dict= filter_condition)

filter_condition = {'Cryopreserved': ['No']}
ds_pclamp = filter_samplesheet(ds_norm, filter_dict= filter_condition)

#%%

filter_condition1 = {'cell_type': ['beta'], 
                    'DiabetesStatus': ['healthy'],
                   'Patched': ['Yes'],
                    'Glucose_mM': [5,10],
                    'TimefromDispersion_days': [1,2,3,4],
                              'preincubation': ['No', 'Yes']}

filter_condition2 = {'cell_type': ['beta'], 'Patched': ['Yes'], 'DiabetesStatus': ['healthy']}

filter_condition3 = {'cell_type': ['beta'], 'Patched': ['Yes'], 'Glucose_mM': [10], 
                    'DiabetesStatus': ['healthy']
                   }


filtered_dataset = filter_samplesheet(ds_pclamp, filter_dict= filter_condition2)
#filtered_dataset.samplesheet.to_csv('beta_glucose180.csv')

#%%
#where to output figures and plot names
folder_output_figs = '/Users/joan/Desktop/FIG2/'
plots_ID = 'beta_cell_main'

ephys_dict= {'Size': ['Cell size'],
         'Exocytosis': ['Total Exocitosis','Early exocytosis','Late exocytosis','Exocytosis norm Ca2+'],
         'Calcium channels': ['Ca2+ entry','Early Ca2+ current','Late Ca2+ current','Late Ca2+ Conductance'],
         'Sodium channels'  : ['Peak Na+ current','Na+ conductance']}
ephys_list = []
for key, par in ephys_dict.items():
    ephys_list = ephys_list + par
    
#Calculate correlation
df = filtered_dataset.samplesheet[ephys_list].copy()
#df = df[~(df==0)]
df = df.corr(method='spearman', min_periods=1) #calculate correlation matrix on samples
#keep absolute values and make a copy to label with non signed values
df_wabs = df.copy()
df = df.abs()

lut_gene = make_group_tags(df=df, dict_groups=ephys_dict)

#code to modify ordering of last clustering to make small plots consistent
from scipy.cluster import hierarchy
z = hierarchy.linkage(df, method='average', metric='correlation', optimal_ordering=False)
z2 = hierarchy.dendrogram(z)
reordered_index = z2['leaves']


ds_genetag = pd.DataFrame(index=df.index)
for group, genes in ephys_dict.items():
    for gene in genes:
        if gene in ds_genetag.index:
            ds_genetag.loc[gene,'group'] = group
ds_genetag = ds_genetag[ds_genetag['group'].notnull()]

#make legend
legend_plot = [mpatches.Patch(color=lut_gene[key], label=key) for key in lut_gene]
#add values onny on values higher than 0,3
annot= df_wabs[(df>0.3) & (df<1)].round(2).astype(str).replace({'nan':''})
annot = annot.iloc[reordered_index, reordered_index]

#plot
g =sns.clustermap(data=df,
                  row_linkage=z,
                  col_linkage=z,
               row_colors=ds_genetag['group'].map(lut_gene) , 
                   yticklabels=1,
               cmap='YlGnBu',
                 metric='correlation',
                 annot=annot,fmt = '',
                 xticklabels='',
                 cbar_kws={'label':'spearman r'})

#labels for gene groups
for key, color in lut_gene.items():
    g.ax_col_dendrogram.bar(0, 0, color=lut_gene[key],
                            label=key, linewidth=0)
    
g.ax_col_dendrogram.legend(loc="right", ncol=1, bbox_to_anchor=(1.3,0.6))
g.savefig(fig_folder +'Fig2A.pdf',dpi=300)