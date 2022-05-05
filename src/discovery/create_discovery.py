
import pandas as pd
from glob import glob
import os
import json
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from collections import defaultdict
import numpy as np

# config for matplotlib
def config_params(font_size=7):
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'

#------------------------------------------
# CREATE LIST DRIVER DISCOVERY WITH FILTERS
#------------------------------------------

# download this dataset from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96811 
expression_data = '/Users/picho/projects/collaborations/bbglab/CH/mutpro/data/expression/GSE96811_RNAseq_counts.txt.gz'
path_to_intogen_res = '../../data/intogen/'
path_to_CH_json = '../../data/tables/list_genes_CH_category.json'

# category CH
if not os.path.isfile(path_to_CH_json):
    print('create the list of genes in CH json by running  driver_CH_list.py')

dic_CH_genes = json.load(open(path_to_CH_json))

# result of discovery
all_f = {"OTHER_WGS_HMF_FULL":"{}/OTHER_WGS_HMF_FULL/all_drivers.tsv".format(path_to_intogen_res),
         "OTHER_WGS_HMF_MOSAIC":"{}/OTHER_WGS_HMF_MOSAIC/all_drivers.tsv".format(path_to_intogen_res),
         "OTHER_WGS_HMF_MUTECT":"{}/OTHER_WGS_HMF_MUTECT/all_drivers.tsv".format(path_to_intogen_res),
         "OTHER_WXS_TCGA_MOSAIC":"{}/OTHER_WXS_TCGA_MOSAIC/all_drivers.tsv".format(path_to_intogen_res),
         "OTHER_WXS_TCGA_FULL":"{}/OTHER_WXS_TCGA_FULL/all_drivers.tsv".format(path_to_intogen_res)
         }

# RNAseq counts
df = pd.read_csv(expression_data, sep ='\t', 
                index_col = 0)

dic_count_fpkm = df.max(axis = 1).to_dict()

toconcat = []
for k,f in all_f.items():
    df = pd.read_csv(f, sep ='\t')
    df['Expression_Healthy'] = df['SYMBOL'].map(dic_count_fpkm)
    df['COHORT'] = k
    toconcat.append(df)
    
merged = pd.concat(toconcat)

forbidden_genes = set()
for gene, data in merged.groupby(by='SYMBOL'):
    if len(data)==1:
        if 'FULL' in data['COHORT'].tolist()[0]:
            if len(data[data['CGC_GENE']==True])==0:
                forbidden_genes.add(gene)
                
merged = merged[(~merged['SYMBOL'].isin(forbidden_genes))]
germline_forbidden = merged[merged['FILTER']=='Germline Warning']['SYMBOL'].tolist()

set_preexpression = merged['SYMBOL'].unique()
expressed_removed = set(merged[(merged['Expression_Healthy']<=15)]['SYMBOL'].tolist())

list_to_remove_afterwards = merged[(merged['Expression_Healthy']>15) & 
(merged['FILTER'].isin(['Samples with more than 3 mutations', 'Known artifact', 
                      'Germline Warning', 'Warning Signature9']))]['SYMBOL'].tolist()

merged = merged[(merged['Expression_Healthy']>15)&
                (~merged['SYMBOL'].isin(forbidden_genes))&
                (~merged['SYMBOL'].isin(germline_forbidden))]

merged = merged[merged['FILTER'].isin(['PASS', 'Lack of literature evidence'])]

pass_filter = merged[(merged['FILTER']=='PASS')|((merged['FILTER']=='Lack of literature evidence') & (merged['SYMBOL'].isin(dic_CH_genes['CH'])))]
pass_filter.to_csv("../../data/tables/pass_filter.tsv", sep ="\t", index = False, header = True)

second_tier = merged[((merged['FILTER']=='Lack of literature evidence') & ~(merged['SYMBOL'].isin(dic_CH_genes['CH'])))]

pass_filter["DECISION"] = "IN COMPENDIUM"
second_tier["DECISION"] = "DISCARDED"
all_res = pd.concat([pass_filter, second_tier])
all_res["COHORT_NAME"] = all_res["COHORT"].apply(lambda x : x.split("_")[-2])
all_res["SET"] = all_res["COHORT"].apply(lambda x : x.split("_")[-1])
all_res[["COHORT_NAME", "SET", "SYMBOL", "DECISION"]].to_csv("../../data/tables/table2_supp.tsv", 
sep ="\t", index = False, header = True)

second_tier_genes = second_tier['SYMBOL'].unique()
list_drivers = set(pass_filter['SYMBOL'].tolist())

json.dump(list(list_drivers), open('../../data/tables/list_driver_discovery.json', 'wt'))
wanted_cols = ['SYMBOL', 'METHODS', 'QVALUE_COMBINATION', 'ROLE', 'CGC_GENE', 'COHORT']
pass_filter[wanted_cols].to_csv('../../data/tables/driver_tier1.tsv', sep ='\t', index = False, header = True)

#-------------------
#SUPPLEMENTARY PLOTS
#-------------------
# FIRST PLOT 
df = pd.read_csv(expression_data, sep ='\t', 
                index_col = 0)

good_genes = list(pass_filter['SYMBOL'].unique())

lgene = [s for s in set_preexpression if s in df.index.tolist()]
dic_exp = df.loc[lgene].T.max().to_dict()
sorted_g = sorted(dic_exp, key = dic_exp.get, reverse=True)

config_params(6)
fig, ax = plt.subplots(1, 1, figsize = (7, 2))
plt.yscale('log')
sns.boxplot(data = df.loc[sorted_g].T, color = '#FFE3B2', fliersize = 1,) 

plt.xticks(np.arange(len(lgene)), ["$\it{0}$".format(lab) if lab in good_genes 
                                 else  "#$\it{0}$".format(lab) if lab in expressed_removed 
                                   else "*$\it{0}$".format(lab) if lab in list_to_remove_afterwards 
                                   else  "^$\it{0}$".format(lab) if lab in second_tier_genes 
                                   else lab for lab in sorted_g], rotation = 90)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.ylabel('FPKM')
plt.xlabel('')
plt.hlines(15, 0, len(sorted_g), ls = '--', color = 'grey', lw = 0.25)
plt.savefig('../../data/supplementary/expr_genes.svg')
plt.close()

