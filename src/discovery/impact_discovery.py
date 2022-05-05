from collections import defaultdict
import seaborn as sns
import pandas as pd
import json
import pickle

path_discovery_impact = '/Users/picho/projects/collaborations/bbglab/CH/ch-drivers/data/intogen/discovery_IMPACT/'

df_oncodr = pd.read_csv('{}/OTHER_WXS_CH_IMPACT_PANEL-oncodrivefml.tsv'.format(path_discovery_impact), 
                sep ='\t')
df_clustl = pd.read_csv('{}/oncodriveclustl_results.tsv'.format(path_discovery_impact), 
                sep ='\t')
df_dnds = pd.read_csv('{}/dndscsv_res.tsv'.format(path_discovery_impact), 
                sep ='\t')
df_smregions = pd.read_csv('{}/OTHER_WXS_CH_IMPACT_PANEL.smregions.tsv.gz'.format(path_discovery_impact), 
                sep ='\t')

list_smreg = df_smregions[df_smregions['Q_VALUE']<0.05]['HUGO_SYMBOL'].tolist()
list_clustl = df_clustl[df_clustl['Q_ANALYTICAL']<0.05]['SYMBOL'].tolist()
list_dndns = df_dnds[df_dnds['qglobal_cv']<0.05]['gene_name'].tolist()
list_fml = df_oncodr[df_oncodr['Q_VALUE']<0.05]['SYMBOL'].tolist()
total_disc = set(list_smreg + list_clustl+ list_dndns + list_fml)

dic_res = defaultdict(dict)
for g in total_disc:
    dic_res[g]['smRegions'] = 0
    dic_res[g]['dndsCV'] = 0
    dic_res[g]['oncodriveFML'] = 0
    dic_res[g]['oncodriveClustL'] = 0

    if g in list_smreg:
        dic_res[g]['smRegions'] = 1
    if g in list_dndns:
        dic_res[g]['dndsCV'] = 1
    if g in list_fml:
        dic_res[g]['oncodriveFML'] = 1
    if g in list_clustl:
        dic_res[g]['oncodriveClustL'] = 1


pickle.dump(dic_res, open("../../data/tables/driver_disc_impact_toweb.pckl", "wb"))

# filter genes by expression
df = pd.read_csv('~/projects/collaborations/bbglab/CH/mutpro/data/expression/GSE96811_RNAseq_counts.txt.gz', sep ='\t', 
                index_col = 0)

dic_count_fpkm = df.max(axis = 1).to_dict()
exp_g = [g for g,v in dic_count_fpkm.items() if v>15]
res_counts = pd.DataFrame(dic_res)
expressed_genes = [c for c in res_counts.columns.tolist() if c in exp_g]
res_counts = res_counts[expressed_genes]
total_disc = set(expressed_genes)
dic_res_counts = res_counts.sum().to_dict()

sorted_genes = sorted(dic_res_counts, key=dic_res_counts.get, reverse=True)

json.dump(list(res_counts.columns), open('../../data/tables/list_genes_impact.json', 'wt'))