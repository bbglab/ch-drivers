import pandas as pd
import json
import numpy as np
from pyliftover import LiftOver
from collections import defaultdict
from tqdm import tqdm
import pickle 

tqdm.pandas() 

def do_liftover_test(df, lo):
    try:
        lov = lo.convert_coordinate('chr'+df['CHROM'], int(df['POS']))[0]
        df['CHROMOSOME'] = lov[0]
        df['POSITION'] = int(lov[1])
    except:
        df['CHROMOSOME'] = None
        df['POSITION'] = None
    return df

def do_liftover_test_impact(df, lo):
    
    try:
        lov = lo.convert_coordinate('chr'+df['Chromosome'], int(df['Start_Position']))[0]
        df['CHROMOSOME'] = lov[0]
        df['POSITION'] = round(int(lov[1]))
    except:
        df['CHROMOSOME'] = None
        df['POSITION'] = None

    return df

# -----
# Path to files
# -----
#  
outpath = "../../data/tables/"

# update with latest run 
hmf_vep = "/Users/picho/projects/collaborations/bbglab/CH/clonal_hematopoiesis_code/data_VAF/HMF/vaf_filter/HMF_FULL/HMF_full.vep.sorted.tsv.filtered.toUSE.gz"
tcga_vep = "/Users/picho/projects/collaborations/bbglab/CH/clonal_hematopoiesis_code/data_VAF/TCGA/vaf_filter/TCGA_FULL/TCGA_full.vep.sorted.tsv.filtered.toUSE.gz"

#clonal_hematopoiesis_run_20210410
path_intogen_run = "~/projects/collaborations/bbglab/CH/ch-drivers/data/intogen/merged_run/"

cgc_parsed = "../../data/misc/cancer_gene_census_parsed.tsv"

category_CH = "../../data/tables/list_genes_CH_category.json"
list_discovery = '../../data/tables/list_driver_discovery.json'
list_impact = "../../data/tables/list_genes_impact.json"

expression_data = '/Users/picho/projects/collaborations/bbglab/CH/mutpro/data/expression/GSE96811_RNAseq_counts.txt.gz'

biomart_cds = "../../data/misc/cds_biomart.tsv"

# the TCGA clinical table can be downloaded from their website
clincal_data_tcga = '/Users/picho/projects/collaborations/bbglab/CH/cluster/clinical_PANCAN_patient_with_followup_151210.tsv'


ch_cbio = '~/projects/collaborations/bbglab/CH/mutpro/data/other_CH_datasets/data_mutations_extended.txt'
ch_into = "~/projects/collaborations/bbglab/CH/mutpro/data/other_CH_datasets/mutations.tsv"
#---------------------------------------------------------------------------------------------------------------------
# 1 - Creating the mutational table
#---------------------------------------------------------------------------------------------------------------------

#---------------------------------------
# HARTWIG has to be liftovered to grch38
#---------------------------------------
hartwig = pd.read_csv(hmf_vep, sep ='\t')
hartwig['CHROM'] = hartwig['ID'].apply(lambda x : x.split('_')[0])
hartwig['POS'] = hartwig['ID'].apply(lambda x : int(x.split('_')[1]))

lo = LiftOver('hg19', 'hg38')
df2 = hartwig.progress_apply(do_liftover_test, args = (lo, ), axis = 1)
df2 = df2[~(df2['CHROMOSOME'].isnull())]

dic_hmf_res = defaultdict(dict)

for ix, (idv, data) in enumerate(df2.groupby(by='ID')):

    dic_hmf_res[ix]['MUTATION'] = '{}:{}:{}>{}'.format(data['CHROMOSOME'].tolist()[0].replace('chr', ''),
                                                   int(data['POSITION'].tolist()[0]),  
                                                  data['REF'].tolist()[0],
                                                  data['ALT'].tolist()[0],)
    dic_hmf_res[ix]['COHORT'] = 'OTHER_WGS_HMF_FULL'
    dic_hmf_res[ix]['CONSEQUENCE'] = data['Consequence'].tolist()[0]
    dic_hmf_res[ix]['CHR'] = data['CHROMOSOME'].tolist()[0].replace('chr', '')
    dic_hmf_res[ix]['POS'] = int(data['POSITION'].tolist()[0])
    dic_hmf_res[ix]['REF'] = data['REF'].tolist()[0]
    dic_hmf_res[ix]['ALT'] = data['ALT'].tolist()[0]
    dic_hmf_res[ix]['TRANSCRIPT'] = data['Feature'].tolist()[0]
    dic_hmf_res[ix]['Protein_position'] = data['Protein_position'].tolist()[0]
    dic_hmf_res[ix]['SYMBOL'] = data['SYMBOL'].tolist()[0]
    dic_hmf_res[ix]['SAMPLES'] = len(set(data['SAMPLE'].tolist()))

df_res_HMF = pd.DataFrame(dic_hmf_res).T

#-------
# TCGA
#-------
tcga = pd.read_csv(tcga_vep, sep ='\t')

tcga['CHROM'] = tcga['ID'].apply(lambda x : x.split('_')[0])
tcga['POS'] = tcga['ID'].apply(lambda x : int(x.split('_')[1]))

dic_hmf_res = defaultdict(dict)

for ix, (idv, data) in enumerate(tcga.groupby(by='ID')):
    dic_hmf_res[ix]['MUTATION'] = '{}:{}:{}>{}'.format(data['CHROM'].tolist()[0].replace('chr', ''),
                                                   int(data['POS'].tolist()[0]),  
                                                  data['REF'].tolist()[0],
                                                  data['ALT'].tolist()[0],)
    dic_hmf_res[ix]['COHORT'] = 'OTHER_WXS_TCGA_FULL'
    dic_hmf_res[ix]['CONSEQUENCE'] = data['Consequence'].tolist()[0]
    dic_hmf_res[ix]['CHR'] = data['CHROM'].tolist()[0].replace('chr', '')
    dic_hmf_res[ix]['POS'] = int(data['POS'].tolist()[0])
    dic_hmf_res[ix]['REF'] = data['REF'].tolist()[0]
    dic_hmf_res[ix]['ALT'] = data['ALT'].tolist()[0]
    dic_hmf_res[ix]['TRANSCRIPT'] = data['Feature'].tolist()[0]
    dic_hmf_res[ix]['Protein_position'] = data['Protein_position'].tolist()[0]
    dic_hmf_res[ix]['SYMBOL'] = data['SYMBOL'].tolist()[0]
    dic_hmf_res[ix]['SAMPLES'] = len(set(data['SAMPLE'].tolist()))
    
df_res_tcga = pd.DataFrame(dic_hmf_res).T

# MERGE COHORT
df_cohort_merged = pd.concat([df_res_HMF, df_res_tcga])

#---------------------------------------
# RECUPERATE MUTATIONS THAT WERE FILTERED BECAUSE OF INTOGEN LIFTOVER to 38!
mutations_tcga = '{}/mutations.tsv'.format(path_intogen_run)
df_mut_allinto = pd.read_csv(mutations_tcga, sep ='\t')
wanted_genes = ['FOXP1', 'GNB1']
muts_missing = df_mut_allinto[(df_mut_allinto['SYMBOL'].isin(wanted_genes))&
                                (df_mut_allinto['COHORT']=='OTHER_WGS_HMF_FULL')]

merged_missing_cohorts = pd.concat([muts_missing, df_cohort_merged ])
merged_missing_cohorts.to_csv(outpath + 'mutations.tsv', sep ='\t', index = False, header = True)

#---------------------------------------
# 2 - Creating the gene table
#---------------------------------------

# load our list discovery
our_list_drivers = json.load(open(list_discovery))

# PATH TO INTOGEN RUN
driver_cohort_tcga = '{}/drivers.tsv'.format(path_intogen_run)
cohort_tcga = '{}/cohorts.tsv'.format(path_intogen_run)
mutations_tcga = '{}/mutations.tsv'.format(path_intogen_run)
unique_tcga = '{}/unique_drivers.tsv'.format(path_intogen_run)

df_tcga_unique = pd.read_csv(unique_tcga, sep ='\t')
df_cohort_tcga = pd.read_csv(cohort_tcga, sep ='\t')
df_mut_tcga = pd.read_csv(mutations_tcga, sep ='\t')
df_tcga = pd.read_csv(driver_cohort_tcga, sep ='\t')
df_cohort_tcga['CANCER_TYPE']= 'Clonal_Hematopoiesis'

#-----------------------
# HMF side of the table
#-----------------------

#hartwig = pd.read_csv(hmf_vep, sep ='\t')
hartwig = df2
hartwig_grhc38 = df_mut_tcga[df_mut_tcga['COHORT'] == 'OTHER_WGS_HMF_FULL']

total_samples_cohort_HMF = len(hartwig['SAMPLE'].unique())

dic_f_hmf = defaultdict(dict)
for gene, data in hartwig.groupby(by='SYMBOL'):
    total_samples = len(data['SAMPLE'].unique())
    mutations = len(data)
    percentage_samples = total_samples/total_samples_cohort_HMF
    dic_f_hmf[gene]['MUTATIONS'] = mutations
    dic_f_hmf[gene]['SAMPLES'] = total_samples
    dic_f_hmf[gene]['%_SAMPLES_COHORT'] = percentage_samples
 
cohorts_hmf = ['OTHER_WGS_HMF_FULL', 'OTHER_WGS_HMF_MOSAIC', 'OTHER_WGS_HMF_MUTECT']
cohorts_hmf_data = df_tcga[df_tcga['COHORT'].isin(cohorts_hmf)]
dic_final_res_hmf = defaultdict(dict)
for ix, (gene, data) in enumerate(cohorts_hmf_data.groupby(by='SYMBOL')):
    
    all_methods = [s.split(',') for s in data['METHODS'].tolist()]
    flattened = ','.join(list(set([val for sublist in all_methods for val in sublist])))
    role = data['ROLE'].value_counts().to_dict()
    role_cons = max(role, key=role.get)
    cgc = data['CGC_GENE'].tolist()[0]
    
    dic_final_res_hmf[ix]['SYMBOL'] = gene
    dic_final_res_hmf[ix]['TRANSCRIPT'] = data['TRANSCRIPT'].tolist()[0]
    dic_final_res_hmf[ix]['COHORT'] = 'OTHER_WGS_HMF_FULL'
    dic_final_res_hmf[ix]['CANCER_TYPE']= 'Clonal_Hematopoiesis'
    dic_final_res_hmf[ix]['CGC_GENE']= cgc
    cgc_t = data['CGC_CANCER_GENE'].tolist()[0]
    dic_final_res_hmf[ix]['CGC_CANCER_GENE']= cgc_t
    dic_final_res_hmf[ix]['METHODS']=flattened

    if gene in dic_f_hmf:
        dic_final_res_hmf[ix]['MUTATIONS']=dic_f_hmf[gene]['MUTATIONS'] 
        dic_final_res_hmf[ix]['SAMPLES']=dic_f_hmf[gene]['SAMPLES'] 
        dic_final_res_hmf[ix]['%_SAMPLES_COHORT']=dic_f_hmf[gene]['%_SAMPLES_COHORT']
    else:
        print(gene)
        dic_final_res_hmf[ix]['MUTATIONS'] = len(hartwig_grhc38[hartwig_grhc38['SYMBOL']==gene]['SAMPLES'])
        dic_final_res_hmf[ix]['SAMPLES']=hartwig_grhc38[hartwig_grhc38['SYMBOL']==gene]['SAMPLES'].sum()
        dic_final_res_hmf[ix]['%_SAMPLES_COHORT']=hartwig_grhc38[hartwig_grhc38['SYMBOL']==gene]['SAMPLES'].sum()/total_samples_cohort_HMF
    
    try:
        domain = data[~data['DOMAIN'].isnull()]['DOMAIN'].tolist()[0]
    except:
        domain = np.nan
    try:
        clust2 = data[~data['2D_CLUSTERS'].isnull()]['2D_CLUSTERS'].tolist()[0]
    except:
        clust2 = np.nan
    try:
        clust3 = data[~data['3D_CLUSTERS'].isnull()]['3D_CLUSTERS'].tolist()[0]
    except:
        clust3 = np.nan

    exc_mis = np.max(data[~data['EXCESS_MIS'].isnull()]['EXCESS_MIS'].tolist())
    exc_non = np.max(data[~data['EXCESS_NON'].isnull()]['EXCESS_NON'].tolist())
    exc_spl = np.max(data[~data['EXCESS_SPL'].isnull()]['EXCESS_SPL'].tolist())
    combi = np.min(data[~data['QVALUE_COMBINATION'].isnull()]['QVALUE_COMBINATION'].tolist())
    
    dic_final_res_hmf[ix]['QVALUE_COMBINATION']=combi
    dic_final_res_hmf[ix]['EXCESS_SPL']=exc_spl
    dic_final_res_hmf[ix]['EXCESS_NON']=exc_non
    dic_final_res_hmf[ix]['EXCESS_MIS']=exc_mis
    dic_final_res_hmf[ix]['3D_CLUSTERS']=clust3
    dic_final_res_hmf[ix]['2D_CLUSTERS']=clust2
    dic_final_res_hmf[ix]['DOMAIN']=domain
    
hmf_finlal = pd.DataFrame(dic_final_res_hmf).T

#-----------------------
# TCGA side of the table
#-----------------------
hartwig = pd.read_csv(tcga_vep, sep ='\t')
total_samples_cohort_TCGA = len(hartwig['SAMPLE'].unique())

dic_f = defaultdict(dict)
for gene, data in hartwig.groupby(by='SYMBOL'):
    total_samples = len(data['SAMPLE'].unique())
    mutations = len(data)
    percentage_samples = total_samples/total_samples_cohort_TCGA
    dic_f[gene]['MUTATIONS'] = mutations
    dic_f[gene]['SAMPLES'] = total_samples
    dic_f[gene]['%_SAMPLES_COHORT'] = percentage_samples
 
cohorts_hmf = ['OTHER_WXS_TCGA_FULL', 'OTHER_WXS_TCGA_MOSAIC']
cohorts_hmf_data = df_tcga[df_tcga['COHORT'].isin(cohorts_hmf)]
dic_final_res_hmf = defaultdict(dict)

for ix, (gene, data) in enumerate(cohorts_hmf_data.groupby(by='SYMBOL')):
    
    all_methods = [s.split(',') for s in data['METHODS'].tolist()]
    flattened = ','.join(list(set([val for sublist in all_methods for val in sublist])))
    role = data['ROLE'].value_counts().to_dict()
    role_cons = max(role, key=role.get)
    cgc = data['CGC_GENE'].tolist()[0]

    dic_final_res_hmf[ix]['SYMBOL'] = gene
    dic_final_res_hmf[ix]['TRANSCRIPT'] = data['TRANSCRIPT'].tolist()[0]
    dic_final_res_hmf[ix]['COHORT'] = 'OTHER_WXS_TCGA_FULL'
    dic_final_res_hmf[ix]['CANCER_TYPE']= 'Clonal_Hematopoiesis'
    dic_final_res_hmf[ix]['CGC_GENE']= cgc
    cgc_t = data['CGC_CANCER_GENE'].tolist()[0]
    dic_final_res_hmf[ix]['CGC_CANCER_GENE']= cgc_t

    dic_final_res_hmf[ix]['METHODS']=flattened
    dic_final_res_hmf[ix]['MUTATIONS']=dic_f[gene]['MUTATIONS'] 
    dic_final_res_hmf[ix]['SAMPLES']=dic_f[gene]['SAMPLES'] 
    dic_final_res_hmf[ix]['%_SAMPLES_COHORT']=dic_f[gene]['%_SAMPLES_COHORT']
    
    try:
        domain = data[~data['DOMAIN'].isnull()]['DOMAIN'].tolist()[0]
    except:
        domain = np.nan
    try:
        clust2 = data[~data['2D_CLUSTERS'].isnull()]['2D_CLUSTERS'].tolist()[0]
    except:
        clust2 = np.nan
    try:
        clust3 = data[~data['3D_CLUSTERS'].isnull()]['3D_CLUSTERS'].tolist()[0]
    except:
        clust3 = np.nan

    exc_mis = np.max(data[~data['EXCESS_MIS'].isnull()]['EXCESS_MIS'].tolist())
    exc_non = np.max(data[~data['EXCESS_NON'].isnull()]['EXCESS_NON'].tolist())
    exc_spl = np.max(data[~data['EXCESS_SPL'].isnull()]['EXCESS_SPL'].tolist())
    combi = np.min(data[~data['QVALUE_COMBINATION'].isnull()]['QVALUE_COMBINATION'].tolist())
    
    dic_final_res_hmf[ix]['QVALUE_COMBINATION']=combi
    dic_final_res_hmf[ix]['EXCESS_SPL']=exc_spl
    dic_final_res_hmf[ix]['EXCESS_NON']=exc_non
    dic_final_res_hmf[ix]['EXCESS_MIS']=exc_mis
    dic_final_res_hmf[ix]['3D_CLUSTERS']=clust3
    dic_final_res_hmf[ix]['2D_CLUSTERS']=clust2
    dic_final_res_hmf[ix]['DOMAIN']=domain

tcga_final = pd.DataFrame(dic_final_res_hmf).T
allowed_cohorts = ['OTHER_WGS_HMF_FULL', 'OTHER_WXS_TCGA_FULL']
merged_discovery_driv = pd.concat([tcga_final, hmf_finlal], sort=False)
merged_discovery_driv = merged_discovery_driv[merged_discovery_driv['SYMBOL'].isin(our_list_drivers)]
df_tcga_unique = df_tcga_unique[df_tcga_unique['SYMBOL'].isin(our_list_drivers)]

#---
# ADD IMPACT HERE
#---
mskcc = pd.read_csv(ch_cbio,
                    skiprows = 2, sep ='\t')
unique_samples = mskcc["Tumor_Sample_Barcode"].unique()


# this is the file from intogen processing from CH
mutations = pd.read_csv(ch_into, 
                       sep = "\t")

mutations["COHORT"] = "MSK_IMPACT"
dic_trans = dict(zip(mutations["SYMBOL"], mutations["TRANSCRIPT"]))

dic_res = pickle.load(open("../../data/tables/driver_disc_impact_toweb.pckl", "rb"))
dic_impact_method = defaultdict(list)
for g, d in dic_res.items():
    add_meth = set()
    for method, v in d.items():
        if v==1:
            add_meth.add(method.lower())
    dic_impact_method[g] = ','.join(add_meth)
    
res = pd.DataFrame(dic_res).T

res.columns = ["smregions", "dndscv", "oncodrivefml", "oncodriveclustl"]
cgc = pd.read_csv(cgc_parsed, 
                 sep ="\t")

cgct = cgc[~cgc["cancer_type"].isnull()]
cgct_aml = cgct[cgct["cancer_type"].str.contains("AML")]
set_aml = cgct_aml["Gene Symbol"].tolist()
cgc_list = cgc["Gene Symbol"].tolist()

remove_bad = ["5'Flank", "Intron", "5'UTR", "3'Flank", "3'UTR" ]
mskcc_m = mskcc[~mskcc["Variant_Classification"].isin(remove_bad)]

new_table = defaultdict(str)
tores = []

for i, row in res.iterrows():
    symbol = i
    trans = dic_trans[i]
    cohort = "MSK_IMPACT"
    typ = "Clonal_Hematopoiesis"
    cgc_value = "False"
    if i in cgc_list:
        cgc_value = "True"
    cgc_value2 = "False"
    if i in set_aml:
        cgc_value2 = "True"   
    
    meth = []
    for c in res.columns:
        if row[c]==1:
            meth.append(c)
    methods = ",".join(meth)
    
    number_mutations = len(mskcc_m[mskcc_m["Hugo_Symbol"]==i])
    sample_mut = len(mskcc_m[mskcc_m["Hugo_Symbol"]==i]["Tumor_Sample_Barcode"].unique())
    sample_per = sample_mut/len(unique_samples)
    
    tores.append((symbol, trans ,cohort, typ, cgc_value, cgc_value2, methods, number_mutations, 
                  sample_mut, sample_per ,
                1,1,1,1,np.nan,np.nan,np.nan))
        
res_table = pd.DataFrame(tores)
res_table.columns = merged_discovery_driv.columns

ensembl = pd.read_csv(biomart_cds, sep ="\t", 
                     names = ["ENSEMBL_GENE", "SYMBOL", "ENSEMBL_PROTEIN", 1,2,3,4,5,6,7, "ENSEMBL_TRANSCRIPT"])

msk_genes = res_table["SYMBOL"].tolist()
genes_to_add = [g for g in msk_genes if g not in df_tcga_unique["SYMBOL"].tolist()]
toadd = []
for g in genes_to_add:
    toadd.append(ensembl[ensembl["SYMBOL"]==g][["ENSEMBL_GENE", "SYMBOL", 
                                               "ENSEMBL_PROTEIN", "ENSEMBL_TRANSCRIPT"]])
    
all_new_d = pd.concat(toadd).drop_duplicates()
final_l_d = pd.concat([df_tcga_unique, all_new_d])
final_driver_table = pd.concat([merged_discovery_driv, res_table])
all_muts = pd.concat([df_mut_tcga[df_mut_tcga['COHORT'].isin(allowed_cohorts)], mutations])

coh = df_cohort_tcga[df_cohort_tcga['COHORT'].isin(allowed_cohorts)]
coh.loc[-1] = ["MSK_IMPACT", "Clonal_Hematopoiesis", "PANEL", "11011", 	"7215"]

# save result
coh.to_csv(outpath + 'cohorts.tsv', sep ='\t', index = False, header = True)
final_driver_table.to_csv(outpath + 'drivers.tsv', sep ='\t', index = False, header = True)
# TODO we should change this file as mutations input
all_muts.to_csv(outpath + 'mutations.tsv', sep ='\t', index = False, header = True)
final_l_d.to_csv(outpath + 'unique_drivers.tsv', sep ='\t', index = False, header = True)


########################
# Supplementary Table 2
########################

print("GENERATING THE TABLE...")
list_impact = "../../data/tables/list_genes_impact.json"

our_list_drivers = json.load(open(list_discovery))
res_list = json.load(open(category_CH))
list_impact =json.load(open(list_impact))

driver_cohort_meth = '../../data/tables/driver_tier1.tsv'
df_methods = pd.read_csv(driver_cohort_meth, sep ='\t')

list_genes_in_impact = []
with open('../../data/tables/list_genes_dndnsc.txt', 'rt') as infile:
    for line in infile:
        list_genes_in_impact.append(line.rstrip())

df_expression = pd.read_csv(expression_data, sep ='\t', 
                index_col = 0)

disc_tcga = df_methods[df_methods['COHORT'].str.contains('TCGA')]['SYMBOL'].tolist()
disc_hmf = df_methods[df_methods['COHORT'].str.contains('HMF')]['SYMBOL'].tolist()

towrite = []
keep_res = defaultdict(dict)

for l in set(our_list_drivers + list_impact):
    
    if l in dic_f:
        keep_res[l]['Number_mutations_TCGA'] = dic_f[l]['MUTATIONS']
        keep_res[l]['Number_samples_TCGA'] = dic_f[l]['SAMPLES']
        keep_res[l]['Proportion_samples_TCGA'] = dic_f[l]['%_SAMPLES_COHORT']
    else:
        keep_res[l]['Number_mutations_TCGA'] =0
        keep_res[l]['Number_samples_TCGA'] =0
        keep_res[l]['Proportion_samples_TCGA'] =0
    
    if l in dic_f_hmf:
        keep_res[l]['Number_mutations_HMF'] = dic_f_hmf[l]['MUTATIONS']
        keep_res[l]['Number_samples_HMF'] = dic_f_hmf[l]['SAMPLES']
        keep_res[l]['Proportion_samples_HMF'] = dic_f_hmf[l]['%_SAMPLES_COHORT']
    else:
        keep_res[l]['Number_mutations_HMF'] = hartwig_grhc38[hartwig_grhc38['SYMBOL']==l]['SAMPLES'].sum()
        keep_res[l]['Number_samples_HMF'] = hartwig_grhc38[hartwig_grhc38['SYMBOL']==l]['SAMPLES'].sum()
        keep_res[l]['Proportion_samples_HMF'] = hartwig_grhc38[hartwig_grhc38['SYMBOL']==l]['SAMPLES'].sum()/total_samples_cohort_HMF        
    
    for k in res_list:
        keep_res[l][k] = 'NO'
    for k in res_list:
        if l in res_list[k]:
            keep_res[l][k] = 'YES'
    keep_res[l]['Gene_in_IMPACT_panel'] = 'NO'
    if l in list_genes_in_impact:
        keep_res[l]['Gene_in_IMPACT_panel'] = 'YES'
            
    me = df_methods[df_methods['SYMBOL']==l]
    
    keep_res[l]['Discovery_TCGA'] = 'NO'
    keep_res[l]['Discovery_HMF'] = 'NO'
    keep_res[l]['Discovery_IMPACT'] = 'NO'

    if l in disc_tcga:
        keep_res[l]['Discovery_TCGA'] = 'YES'
        out = '{}\t{}\n'.format(l, 'TCGA_BLOOD')
        towrite.append(out)
    if l in disc_hmf:
        keep_res[l]['Discovery_HMF'] = 'YES'
        out = '{}\t{}\n'.format(l, 'HMF_BLOOD')
        towrite.append(out)
    
    if l in list_genes_in_impact:
        if l in list_impact:
            keep_res[l]['Discovery_IMPACT'] = 'YES'
            out = '{}\t{}\n'.format(l, 'IMPACT_BLOOD')
            towrite.append(out)
    
    mean_exp = np.mean(df_expression.loc[l])
    keep_res[l]['Mean_FPKM'] = mean_exp

# add discovery methods as a column
res_table = pd.DataFrame(keep_res).T
dic_methos = defaultdict(dict)
for cohort, data in final_driver_table.groupby(by="COHORT"):
    if "IMPACT" not in cohort:
        cohort = cohort.split("_")[2]
    dic_methos[cohort]= dict(zip(data["SYMBOL"], data["METHODS"]))

res_table["Gene"] = res_table.index

res_table["DISCOVERY_METHODS_HMF"] = res_table["Gene"].apply(lambda x : dic_methos["HMF"].get(x, ""))
res_table["DISCOVERY_METHODS_TCGA"] = res_table["Gene"].apply(lambda x : dic_methos["TCGA"].get(x, ""))
res_table["DISCOVERY_METHODS_IMPACT"] =  res_table["Gene"].apply(lambda x : dic_impact_method.get(x, ""))


###### 
# Add mean aging TCGA
######

df_t = pd.read_csv(clincal_data_tcga, sep ='\t', encoding = "ISO-8859-1", 
                  low_memory=False )

aged = df_t[df_t['age_at_initial_pathologic_diagnosis']!='[Not Available]']
aged = aged[aged['gender']!='[Not Available]']
aged['age_at_initial_pathologic_diagnosis'] = aged['age_at_initial_pathologic_diagnosis'].astype(int)
dic_aging = dict(zip(aged['bcr_patient_barcode'], aged['age_at_initial_pathologic_diagnosis']))

dic_mean_age = {}
for gene, data in tcga.groupby(by="SYMBOL"):
    dic_mean_age[gene] = np.nanmean([dic_aging.get(s, np.nan) for s in data["SAMPLE"].unique()])
res_table["Mean_Age"] = res_table["Gene"].map(dic_mean_age)

#
# Final touches
#
res_table.index.name = 'Gene'
res_table.drop("Gene", axis = 1, inplace = True)
res_table.sort_values(by="Proportion_samples_HMF", ascending = False, inplace = True)

res_table.to_csv('{}/supplementary_table_2.tsv'.format(outpath), 
                sep ='\t', index = True, header = True)