import pandas as pd 
import numpy as np 

outpath = "../../data/tables/"

# data obtained from Cbioportal 
path_to_mskcc = '~/projects/collaborations/bbglab/CH/mutpro/data/other_CH_datasets/data_mutations_extended.txt'


mskcc = pd.read_csv(path_to_mskcc,
                    skiprows = 2, sep ='\t')
                
list_impact = mskcc["Hugo_Symbol"].unique()
with open("{}/impact_list_genes.tsv".format(outpath), "wt") as outfile:
    for l in list_impact:
        outfile.write(l + "\n")

gene_set = set(mskcc['Hugo_Symbol'].tolist())
with open('{}/list_genes_dndnsc.txt'.format(outpath), 'wt') as outf:
    for g in gene_set:
        out = '{}\n'.format(g)
        outf.write(out)

set_muts = mskcc[['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode']].drop_duplicates()
set_muts.columns = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE']

set_muts.to_csv('{}/ch_muts.muts.gz'.format(outpath), 
               sep ='\t', index = False, header = True, compression = 'gzip')

set_muts.columns = ["chr", "pos",  "ref",  "mut", "sampleID"]
set_muts = set_muts[set_muts['chr'].isin([str(s) for s in list(np.arange(23))] + ['X', 'Y'])]

set_muts.dropna()[['sampleID', 'chr', 'pos', 'ref', 'mut']].to_csv('{}/ch_muts.muts.gz'.format(outpath), 
               sep ='\t', index = False, header = True, compression = 'gzip')