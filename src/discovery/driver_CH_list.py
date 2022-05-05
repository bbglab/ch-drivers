import pandas as pd
import json

# previous list of drivers according to DECODE paper
CH_list = set(['DNMT3A', 'TET2','ASXL1', 'JAK2', 'TP53','PPM1D', 'IDH2', 
             'CBL', 'SF3B1', 'SRSF2', 'GNAS', 'KRAS', 'GNB1', 'NRAS', 
             'MYD88'])
             
# from Ellis paper
myeloid_genes = set(['ABL1','ALK','ARID1A','ARID2','ASXL1','ASXL2','ATRX','BAP1',
             'BCL10','BCL2','BCOR','BRAF','CALR','CBL','CDK4','CDKN1B',
             'CDKN2A','CDKN2B','CDKN2C','CEBPA','CHEK2','CREBBP','CRLF2','CSF1R',
             'CSF3R','CTCF','CTNNB1','DICER1','DNMT3A','DNMT3B','EED','EGFR',
             'EP300','ETV6','EZH2','FAM175A','FBXW7','FGFR2','FLT3','GATA1',
             'GATA2','GNAS','H3F3A','H3F3B','HRAS','IDH1','IDH2','IRF4',
             'JAK2','JAK3','KDM5C','KDM6A','KIT','KRAS','MGA','MPL','MYC',
             'NF1','NF2','NFE2L2','NOTCH1','NOTCH2','NPM1','NRAS','PAX5',
             'PIK3CA','PPM1D','PTEN','PTPN11','RAC1','RAD21','RAD50','RAD51',
             'RB1','RHOA','RRAS','RUNX1','SETD2','SF3B1','SH2B3','SRSF2','STAG2',
             'STAT3','STAT5A','SUZ12','TERT','TET2','TP53','U2AF1','WHSC1','WT1','ZRSR2',
              'NF1','DNMT3A','TET2','IKZF1','RAD21','WT1','KMT2D','SH2B3','TP53',
              'CEBPA','ASXL1','RUNX1','BCOR','KDM6A','STAG2','PHF6','KMT2C','PPM1D',
              'ATM','ARID1A','ARID2','ASXL2','CHEK2','CREBBP','ETV6',
              'EZH2','FBXW7','MGA','MPL','RB1','SETD2','SUZ12','ZRSR2'])

# Cancer gene census list

# this file is downloaded from the Cancer Gene Census website
# Download your own copy from CGC (Sanger website)
cgc_df = pd.read_csv('../../data/intogen/cancer_gene_census.csv')
set_cgc = set(cgc_df['Gene Symbol'].tolist())

res_list = {}
res_list['CH'] = list(CH_list)
res_list['Myeloid'] = list(myeloid_genes)
res_list['CGC'] =  list(set_cgc)

# dump categories
json.dump(res_list, open('../../data/tables/list_genes_CH_category.json', 
                        'wt'))