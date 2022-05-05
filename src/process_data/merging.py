from glob import glob 
import pandas as pd
from tqdm import tqdm 

# path to the output files
path_out = '/workspace/projects/reverse_calling/data/results_calling/'

# path HMF
path_res = '/workspace/projects/reverse_calling/data/reversed_HMF_gnomad_20200918/'
path_files = glob('{}/*.muts.gz'.format(path_res))
toconcat = []
for f in tqdm(path_files):
    df = pd.read_csv(f, sep ='\t')
    toconcat.append(df)

concat = pd.concat(toconcat)
concat.to_csv('{}/HMF_rev_muts_gnomad.hg19.tsv.gz'.format(path_out), sep ='\t', index = False,
             header = True, compression = 'gzip')

# path TCGA
path_res = '/workspace/datasets/reverse_calling/TCGA/filtered_20200918/'
path_files = glob('{}/*.muts.gz'.format(path_res))
toconcat=[]
for f in tqdm(path_files):
    df = pd.read_csv(f, sep ='\t')
    toconcat.append(df)
concat = pd.concat(toconcat)

concat.to_csv('{}/TCGA_rev_muts_gnomad.hg38.tsv.gz'.format(path_out), sep ='\t', index = False,
             header = True, compression = 'gzip')