import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import tabix 
from pybedtools import BedTool
import os
from tqdm import tqdm
import sys
import matplotlib as mpl

pd.options.mode.chained_assignment = None

tqdm.pandas()

def config_params(font_size=7):

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'

def order_muts(type_mut):

    order = []
    if type_mut == 'snv':
        first = ['A', 'C', 'G', 'T']
        pyr = ['C', 'T']
        for p in pyr:
            for mut in first:
                if mut != p:
                    for f in first:
                        for f2 in first:
                            comb = '{}[{}>{}]{}'.format(f, p, mut, f2)
                            order.append(comb)

    return order

def process_variants(df, min_muts):
    order_m = order_muts('snv')
    dic_matrix = defaultdict(dict)
    for sample, data in df.groupby(by='SAMPLE'):
        if len(data)>min_muts:
            d = data['VARIANT_CLASS'].value_counts().to_dict()
            for order in order_m:
                dic_matrix[sample][order] = int(d.get(order, 0))

    ordered_df = pd.DataFrame(dic_matrix).loc[order_m]
    
    return ordered_df

def get_DBS(df):

    forbidden_ids = []
    keep_DBS_records = []
    df2 = df.copy()
    
    for chrom, data in df2[df2['TYPE']=='SNV'].groupby(by='CHROM'):
        
        data.sort_values(by='POS', inplace = True)
        index_done = data.index.tolist()
        data.reset_index(inplace = True)
        new_index = data.index.tolist()

        dic_ix = dict(zip(new_index, index_done))

        data['DIST'] = data['POS'].diff()

        closest_list = data[data['DIST']==1].index.tolist()
        forbidden_in_chrom = []

        for i in closest_list:
            row = data.loc[i]
            next_pos = i+1
            if next_pos in closest_list:

                forbidden_in_chrom.append(i)
                forbidden_in_chrom.append(next_pos)
                forbidden_ids.append(dic_ix[i])
                forbidden_ids.append(dic_ix[i-1])
                forbidden_ids.append(dic_ix[i+1])

            if i not in forbidden_in_chrom:
                previous_mut = data.loc[i-1]
                previous_mut['REF'] = '{}{}'.format(previous_mut['REF'], row['REF'])
                previous_mut['TYPE'] = 'DBS'
                previous_mut['ALT'] = '{}{}'.format(previous_mut['ALT'], row['ALT'])
                previous_mut['ID'] = 'chr{}_{}_{}_{}'.format(previous_mut['CHROM'], 
                                                            previous_mut['POS'], 
                                                            previous_mut['REF'], 
                                                            previous_mut['ALT'])
                previous_mut['VARIANT_CLASS'] = '{}_{}'.format(previous_mut['REF'], previous_mut['ALT'])

                keep_DBS_records.append(previous_mut)
                forbidden_ids.append(dic_ix[i])
                forbidden_ids.append(dic_ix[i-1])
                
    all_ix = df.index.tolist()
    good_ix = [ix for ix in all_ix if ix not in forbidden_ids]
    SNVs_df = df.loc[good_ix]
    DBS_df = pd.DataFrame(keep_DBS_records)
    if len(DBS_df):
        DBS_df.drop(['DIST', 'index'], axis = 1, inplace = True)
    df = pd.concat([SNVs_df, DBS_df], sort=False)
    
    return df

def annotate_in_PON(df, tb):

    chrom = df['CHROM']
    pos = df['POS']
    refv = df['REF']
    altv = df['ALT']

    records = tb.querys("{}:{}-{}".format(chrom, pos, pos))
    val_found = 0

    for l in records:
        chrom, pos, d1, ref, alt, d2, d3, info = l
        if ',' in alt:
            for ix, al in enumerate(alt.split(',')):
                if (ref == refv)&(al == altv):
                    info_d = {d.split('=')[0]:d.split('=')[1] for d in info.split(';') if '=' in d}
                    val_found = int(info_d['PON_COUNT'].split(',')[ix])

        else:
            if (ref == refv)&(alt == altv):
                info_d = {d.split('=')[0]:d.split('=')[1] for d in info.split(';') if '=' in d}
                val_found = int(info_d['PON_COUNT'])


    return val_found

def remove_common_variants(df, tb):
    
    chrom = df['CHROM']
    pos = df['POS']
    refv = df['REF']
    altv = df['ALT']

    records = tb.querys("{}:{}-{}".format(chrom, pos, pos))
    found = 'NO'
    
    for l in records:
        chrom, pos, d1, ref, alt, d2, d3, info = l
        if ',' in alt:
            for ix, al in enumerate(alt.split(',')):
                if (ref == refv)&(al == altv):
                    found = 'YES'

        else:
            if (ref == refv)&(alt == altv):
                found = 'YES'

    return found      

def process_segdups_grch38():
    '''
    Get all these files from UCSC and execute the commands so that you generate the all_repeats_to_remove files
    '''

    ## repeatmasker from UCSC
    # zcat repeat_masker_hg38.gz | tail -n +2 | cut -f6,7,8 | sed 's/chr//g' | sort -k1,1 -k2,2n | gzip > repeat_masker_hg38.sorted.gz
    df_repeats = '/home/opich/bg/mutpro/data/noncoding/repeat_masker_hg38.sorted.gz'
    
    
    #zcat segmental_duplications_hg38.gz | tail -n +2 |  cut -f2,3,4 | sed 's/chr//g' | sort -k1,1 -k2,2n | gzip >segmental_duplications_hg38.sorted.gz
    df_segmental_duplication = '/home/opich/bg/mutpro/data/noncoding/segmental_duplications_hg38.sorted.gz'
    
    #zcat simple_repeats_hg38.gz | tail -n +2 |  cut -f2,3,4 | sed 's/chr//g' | sort -k1,1 -k2,2n | gzip >simple_repeats_hg38.sorted.gz
    df_short_repeats = '/home/opich/bg/mutpro/data/noncoding/simple_repeats_hg38.sorted.gz'

    # merge all filters
    cmd = 'zcat /home/opich/bg/mutpro/data/noncoding/repeat_masker_hg38.sorted.gz /home/opich/bg/mutpro/data/noncoding/segmental_duplications_hg38.sorted.gz /home/opich/bg/mutpro/data/noncoding/simple_repeats_hg19.sorted.bed.gz | cut -f1,2,3 | sort -k1,1n -k2,2n | gzip > /home/opich/bg/mutpro/data/noncoding/all_repeats_to_remove.hg38.sorted.bed.gz'
    all_merged = '/home/opich/bg/mutpro/data/noncoding/all_repeats_to_remove.hg38.sorted.bed.gz'

    mut_bed = BedTool(all_merged)
    
    return mut_bed

def process_segdups_grch37():
    '''
    Get all these files from UCSC and execute the commands so that you generate the all_repeats_to_remove files
    '''
    
    ## repeatmasker from UCSC
    # zcat repeat_masker_hg38.gz | tail -n +2 | cut -f6,7,8 | sed 's/chr//g' | sort -k1,1 -k2,2n | gzip > repeat_masker_hg38.sorted.gz
    df_repeats = '/home/opich/bg/mutpro/data/noncoding/repeat_masker_hg19.gz'
    
    
    #zcat segmental_duplications_hg38.gz | tail -n +2 |  cut -f2,3,4 | sed 's/chr//g' | sort -k1,1 -k2,2n | gzip >segmental_duplications_hg38.sorted.gz
    df_segmental_duplication = '/home/opich/bg/mutpro/data/noncoding/segmental_duplications_hg19.sorted.gz'
    
    #zcat simple_repeats_hg38.gz | tail -n +2 |  cut -f2,3,4 | sed 's/chr//g' | sort -k1,1 -k2,2n | gzip >simple_repeats_hg38.sorted.gz
    df_short_repeats = '/home/opich/bg/mutpro/data/noncoding/simple_repeats_hg19.sorted.bed.gz'

    # merge all filters 
    cmd = 'zcat repeat_masker_sorted.bed.gz segmental_duplications_hg19.sorted.gz simple_repeats_hg19.sorted.bed.gz | cut -f1,2,3 | sort -k1,1n -k2,2n | gzip > all_repeats_to_remove.sorted.bed.gz'

    all_merged = '/home/opich/bg/mutpro/data/noncoding/all_repeats_to_remove.sorted.bed.gz'

    mut_bed = BedTool(all_merged)
    
    return mut_bed

def remove_duplicated_HMF(tcga, path_metadata):

    pat = pd.read_csv(path_metadata, sep='\t', encoding='latin-1')
    keep_samples = {}

    for hmf, data in pat.groupby(by='hmfPatientId'):
        l = data['sampleId'].tolist()
        for s in l:
            keep_samples[s] = hmf
            
    tcga['HEALTHY_SAMPLE'] = tcga['SAMPLE'].apply(lambda x : x.split('_')[1])
    tcga['TUMORAL_SAMPLE'] = tcga['SAMPLE'].apply(lambda x : x.split('_')[0])
        
    # keep only one sample per HMF patient
    dic_mutations_per_sample = tcga[tcga['count_repeat']==0]['TUMORAL_SAMPLE'].value_counts().to_dict()
    patient_done = set()
    sorted_samples = sorted(dic_mutations_per_sample, key=dic_mutations_per_sample.get, reverse=True)
    good_samples = set()

    for s in sorted_samples:
        if s in keep_samples:
            if keep_samples[s] not in patient_done:
                good_samples.add(s)
                patient_done.add(keep_samples[s])
        else:
            good_samples.add(s)

    tcga = tcga[tcga['TUMORAL_SAMPLE'].isin(good_samples)]
    # check additionaly if we have other samples with similar mutations (missanotation from HMF)
    dic_res = defaultdict(lambda : defaultdict(int))
    for i, data in tqdm(tcga.groupby(by='ID')):
        sample_list = data['TUMORAL_SAMPLE'].tolist()
        for s in sample_list:
            for s2 in sample_list:
                if s2!=s:
                    dic_res[s][s2]+=1
                    

    blacklisted = set()
    for k, v in dic_res.items():
        for s, count in v.items():
            perc = count/dic_mutations_per_sample[k.split('_')[0]]
            if perc > 0.5:
                # add to the blacklist the one with less mutations
                s1 = dic_mutations_per_sample[k.split('_')[0]]
                s2 = dic_mutations_per_sample[s.split('_')[0]]
                sel = min(s1, s2)
                if s1==sel:
                    blacklisted.add(k)
                else:
                    blacklisted.add(s)
                    
    tcga = tcga[~tcga['TUMORAL_SAMPLE'].isin(blacklisted)]

    return tcga

def filter_DBS(tcga, run, sage):

    if run =="TCGA":
        toconcat = []
        for sample, data in tqdm(tcga.groupby(by='SAMPLE')):
            out = get_DBS(data)
            toconcat.append(out)
        tcga = pd.concat(toconcat, sort = False, ignore_index = True)

    # if its HMF, we can check the DBS in the PoN too
    else:
        toconcat = []
        # select whether we have SNVs or others
        tcga['len_alt'] = tcga['ALT'].str.len()

        # number of characters in ref
        tcga['len_ref'] = tcga['REF'].str.len()

        # first classification between SNV and others
        tcga['TYPE'] = tcga.apply(
            lambda x: 'SNV' if ((x['len_alt'] == 1) and (x['len_ref'] == 1) and (x['ALT'] != '-') and (x['REF'] != '-')) else 'INDEL', axis=1
        )

        for sample, data in tqdm(tcga.groupby(by='SAMPLE')):
            out = get_DBS(data)
            toconcat.append(out)
            
        tcga = pd.concat(toconcat, sort = False, ignore_index = True)

        dbs_muts = tcga[tcga['TYPE']=='DBS']
        no_dbs_muts = tcga[tcga['TYPE']!='DBS']
        
        tb = tabix.open(sage)
        dbs_muts['PON'] = dbs_muts.progress_apply(annotate_in_PON, args =(tb, ), axis = 1)
        dbs_muts['NORMAL_PON'] = dbs_muts['PON']/1000
        dnmt3a_pon_hmf = 8/1000

        dbs_muts = dbs_muts[(dbs_muts['NORMAL_PON']<= dnmt3a_pon_hmf)]

        dbs_muts.drop(['NORMAL_PON', 'PON'], axis =1, inplace = True)

        tcga = pd.concat([no_dbs_muts, dbs_muts ], sort = False, ignore_index = True)
    return tcga

def get_males_TCGA():

    clincal_data_tcga = '/home/opich/bg/clonal_hematopoiesis_code/data/TCGA/clinical_PANCAN_patient_with_followup_151210.tsv'
    df_t = pd.read_csv(clincal_data_tcga, sep ='\t', encoding = "ISO-8859-1", 
                    low_memory=False )
    gender_patients = dict(zip(df_t['bcr_patient_barcode'], df_t['gender']))
    males = [s for s, v in gender_patients.items() if v == "MALE"]

    return males

def get_males_HMF(path_metadata):

    df_t = pd.read_csv(path_metadata, sep='\t', encoding='latin-1')
    gender_patients = dict(zip(df_t['sampleId'], df_t['gender']))
    males = [s for s, v in gender_patients.items() if v == "male"]

    return males

def main_processing(file, run, outpath, comprehensive = 0, metadata = "/home/opich/bg/mutpro/data/CH/HMF/metadata.tsv", sage = '/home/opich/bg/mutpro/data/noncoding/common_variants/SageGermlinePon.hg19.vcf.gz'):
    """
    This will generate a set of filtered mutations for all downstream analyses, both driver discovery and regressions.
    Briefly, it will:

    - remove low mappability defined by SDUST (same as in gnomAD publication)
    - filter mutations below a PoN treshold and gnomAD thresholds.
    - filter mutations by minimum number of reads
    - filter by VAF (unless comprehensive mode)
    - Extra step of removing common variants in dbSNP
    - remove hypermutant samples, likely to be enriched by artifacts
    - make sure that in HMF we dont have related samples with different IDs, which has happened in the past.

    Three levels of comprehensiveness:
    0: remove everything below 0.5 VAF
    1: remove everything below 0.5 except sexual chromosomes in males
    2: no VAF filters at all

    """

    os.makedirs(outpath, exist_ok=True)
    
    # get dbSNP info 
    if run =="TCGA":

        # TCGA is hg38
        file_common_ncbi_grch38 = tabix.open('/home/opich/bg/mutpro/data/noncoding/common_variants/00-common_all.grch38.vcf.gz')
        file_common_UCSC_grch38 = tabix.open('/home/opich/bg/mutpro/data/noncoding/common_variants/snp151Common.hg38.vep.txt.sorted.bgz')

    else:
        # HMF is hg19
        file_common_ncbi_grch38 = tabix.open('/home/opich/bg/mutpro/data/noncoding/common_variants/00-common_all.grch37.vcf.gz')
        file_common_UCSC_grch38 = tabix.open('/home/opich/bg/mutpro/data/noncoding/common_variants/snp151Common.hg19.vep.txt.sorted.bgz')

    jak2_af_exome_hg19 = 88/250262
    jak2_af_genome_hg19 = 9/31364
    dnmt3a_af_exome_hg19 = 55/251040
    dnmt3a_af_genome_hg19 = 9/31364

    jak2_pon = 12/4136
    dnmt3a_pon_hmf = 8/1000

    tcga_pon = dnmt3a_pon_hmf

    gnomad_filter_exome = np.max([jak2_af_exome_hg19, dnmt3a_af_exome_hg19 ])
    gnomad_filter_genome = np.max([jak2_af_genome_hg19, dnmt3a_af_genome_hg19 ])

    # Filter by SDUST and PON in TCGA, HMF was coded in a previous step
    if run =="TCGA":

        # read file from the calling merged
        tcga = pd.read_csv(file, sep ='\t')

        tcga_pon = jak2_pon
        ponsize = 4136   
        # create proportion of PoN. The size for TCGA is 4136 individuals
        tcga['NORMAL_PON'] = tcga['PON']/ponsize

        # remove mappability and proportion of tcga pon
        tcga = tcga[(tcga['SDUST']==0)&(tcga['NORMAL_PON']<= tcga_pon)]

    else:
        # read file from the calling merging
        list_cols = ["CHROM","POS-1","POS","VAR_COUNTS","VAF","REF","ALT","SAMPLE","VARIANT_CLASS","TRIPLET",
        "FOUND_GERMLINE","FOUND_SOMATIC","ID","PHASING","PHASING_PREDICTED","AC_genome","AF_genome","AC_exome","AF_exome","MUTECT","NEW_S", "minorCN"]

        # this will have a huge impact on RAM
        tcga = pd.read_csv(file, sep ='\t', low_memory = False, usecols = list_cols)
        tcga.rename(columns={"VAR_COUNTS": "var_reads"}, inplace = True)

    # remove unknown sites. Be aware this removes Y chromosome because it is missing from gnomAD.
    print("fixing gnomad annotations...")
    if run =="TCGA":
        if comprehensive == 0:
            tcga = tcga[tcga['AF_genome']!='unk']
        else:
            tcga["AF_genome"] = tcga["AF_genome"].apply(lambda x : float(str(x).replace("unk", "0")))
    
    else:
        no_found = tcga[(tcga['AF_genome']=='.')|(tcga['AF_exome']=='.')]
        list_replace = no_found.index.tolist()

        if len(list_replace):
            tcga.loc[list_replace, 'AF_exome'] = tcga.loc[list_replace]['AF_genome'].tolist()
            
        if comprehensive == 0:
            tcga = tcga[tcga['AF_exome']!='.']
            tcga = tcga[tcga['AF_genome']!='.']
            tcga = tcga[tcga['AF_exome']!='unk']
            tcga = tcga[tcga['AF_genome']!='unk']
        else:
            tcga["AF_exome"] = tcga["AF_exome"].apply(lambda x : float(str(x).replace("unk", "0")))
            tcga["AF_exome"] = tcga["AF_exome"].apply(lambda x : 0 if x == "." else x)
            tcga["AF_genome"] = tcga["AF_genome"].apply(lambda x : float(str(x).replace("unk", "0")))
            tcga["AF_genome"] = tcga["AF_genome"].apply(lambda x : 0 if x == "." else x)

    tcga['AF_exome'] = tcga['AF_exome'].astype(float)
    tcga['AF_genome'] = tcga['AF_genome'].astype(float)

    len_full_tcga = len(tcga)

    # select variants with more than 2 read supporting
    tcga = tcga[tcga['var_reads']>=2]
    len_reads_tcga = len(tcga)

    print("apply VAF corrections...")
    # select variants with less than 0.5 VAF
    tcga['VAF'] = tcga['VAF'].astype(float)

    # same as before, remove filter VAF unless we do a comprehensive run (for sexual chroms)
    if comprehensive == 0:
        tcga = tcga[tcga['VAF']<0.5]

    # First level of comprehensive.
    # Here remove the VAF filter in males in chromX and Y.
    elif comprehensive == 1:

        if run=="TCGA":
            males = get_males_TCGA()
            sam = "SAMPLE"
        else:
            males = get_males_HMF(metadata)
            sam = "NEW_S"
            
        female_samples = tcga[(~tcga[sam].isin(males))&(tcga["VAF"]<0.5)]
        male_samples =  tcga[(tcga[sam].isin(males))]
        nonsexual_males = male_samples[(~male_samples["CHROM"].isin(["X", "Y"]))&(male_samples["VAF"]<0.5)]
        sexual_males = male_samples[male_samples["CHROM"].isin(["X", "Y"])]

        tcga = pd.concat([female_samples, nonsexual_males, sexual_males], 
                        sort = False, ignore_index = True)

    len_vaf_tcga = len(tcga)

    print("len before gnomad", len(tcga))

    # remove variants in gnomAD
    tcga = tcga[(tcga['AF_exome']<=gnomad_filter_exome)&(tcga['AF_genome']<=gnomad_filter_genome)]

    print("len after gnomad", len(tcga))

    #--------------------------
    # get extra common variants
    #-------------------------
    print("labelling common SNPs...")
    tcga['COMMON_dbSNP'] = tcga.progress_apply(remove_common_variants, args =(file_common_ncbi_grch38, ), axis = 1)
    tcga['COMMON_UCSC'] = tcga.progress_apply(remove_common_variants, args =(file_common_UCSC_grch38, ), axis = 1)
    tcga = tcga[(tcga['COMMON_dbSNP']=='NO')&(tcga['COMMON_UCSC']=='NO')]

    print("len after COMMON_dbSNP", len(tcga))

    tcga.sort_values(by=['CHROM', 'POS'], inplace = True)

    #--------------------------
    # annotate repeated regions
    #-------------------------
    print("removing duplicated regions...")
    len_gnomad_tcga = len(tcga)

    if run =="TCGA":
        bed_rep_hg38 = process_segdups_grch38()
    else:
        bed_rep_hg38 = process_segdups_grch37()

    bed_TCGA = BedTool.from_dataframe(tcga)
    cols = list(tcga.columns) + ['count_repeat']
    out = bed_TCGA.intersect(bed_rep_hg38, c = True)
    tcga = out.to_dataframe(names = cols)
    len_tcga_mapp = len(tcga[tcga['count_repeat']==0])

    #------------------------------------
    # Check for duplicated samples in HMF
    #------------------------------------
    if run !="TCGA":
        tcga = remove_duplicated_HMF(tcga, metadata)

    #--------------
    # remove hypermutants (likely artifacts) after controlling for repeats etc
    #---------------
    dict_samples = tcga[tcga['count_repeat']==0]['SAMPLE'].value_counts().to_dict()
    percentile_to_cut =  np.percentile(list(dict_samples.values()), 97.5)
    samples_to_remove = [s for s, c in dict_samples.items() if c >= percentile_to_cut]
    tcga = tcga[~tcga['SAMPLE'].isin(samples_to_remove)]

    tcga['CHROMOSOME'] = tcga['CHROM']
    tcga['POSITION'] = tcga['POS']

    tcga.sort_values(by=['CHROM', 'POS'], inplace = True)

    tcga.to_csv("temp.gz", sep ="\t", index = False, header = True)

    tcga = filter_DBS(tcga, run, sage)
    len_hypermutants_tcga = len(tcga[tcga['count_repeat']==0])

    #---------------------------------
    # UNTIL HERE THE NORMAL PROCESSING, now saving the outputs
    #---------------------------------
    tcga_mosaic = tcga[(tcga['PHASING_PREDICTED'].astype(str).str.contains('mosaic'))&
        (~tcga['PHASING_PREDICTED'].astype(str).str.contains('low-confidence'))]

    # full
    path_out = '{}/{}_FULL'.format(outpath, run)

    print(path_out)
    os.makedirs(path_out, exist_ok=True)

    tcga[(tcga['count_repeat']==0)&(tcga["VAF"]<0.4)][['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE']].to_csv('{}/{}_full.muts.gz'.format(path_out, run), 
                                                                    sep ='\t', index = False, header = True, 
                                                                    compression = 'gzip')
                                                                    

    tcga.to_csv('{}/{}_full.tsv.gz'.format(path_out, run), 
            sep ='\t', index = False, header = True, compression = 'gzip')

    # -- MOSAIC
    path_out = '{}/{}_MOSAIC'.format(outpath, run)
    os.makedirs(path_out, exist_ok=True)
    tcga_mosaic = tcga[(tcga['PHASING_PREDICTED'].astype(str).str.contains('mosaic'))&
        (~tcga['PHASING_PREDICTED'].astype(str).str.contains('low-confidence'))]

    tcga_mosaic[(tcga_mosaic['count_repeat']==0)&(tcga["VAF"]<0.4)][['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE']].to_csv('{}/{}_mosaic.muts.gz'.format(path_out, run), 
                                                                    sep ='\t', index = False, header = True, 
                                                                    compression = 'gzip')
    tcga_mosaic.to_csv('{}/{}_mosaic.tsv.gz'.format(path_out, run), 
            sep ='\t', index = False, header = True, compression = 'gzip')

    #----------------
    # length mosaic
    len_hypermutants_tcga_mosaic = len(tcga_mosaic[tcga_mosaic['count_repeat']==0])
    #-----------------

    if run !="TCGA":
        # -- MUTECT
        tcga_mutect = tcga[tcga['MUTECT']=='YES']
        path_out = '{}/{}_MUTECT'.format(outpath, run)
        os.makedirs(path_out, exist_ok=True)
        tcga_mutect[(tcga_mutect['count_repeat']==0)&(tcga_mutect["VAF"]<0.4)][['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE']].to_csv('{}/{}_mutect.muts.gz'.format(path_out, run), 
                                                                    sep ='\t', index = False, header = True, 
                                                                    compression = 'gzip')

        #---------------
        # length mutect
        len_hypermutants_tcga_mutect = len(tcga_mutect[tcga_mutect['count_repeat']==0])
        #---------------

        tcga_mutect.to_csv('{}/{}_mutect.tsv.gz'.format(path_out, run),
                sep ='\t', index = False, header = True, compression = 'gzip')

        with open('{}/numbers_HMF.txt'.format(path_out), 'wt') as outfile:

            outfile.write(str(len_full_tcga)+'\n')
            outfile.write(str(len_reads_tcga)+'\n')
            outfile.write(str(len_vaf_tcga) + '\n')
            outfile.write(str(len_gnomad_tcga) + '\n') 
            outfile.write(str(len_tcga_mapp) + '\n') 
            outfile.write(str(len_hypermutants_tcga) + '\n')    
            outfile.write(str(len_hypermutants_tcga_mosaic) + '\n')   
            outfile.write(str(len_hypermutants_tcga_mutect) + '\n')   

file = sys.argv[1]
run = sys.argv[2]
outpath = sys.argv[3]
# the comprehensiveness was implemented to answer some of the reviewers questions
comprehensive = int(sys.argv[4])

main_processing(file, run, outpath, comprehensive)