import io
import os
import sys
import gzip
import tabix
import pandas as pd
from pybedtools import BedTool
import pybedtools
import numpy as np
from bgreference import hg19
import pandas as pd
import json

# vcf parser
def vcf_reader(path):
    with gzip.open(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(str.join(os.linesep, lines)),
        dtype={
            '#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str
        }, sep='\t', low_memory=False,
    ).rename(columns={'#CHROM': 'CHROM'})



def return_reads(row):
    
    list_format = row['FORMAT'].split(':')
    d = {list_format[i]: i for i in range(len(list_format))}
    colvaf_list = row["TUMOR"].split(':')
    
    reference_v = row['REF']
    alternate_v = row['ALT']
    
    # select only tier1
    reference = int(colvaf_list[d['{}U'.format(reference_v)]].split(',')[0])
    variant = int(colvaf_list[d['{}U'.format(alternate_v)]].split(',')[0])

    total_reads = reference + variant
    row['total_reads'] = total_reads
    row['ref_reads'] = reference
    row['var_reads'] = variant

    if total_reads > 0:
        row['VAF'] = row['var_reads'] / row['total_reads']
    else:
        row['VAF'] = np.nan
            
    return row

def return_indels_reads(row):
        
    list_format = row['FORMAT'].split(':')
    d = {list_format[i]: i for i in range(len(list_format))}
    colvaf_list = row["TUMOR"].split(':')
    
    total_depth = int(colvaf_list[d['DP']])
    #select only tier 1
    variant = int(colvaf_list[d['TIR']].split(',')[0])
    reference = total_depth-variant

    total_reads = reference + variant
    row['total_reads'] = total_reads
    row['ref_reads'] = reference
    row['var_reads'] = variant

    if total_reads > 0:
        row['VAF'] = row['var_reads'] / row['total_reads']
    else:
        row['VAF'] = np.nan
            
    return row


def get_mappable_regions(bed):
    
    path_mappable = '/workspace/projects/reverse_calling/data/mappability/hg19/hg19/k36.umap.bed.noheader.nochr.gz'
    mapp_bed = BedTool(path_mappable)
    mappable_mut = bed.intersect(mapp_bed, u=True)

    return mappable_mut

def SDUST_annot(bed):

    file_sdust = '/workspace/datasets/reverse_calling/mutect/refgenomes/Homo_sapiens.GRCh37.GATK.illumina.fasta.sdust.30.gz'
    mapp_bed = BedTool(file_sdust)
    mappable_mut = bed.intersect(mapp_bed, c=True)
    
    return mappable_mut

def get_DBS(df):

    forbidden_ids = []
    keep_DBS_records = []

    for chrom, data in df.groupby(by='CHROM'):
        
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
                previous_mut['ALT'] = '{}{}'.format(previous_mut['ALT'], row['ALT'])
                keep_DBS_records.append(previous_mut)
                forbidden_ids.append(dic_ix[i])
                forbidden_ids.append(dic_ix[i-1])
                
    all_ix = df.index.tolist()
    good_ix = [ix for ix in all_ix if ix not in forbidden_ids]
    SNVs_df = df.loc[good_ix]
    DBS_df = pd.DataFrame(keep_DBS_records)
    DBS_df.drop(['DIST', 'index'], axis = 1, inplace = True)
    df = pd.concat([SNVs_df, DBS_df], sort=False)
    
    return df

def get_mutect(df, f, s1, s2):
    
    mutect_file = '{}/output/mutect2/mutect2/{}_{}.filtered.gz'.format(f.split('/output/strelka/')[0], s2, s1)

    if os.path.isfile(mutect_file):
        df_mutect = vcf_reader(mutect_file)
        df_mutect = df_mutect[df_mutect['FILTER']=='PASS']
        #df_mutect['POS'] = df_mutect['POS'].astype(str)

        df_mutect['ID'] = df_mutect.apply(lambda x : '{}_{}_{}_{}'.format(x['CHROM'], x['POS'], x['REF'], 
                                                                       x['ALT']), axis = 1)
        pass_mutect = set(df_mutect['ID'].tolist())

        df['MUTECT'] = df['ID'].apply(lambda x: 'YES' if x in pass_mutect else 'NO')

    else:
        df['MUTECT'] = 'NOT_DONE'
    
    return df
        
# create order format similar to SigProfiler
def create_snv_class(df):
    pyr = ['C', 'T']
    rev = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N':'N'}
    x = df['TRIPLET']
    if x[1] in pyr:
        out = '{}[{}>{}]{}'.format(x[0], x[1], df['ALT'], x[2])
    else:
        out = '{}[{}>{}]{}'.format(rev[x[2]], rev[x[1]], rev[df['ALT']], rev[x[0]])
    return out


def add_vepannot(df, reader):
    forbidden_cons = ["downstream_gene_variant", "upstream_gene_variant", "intron_variant", "3_prime_UTR_variant", 
    "non_coding_transcript_variant","NMD_transcript_variant","intergenic_variant", "non_coding_transcript_exon_variant", 
                    "5_prime_UTR_variant", ]
    
    good_cons=[]
    for data in reader.get('chr'+str(df['CHROM']), df['POS'], df['POS']):
        forb = 0
        if (data[3] == df['ALT'])&(data[-2]=='YES'):
            for cons in data[7].split(','):
                if cons in forbidden_cons:
                    forb = 1
            if forb == 0:
                good_cons.append('.'.join(data))
    if len(good_cons):
        df['INFO'] = '|'.join(good_cons)

    else:
        df['INFO']=''
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


def annotate_gnomad(df, tb):

    try:
        records = tb.querys("{}:{}-{}".format(df['CHROM'], df['POS'], df['POS']))
        AC_AF = '0|0|NO' #1    15143   .   T   A   2142.27 RF  AC
        AC = 0
        AF = 0
        for l in records:
            chrom,pos, dot, ref, alt, qual, v,  info = l
            if ',' in alt:
                for ix,al in enumerate(alt.split(',')):
                    if pos == df["POS"]:
                        AC_AF = '0|0|YES'.format(AC, AF)

                    if (ref == df['REF'])&(al == df['ALT']):
                        info_d = {d.split('=')[0]:d.split('=')[1] for d in info.split(';')  if '=' in d}
                        AC = info_d.get('AC', 0).split(',')[ix]
                        AF = info_d.get('AF', 0).split(',')[ix]
                        AC_AF = '{}|{}|YES'.format(AC, AF)
                        break
            else:
                if pos == df["POS"]:
                    AC_AF = '0|0|YES'.format(AC, AF)

                if (ref == df['REF'])&(alt == df['ALT']):
                    info_d = {d.split('=')[0]:d.split('=')[1] for d in info.split(';')  if '=' in d}
                    AC = info_d.get('AC', 0)
                    AF = info_d.get('AF', 0)
                    AC_AF = '{}|{}|YES'.format(AC, AF)

                    break

    except:
        AC_AF = 'unk|unk|unk'

    return AC_AF


def class_DBS(row):
    
    comp = {'A': 'T', 'G': 'C'}
    complementary = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    dbs_list = {
        'AC_CA', 'AC_CG', 'AC_CT', 'AC_GA', 'AC_GG', 'AC_GT', 'AC_TA', 'AC_TG', 'AC_TT', 'AT_CA', 'AT_CC', 'AT_CG',
        'AT_GA', 'AT_GC', 'AT_TA', 'CC_AA', 'CC_AG', 'CC_AT', 'CC_GA', 'CC_GG', 'CC_GT', 'CC_TA', 'CC_TG', 'CC_TT',
        'CG_AT', 'CG_GC', 'CG_GT', 'CG_TA', 'CG_TC', 'CG_TT', 'CT_AA', 'CT_AC', 'CT_AG', 'CT_GA', 'CT_GC', 'CT_GG',
        'CT_TA', 'CT_TC', 'CT_TG', 'GC_AA', 'GC_AG', 'GC_AT', 'GC_CA', 'GC_CG', 'GC_TA', 'TA_AT', 'TA_CG', 'TA_CT',
        'TA_GC', 'TA_GG', 'TA_GT', 'TC_AA', 'TC_AG', 'TC_AT', 'TC_CA', 'TC_CG', 'TC_CT', 'TC_GA', 'TC_GG', 'TC_GT',
        'TG_AA', 'TG_AC', 'TG_AT', 'TG_CA', 'TG_CC', 'TG_CT', 'TG_GA', 'TG_GC', 'TG_GT', 'TT_AA', 'TT_AC', 'TT_AG',
        'TT_CA', 'TT_CC', 'TT_CG', 'TT_GA', 'TT_GC', 'TT_GG'
    }

    class_in = '{}_{}'.format(row['REF'], row['ALT'])
    if class_in not in dbs_list:

        # AC>CA	GT>TG
        class_in = '{}{}_{}{}'.format(
            complementary[row['REF'][1]], complementary[row['REF'][0]],
            complementary[row['ALT'][1]], complementary[row['ALT'][0]]
        )
    return class_in

def indel_mapp(f, mapped):

	indel_file = f.replace('snvs', 'indels')
	df_indel = vcf_reader(indel_file)
	df_indel = df_indel[df_indel['FILTER']=='PASS']
	df_indel['POS-5'] = df_indel['POS']-6
	df_indel['POS+5'] = df_indel['POS']+5
	bed_indel = BedTool.from_dataframe(df_indel[['CHROM', 'POS-5', 'POS+5']])
	bed_mapp = mapped.intersect(bed_indel, wao = True)

	return bed_mapp

def CNA_mapp(f, mapped):

	df_CNA = pd.read_csv(f, sep ='\t')
	bed_CNA = BedTool.from_dataframe(df_CNA[['chromosome', 'start', 'end', 'minorAllelePloidy', 'majorAllelePloidy']])
	bed_mapp = mapped.intersect(bed_CNA, wao = True)

	return bed_mapp


def flag_somatic_muts(df, tb):
    
    val_found = 'NO'
    try:
        records = tb.querys("{}:{}-{}".format(df['CHROM'], df['POS'], df['POS']))
        for l in records:
            if (l[6]=='PASS'):
                val_found = '{}>{}'.format(l[3], l[4])
    except:
        val_found = 'UNKNOWN'

    return val_found


def flag_germline_muts(df, tb):

    records = tb.querys("{}:{}-{}".format(df['CHROM'], df['POS'], df['POS']))
    val_found = 'NO'
    for l in records:
    	if (l[3] == df['REF'])&(l[4] == df['ALT'])&(l[6]=='PASS'):
            val_found = '{}>{}'.format(l[9], l[10])

    return val_found

# add phased mutations
def flag_phased(f, df):

    df['ID'] = df.apply(lambda x : '{}_{}_{}_{}'.format(x['CHROM'], x['POS'], x['REF'], x['ALT']), axis = 1)

    file_phasing = f.replace('output/strelka/results/variants/somatic.snvs.vcf.gz', 'output/mosaic/all.phasing.gz')
    df_phasing = pd.read_csv(file_phasing, sep ='\t')
    df_phasing['ID'] = df_phasing.apply(lambda x : '{}_{}_{}_{}'.format(x['chr'], x['pos'], x['ref'], x['alt']), axis = 1)
    dict_phasing = dict(zip(df_phasing['ID'], df_phasing['phasing']))
    df['PHASING'] = df['ID'].map(dict_phasing)

    file_phasing_predicted = f.replace('output/strelka/results/variants/somatic.snvs.vcf.gz', 'output/mosaic/predictions_test.tsv.gz')
    df_phasing_predicted = pd.read_csv(file_phasing_predicted, sep ='\t')
    df_phasing_predicted['REF'] = df_phasing_predicted['id'].apply(lambda x : x.split('~')[3])
    df_phasing_predicted['ALT'] = df_phasing_predicted['id'].apply(lambda x : x.split('~')[4])
    df_phasing_predicted['POS'] = df_phasing_predicted['id'].apply(lambda x : int(x.split('~')[2]))
    df_phasing_predicted['CHROM'] = df_phasing_predicted['id'].apply(lambda x : x.split('~')[1])
    df_phasing_predicted['ID'] = df_phasing_predicted.apply(lambda x : '{}_{}_{}_{}'.format(x['CHROM'], x['POS'], x['REF'], x['ALT']), axis = 1)
    dict_phasing_pred = dict(zip(df_phasing_predicted['ID'], df_phasing_predicted['prediction']))

    df['PHASING_PREDICTED'] = df['ID'].map(dict_phasing_pred)

    return df

# main function
def format_reverse(mutation_file, outpath, force = False):
    
    tb = tabix.open('/workspace/datasets/reverse_calling/SAGE/SageGermlinePon.hg19.vcf.gz')
    tb_gnomad = tabix.open('/workspace/datasets/gnomad/data/v2.1.1/genomes/hg19/gnomad.genomes.r2.1.1.sites.vcf.bgz')
    tb_gnomad_exomes = tabix.open('/workspace/datasets/gnomad/data/v2.1.1/exomes/hg19/gnomad.exomes.r2.1.1.sites.vcf.bgz')

    path_sample = mutation_file.split('/')[-6]
    dnmt3a_pon_hmf = 8/1000

    # some HMF samples have different naming conventions
    if 'dedup.realigned' in path_sample:
        sample_tum = path_sample.split('_')[0]
        sample_norm = path_sample.split('_')[2]

        sample = '{}_{}'.format(path_sample.split('_')[0], path_sample.split('_')[2])
    else:
        sample_tum = path_sample.split('.')[0]
        sample_norm = path_sample.split('_')[1].split('.')[0]

        sample = '{}_{}'.format(sample_tum, sample_norm)

    outfile_final = '{}/{}.muts.gz'.format(outpath, sample)

    if (not os.path.isfile(outfile_final)) or (force is True):
        
        # read vcf snvs and indels
        df_mut = vcf_reader(mutation_file)
        df_mut['CHROM'] = df_mut['CHROM'].astype(str)

        # get only canonical chromosomes
        wantedchroms = [str(i) for i in range(1, 23)]
        wantedchroms.append('Y')
        wantedchroms.append('X')

        # select only variants with the PASS filter and in the chromosomes we are interested in
        df_mut = df_mut[df_mut['CHROM'].isin(wantedchroms)]
        df_mut = df_mut[df_mut['FILTER'] == 'PASS']
        df_mut['FCLASS'] = 'SNV'
        df_mut = df_mut.apply(return_reads, axis=1)

        #-------
        #INDELS
        #-------
        indel_file = mutation_file.replace('snvs.vcf.gz', 'indels.vcf.gz')
        indel_vcf = vcf_reader(indel_file)
        indel_vcf['FCLASS'] = 'INDEL'

        indel_vcf = indel_vcf[indel_vcf['CHROM'].isin(wantedchroms)]
        indel_vcf = indel_vcf[indel_vcf['FILTER'] == 'PASS']
        indel_vcf = indel_vcf.apply(return_indels_reads, axis=1)

        # concat dataframes
        df_reads = pd.concat([df_mut, indel_vcf])

        df_reads['PON'] = df_reads.apply(annotate_in_PON, args =(tb, ), axis = 1)

        # this is the sample column
        lastcol = list(df_reads.columns)[-1]

        # select whether we have SNVs or others
        df_reads['len_alt'] = df_reads['ALT'].str.len()

        # number of characters in ref
        df_reads['len_ref'] = df_reads['REF'].str.len()

        # first classification between SNV and others
        df_reads['TYPE'] = df_reads.apply(
            lambda x: 'SNV' if ((x['len_alt'] == 1) and (x['len_ref'] == 1) and (x['ALT'] != '-') and (x['REF'] != '-')) else 'INDEL', axis=1
        )

        df_reads['pos-1'] = df_reads['POS'] - 1

        # get the triplet
        df_reads['TRIPLET'] = df_reads.apply(lambda x: hg19(x['CHROM'], x['pos-1'], 3), axis=1)
        df_reads['EXTENDED'] = df_reads.apply(lambda x: hg19(x['CHROM'], int(x['POS']) - 2, 5), axis=1)

        snv_df = df_reads[df_reads['TYPE'] != 'INDEL']

        snv_df['CLASS'] = 'SNV'
        snv_df['VARIANT_CLASS'] = snv_df.apply(create_snv_class, axis=1)

        indels_df = df_reads[df_reads['TYPE'] == 'INDEL']

        indels_df['CLASS'] = 'NOT_NEEDED'
        indels_df['VARIANT_CLASS'] = 'NOT_NEEDED'

        final_df = pd.concat([snv_df, indels_df])

        # assing the name of the sample
        final_df['sample'] = sample

        # create bed file
        mut_bed = BedTool.from_dataframe(final_df[[
            'CHROM', 'pos-1', 'POS', 'ref_reads', 'var_reads', 'VAF', 'total_reads', 'REF',
            'ALT', 'sample', 'TYPE', 'CLASS', 'VARIANT_CLASS', 'TRIPLET', 'EXTENDED', 'PON', 
        ]])

        # Remove unmappable regions
        mapped_regions = get_mappable_regions(mut_bed)

        mapped = SDUST_annot(mapped_regions)
        
        # label copy number status
        CNA_dict = json.load(open('/workspace/projects/reverse_calling/data/reversed_HMF/CN_HMF.json', 'rt'))
        path_CNAS = CNA_dict[sample_tum]

        mapped = CNA_mapp(path_CNAS, mapped)

        # merge to dataframe
        merge = mapped.to_dataframe(names=[
            'CHROM', 'POS-1', 'POS', 'REF_COUNTS', 'VAR_COUNTS', 'VAF', 'TOTAL_READS', 'REF', 'ALT', 'SAMPLE', 'TYPE',
            'CLASS', 'VARIANT_CLASS', 'TRIPLET', 'EXTENDED','PON', 'SDUST', 
             'CN1', 'CN2', 'CN3','minorCN', 'majorCN','MAP_CN',
        ])

        merge['NORMAL_PON'] = merge['PON']/1000
        merge = merge[(merge['SDUST']==0)&(merge['NORMAL_PON']<= dnmt3a_pon_hmf)]


        #-----
        # The following step is not needed unless you want to check the variants in hartwig. The Json will point towards the HMF files
        #-----

        # label somatics
        # this is an adhoc dictionary where the first key is the name of the sample, and the value is the location of the somatic mutations vcf, for example:
        #{"SAMPLE1": "/workspace/datasets/hartwig/SAMPLE1.purple.somatic.vcf.gz",}

        somatic_dict = json.load(open('/workspace/projects/reverse_calling/data/reversed_HMF/mut_HMF.json', 'rt'))
        path_somatics = somatic_dict[sample_tum]
        tb = tabix.open(path_somatics)
        merge['FOUND_SOMATIC'] = merge.apply(flag_somatic_muts, args =(tb, ), axis = 1)
        
        # label germline
        # dictionary same as before but pointing towards the germline file
        germline_dict = json.load(open('/workspace/projects/reverse_calling/data/reversed_HMF/germline_HMF.json', 'rt'))


        if sample_norm[:-1] in germline_dict:
            path_germline = germline_dict[sample_norm[:-1]]
            tb = tabix.open(path_germline)
            merge['FOUND_GERMLINE'] = merge.apply(flag_germline_muts, args =(tb, ), axis = 1)
        else:
            merge['FOUND_GERMLINE'] = np.nan
            
        merge = flag_phased(mutation_file, merge,)

        merge['gnomAD_genome'] = merge.apply(annotate_gnomad, args =(tb_gnomad, ), axis = 1)

        merge['AC_genome'] = merge['gnomAD_genome'].apply(lambda x : x.split('|')[0])
        merge['AF_genome'] = merge['gnomAD_genome'].apply(lambda x : x.split('|')[1])
        merge['Poly_Genome'] = merge['gnomAD_genome'].apply(lambda x : x.split('|')[2])

        merge['gnomAD_exomes'] = merge.apply(annotate_gnomad, args =(tb_gnomad_exomes, ), axis = 1)

        merge['AC_exome'] = merge['gnomAD_exomes'].apply(lambda x : x.split('|')[0])
        merge['AF_exome'] = merge['gnomAD_exomes'].apply(lambda x : x.split('|')[1])
        merge['Poly_Exome'] = merge['gnomAD_genome'].apply(lambda x : x.split('|')[2])

        # clean BedTools temp files
        pybedtools.cleanup()

        merge = get_mutect(merge, mutation_file, sample_tum, sample_norm)

        merge['NEW_S'] = merge['SAMPLE'].apply(lambda x : x.split('_')[0])

        merge[['CHROM', 'POS-1', 'POS', 'REF_COUNTS', 'VAR_COUNTS', 'VAF', 'TOTAL_READS', 'REF', 'ALT', 'SAMPLE', 'TYPE',
            'CLASS', 'VARIANT_CLASS', 'TRIPLET', 'EXTENDED','PON', 'SDUST',
             'CN1', 'CN2', 'CN3', 'minorCN', 'majorCN','MAP_CN', 'FOUND_GERMLINE', 'FOUND_SOMATIC', 'ID', 'PHASING', 'PHASING_PREDICTED', 
            'AC_genome', 'AF_genome','AC_exome', 'AF_exome', 'MUTECT', 'NEW_S', 'Poly_Genome', 'Poly_Exome']].to_csv(outfile_final, sep ='\t', header = True, index = False, 
        	compression = 'gzip')



# example of file
# /workspace/datasets/hartwig/20200504/oriol_pipeline/DRUP01100007T_dedup.realigned.cram_DRUP01100007R_dedup.realigned.cram/output/strelka/results/variants/somatic.snvs.vcf.gz

file = sys.argv[1]
outpath = sys.argv[2]

format_reverse(file, outpath)

#--
# FILTERS
#--
# A) Remove those found more than once in PON HMF
# B) FLAG mutations near germline indels (within 5 bp), because these are mostly false positives due to mismapping.
# C) FLAG copy number status in the mutations
# D) FLAG if mutation is observed in the somatic calling
# E) FLAG if the mutations is observed in the germline
# F) FLAG if phase
# G) FLAG if mosaic predicted