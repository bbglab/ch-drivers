import tabix
import gzip
import pandas as pd
import io
import sys
import os
import numpy as np
from pybedtools import BedTool
from bgreference import hg38
import pybedtools
from pyliftover import LiftOver


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

def create_snv_class(df):
    pyr = ['C', 'T']
    rev = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N':'N'}
    x = df['TRIPLET']
    if x[1] in pyr:
        out = '{}[{}>{}]{}'.format(x[0], x[1], df['ALT'], x[2])
    else:
        out = '{}[{}>{}]{}'.format(rev[x[2]], rev[x[1]], rev[df['ALT']], rev[x[0]])
    return out


# annotate PoN
def annotate_PoN(chrom, pos, refv, altv, tb):

    val_found = 0

    records = tb.querys("{}:{}-{}".format(chrom, pos, pos))
    for l in records:
        chrom, pos, d1, ref, alt, d2, d3, info = l
        if ',' in alt:
            for ix, al in enumerate(alt.split(',')):
                if (ref == refv)&(al == altv):
                    info_d = {d.split('=')[0]:d.split('=')[1] for d in info.split(';') if '=' in d}
                    val_found = int(info_d['AC'].split(',')[ix])


        else:
            if (ref == refv)&(alt == altv):
                info_d = {d.split('=')[0]:d.split('=')[1] for d in info.split(';') if '=' in d}
                val_found = int(info_d['AC'])

    return val_found


def get_PON(df, tb):

    chrom = df['CHROM']
    pos = df['POS']
    refv = df['REF']
    altv = df['ALT']

    val_found = annotate_PoN(chrom, pos, refv, altv, tb)

    return val_found

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

    reference_v = row['REF']
    alternate_v = row['ALT']

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


def do_liftover_test(df, lo):
    try:
        lov = lo.convert_coordinate("chr"+str(df['CHROM']), int(df['POS']))[0]
        df['CHROMOSOME_HG19'] = lov[0].replace('chr', '')
        df['POSITION_HG19'] = int(lov[1])
    except:
        df['CHROMOSOME_HG19'] = None
        df['POSITION_HG19'] = None

    return df


def get_mappable_regions(bed):

    path_mappable = '/workspace/projects/reverse_calling/data/mappability/hg38/hg38/k36.umap.bed.noheader.gz'
    mapp_bed = BedTool(path_mappable)
    mappable_mut = bed.intersect(mapp_bed, u=True)

    return mappable_mut

def SDUST_annot(bed):

    file = '/workspace/datasets/genomes/iGenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.sdust.30.gz'
    mapp_bed = BedTool(file)
    mappable_mut = bed.intersect(mapp_bed, c=True)

    return mappable_mut

def annotate_gnomad(df, tb):

    try:
        records = tb.querys("chr{}:{}-{}".format(df['CHROM'], df['POS'], df['POS']))
        AC_AF = '0|0' #1    15143   .   T   A   2142.27 RF  AC
        for l in records:
            chrom,pos, dot, ref, alt, qual, v,  info = l
            if ',' in alt:
                for ix,al in enumerate(alt.split(',')):
                    if (ref == df['REF'])&(al == df['ALT']):
                        info_d = {d.split('=')[0]:d.split('=')[1] for d in info.split(';') if '=' in d}
                        AC = info_d.get('AC', 0).split(',')[ix]
                        AF = info_d.get('AF', 0).split(',')[ix]
                        AC_AF = '{}|{}'.format(AC, AF)
                        break
            else:
                if (ref == df['REF'])&(alt == df['ALT']):
                    info_d = {d.split('=')[0]:d.split('=')[1] for d in info.split(';') if '=' in d}
                    AC = info_d.get('AC', 0)
                    AF = info_d.get('AF', 0)
                    AC_AF = '{}|{}'.format(AC, AF)

                    break

    except:
        AC_AF = 'unk|unk'

    return AC_AF

def annotate_gnomad_hg19(df, tb):

    try:
        records = tb.querys("{}:{}-{}".format(df['CHROMOSOME_HG19'], df['POSITION_HG19'], df['POSITION_HG19']))
        AC_AF = '0|0' #1    15143   .   T   A   2142.27 RF  AC
        for l in records:
            chrom,pos, dot, ref, alt, qual, v,  info = l
            if ',' in alt:
                for ix,al in enumerate(alt.split(',')):
                    if (ref == df['REF'])&(al == df['ALT']):
                        info_d = {d.split('=')[0]:d.split('=')[1] for d in info.split(';') if '=' in d}
                        AC = info_d.get('AC', 0).split(',')[ix]
                        AF = info_d.get('AF', 0).split(',')[ix]
                        AC_AF = '{}|{}'.format(AC, AF)
                        break
            else:
                if (ref == df['REF'])&(alt == df['ALT']):
                    info_d = {d.split('=')[0]:d.split('=')[1] for d in info.split(';') if '=' in d}
                    AC = info_d.get('AC', 0)
                    AF = info_d.get('AF', 0)
                    AC_AF = '{}|{}'.format(AC, AF)

                    break

    except:
        AC_AF = 'unk|unk'

    return AC_AF

# add phased mutations
def flag_phased(f, df):

    df['ID'] = df.apply(lambda x : '{}_{}_{}_{}'.format(x['CHROM'], x['POS'], x['REF'], x['ALT']), axis = 1)
    file_phasing_predicted = f.replace('strelka/results/variants/somatic.snvs.vcf.gz', 'mosaic/predictions.txt')

    if os.path.isfile(file_phasing_predicted):

        df_phasing_predicted = pd.read_csv(file_phasing_predicted, sep ='\t')
        df_phasing_predicted['REF'] = df_phasing_predicted['id'].apply(lambda x : x.split('~')[3])
        df_phasing_predicted['ALT'] = df_phasing_predicted['id'].apply(lambda x : x.split('~')[4])
        df_phasing_predicted['POS'] = df_phasing_predicted['id'].apply(lambda x : int(x.split('~')[2]))
        df_phasing_predicted['CHROM'] = df_phasing_predicted['id'].apply(lambda x : x.split('~')[1])
        df_phasing_predicted['ID'] = df_phasing_predicted.apply(lambda x : '{}_{}_{}_{}'.format(x['CHROM'], x['POS'], x['REF'], x['ALT']), axis = 1)
        dict_phasing_pred = dict(zip(df_phasing_predicted['ID'], df_phasing_predicted['prediction']))

        df['PHASING_PREDICTED'] = df['ID'].map(dict_phasing_pred)
    else:
        df['PHASING_PREDICTED'] = 'NOT_AVAILABLE'

    return df


def format_reverse(file, outpath, force = False):

    name_sample = os.path.dirname(file).split('/')[-4]
    outfile_final = '{}/{}.muts.gz'.format(outpath, name_sample)

    if (not os.path.isfile(outfile_final)) or (force is True):

        # this requires CHR
        tb_gnomad_full = tabix.open('/workspace/datasets/gnomad_2020/data/v2.1.1/exomes/hg38/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz')

        # this is in hg19 so we will do a liftover at the end. This doesn't require chr
        tb_gnomad_genomes = tabix.open('/workspace/datasets/gnomad_2020/data/v2.1.1/genomes/hg19/gnomad.genomes.r2.1.1.sites.vcf.bgz')

        # pon from TCGA. You can download this from TCGA website
        tb_PON = tabix.open('/workspace/datasets/reverse_calling/TCGA/PON/6b45b9f7-893e-4947-83b6-db0402471e23/MuTect2.PON.4136.vcf.gz')

        # get only canonical chromosomes
        wantedchroms = ['chr{}'.format(i) for i in range(1, 23)]
        wantedchroms.append('chrY')
        wantedchroms.append('chrX')

        df = vcf_reader(file)
        pass_df = df[df['FILTER'] =='PASS']
        pass_df = pass_df[pass_df['CHROM'].isin(wantedchroms)]


        pass_df['SAMPLE'] = name_sample
        pass_df['PON'] = pass_df.apply(get_PON, args = (tb_PON, ), axis = 1)
        pass_df = pass_df.apply(return_reads, axis = 1)

        # get indel
        indel_file = file.replace('.snvs.vcf.gz', '.indels.vcf.gz')
        df_ind = vcf_reader(indel_file)
        df_ind = df_ind[df_ind['FILTER'] =='PASS']
        df_ind = df_ind[df_ind['CHROM'].isin(wantedchroms)]

        df_ind['SAMPLE'] = name_sample
        df_ind['PON'] = df_ind.apply(get_PON, args = (tb_PON, ), axis = 1)
        df_ind = df_ind.apply(return_indels_reads, axis = 1)

        out = pd.concat([pass_df, df_ind])
        out['POS-1'] = out['POS']-1
        out_bed = BedTool.from_dataframe(out[
            ['CHROM', 'POS-1', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
            'NORMAL', 'TUMOR', 'SAMPLE', 'PON', 'total_reads', 'ref_reads',
            'var_reads', 'VAF']
        ])

        # Remov
        mapped_regions = get_mappable_regions(out_bed)
        mapped = SDUST_annot(mapped_regions)

        # merge to dataframe
        out = mapped.to_dataframe(names=[
            'CHROM', 'POS-1','POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
            'NORMAL', 'TUMOR', 'SAMPLE', 'PON', 'total_reads', 'ref_reads',
            'var_reads', 'VAF', 'SDUST'])

        out = flag_phased(file, out)

        out['CHROM'] = out['CHROM'].apply(lambda x : x.replace('chr', ''))
        out['gnomAD_exomes'] = out.apply(annotate_gnomad, args =(tb_gnomad_full, ), axis = 1)
        out['AC_exome'] = out['gnomAD_exomes'].apply(lambda x : x.split('|')[0])
        out['AF_exome'] = out['gnomAD_exomes'].apply(lambda x : x.split('|')[1])

        # get the triplet
        out['TRIPLET'] = out.apply(lambda x: hg38(x['CHROM'], x['POS-1'], 3), axis=1)
        out['EXTENDED'] = out.apply(lambda x: hg38(x['CHROM'], int(x['POS']) - 2, 5), axis=1)
        # select whether we have SNVs or others
        out['len_alt'] = out['ALT'].str.len()

        # number of characters in ref
        out['len_ref'] = out['REF'].str.len()

        # first classification between SNV and others
        out['TYPE'] = out.apply(
            lambda x: 'SNV' if ((x['len_alt'] == 1) and (x['len_ref'] == 1) and (x['ALT'] != '-') 
            and (x['REF'] != '-')) else 'INDEL', axis=1
        )

        snv_df = out[out['TYPE'] != 'INDEL']
        snv_df['CLASS'] = 'SNV'
        snv_df['VARIANT_CLASS'] = snv_df.apply(create_snv_class, axis=1)

        indels_df = out[out['TYPE'] == 'INDEL']
        indels_df['CLASS'] = 'NOT_NEEDED'
        indels_df['VARIANT_CLASS'] = 'NOT_NEEDED'

        final_df = pd.concat([snv_df, indels_df])

        lo = LiftOver('hg38', 'hg19')
        final_df = final_df.apply(do_liftover_test, args = (lo, ), axis = 1)

        final_df['POSITION_HG19'] =final_df['POSITION_HG19'].fillna(-1)
        final_df['POSITION_HG19'] = final_df['POSITION_HG19'].astype(int)
        final_df['POSITION_HG19'] = final_df['POSITION_HG19'].replace(-1, np.nan)

        final_df['gnomAD_genome'] = final_df.apply(annotate_gnomad_hg19, args =(tb_gnomad_genomes, ), axis = 1)

        final_df['AC_genome'] = final_df['gnomAD_genome'].apply(lambda x : x.split('|')[0])
        final_df['AF_genome'] = final_df['gnomAD_genome'].apply(lambda x : x.split('|')[1])

        # clean BedTools temp files
        pybedtools.cleanup()

        final_df.to_csv(outfile_final, sep ='\t',
                        index = False, header = True, compression = 'gzip')

file = sys.argv[1]
outpath = sys.argv[2]

format_reverse(file, outpath)