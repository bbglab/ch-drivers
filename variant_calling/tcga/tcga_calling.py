import os
import subprocess
from multiprocessing import Pool
import click
import gzip
import io
import tempfile
import json
import gzip
import pandas as pd
import tempfile

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


def get_index_bam(bam_file):
    cmd = r''' samtools index {} '''.format(bam_file)
    subprocess.check_output(cmd, shell=True)

def get_sample_bam(bam_file):
    cmd = r''' samtools view -H {} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq '''.format(bam_file)
    sample_name = subprocess.getoutput(cmd)
    return sample_name

def download_GCD(id1, id2, outpath, logpath, cpus):

    #for i in range(5):
    cmd = 'gdc-client download {0} -t {1} -d {2} -n {3} --log-file {4}/log.txt '.format(id1, path_TOKEN_TCGA, outpath, cpus, logpath)
    subprocess.check_output(cmd, shell=True)

    #for i in range(5):
    cmd = 'gdc-client download {0} -t {1} -d {2}  -n {3} --log-file {4}/log.txt '.format(id2, path_TOKEN_TCGA, outpath, cpus, logpath)
    subprocess.check_output(cmd, shell=True)


def do_strelka(blood_bam, tumoral_bam, outpath, cpus, sample_name):

    os.makedirs('{}/strelka/'.format(outpath), exist_ok = True)
    cmd = ''' python2 {} \
        --normalBam {} --tumorBam {} \
        --referenceFasta {} \
        --runDir {}/strelka/ --exome '''.format(path_strelka, tumoral_bam, blood_bam, path_genome_ref, outpath)

    # remove old folder if exists
    remove_old = '''rm -rf {}/strelka/ '''.format(outpath)
    subprocess.check_output(remove_old, shell=True)

    subprocess.check_output(cmd, shell=True)

    cmd = ''' python2 {}/strelka/runWorkflow.py -m local -j {}'''.format(outpath, cpus)
    subprocess.check_output(cmd, shell=True)

    mosaic_input = '{}/mosaic_input'.format(os.path.dirname(blood_bam))

    cmd = ''' zcat {}/strelka/results/variants/somatic.snvs.vcf.gz | grep -v "#" | grep -w PASS | cut -f1,2,4,5 > {} '''.format(outpath, mosaic_input)
    subprocess.check_output(cmd, shell=True)

    out_mosaic = '{}/mosaic'.format(outpath)
    os.makedirs(out_mosaic, exist_ok = True)

    muts = vcf_reader('{}/strelka/results/variants/somatic.snvs.vcf.gz'.format(outpath))
    pass_muts = muts[muts['FILTER']=='PASS']
    pass_muts['SAMPLE'] = sample_name

    pass_muts['POS-1'] = pass_muts['POS']-1
    pass_muts[['CHROM', 'POS-1', 'POS', 'REF', 'ALT', 'SAMPLE']].to_csv(mosaic_input, sep ='\t', header = False,
        index = False)

    cmd = 'python {0}/ReadLevel_Features_extraction.py {1} {2}/features.txt {3} /{4} {5} {6} bam'.format(path_mosaicfc, mosaic_input,
        out_mosaic, os.path.dirname(blood_bam), path_genome_ref, k24, cpus)

    subprocess.check_output(cmd, shell=True)

# path genome reference from TCGA. Download from https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files 
path_genome_ref = '/workspace/datasets/reverse_calling/TCGA/reference/GRCh38.d1.vd1.fa'


# these are dictionaries with ID matching files. I provide the dictionaries in the repo
dic_equivalent = json.load(open('/workspace/datasets/reverse_calling/TCGA/sample_2_bigID.json'))
dic_equivalent_file = json.load(open('/workspace/datasets/reverse_calling/TCGA/sample_to_filename.json'))

# this is the TOKEN for TCGA. https://gdc.cancer.gov/access-data/obtaining-access-controlled-data
path_TOKEN_TCGA = '/workspace/datasets/reverse_calling/TCGA/TOKEN.txt'

# this is the path to Strelka. https://github.com/Illumina/strelka
path_strelka = '/workspace/datasets/reverse_calling/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py'

# folder to store results
folder_out = '/workspace/datasets/reverse_calling/TCGA/calling/'

# path to mosaic forecast. https://github.com/parklab/MosaicForecast
path_mosaicfc = '/workspace/datasets/reverse_calling/stjudelife/software/MosaicForecast/'
#k24 for mosaic forecast. Refer to https://github.com/parklab/MosaicForecast to generate this file
k24 = '/workspace/projects/reverse_calling/src/sofware/MosaicForecast/resources/hg38/k24.umap.wg.bw'

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--blood_sample', required=True, help='name of sample', type=str)
@click.option('--tumoral_sample', required=True, help='name of sample', type=str)
@click.option('--cpus', required=True, help='number of parallel processes', type=int)
def run(blood_sample, tumoral_sample, cpus):

    outpath = '{}/{}/'.format(folder_out, dic_equivalent[blood_sample])
    final_sample = 'zcat {}/strelka/results/variants/somatic.snvs.vcf.gz'.format(outpath)

    if not os.path.isfile(final_sample):

        os.makedirs(outpath, exist_ok = True)
        dirpath = tempfile.mkdtemp()

        download_GCD(blood_sample, tumoral_sample, dirpath, outpath, cpus)

        bam_file_blood = '{}/{}/{}'.format(dirpath, blood_sample,  dic_equivalent_file[blood_sample])
        bam_file_tumor = '{}/{}/{}'.format(dirpath, tumoral_sample, dic_equivalent_file[tumoral_sample])

        get_index_bam(bam_file_blood)
        get_index_bam(bam_file_tumor)

        do_strelka(bam_file_blood, bam_file_tumor, outpath, cpus,  dic_equivalent_file[blood_sample].split('.bam')[0])

        cmd = '''rm {}'''.format(bam_file_blood)
        subprocess.check_output(cmd, shell=True)
        cmd = '''rm {}.bai'''.format(bam_file_blood)
        subprocess.check_output(cmd, shell=True)

        cmd = '''rm {}'''.format(bam_file_tumor)
        subprocess.check_output(cmd, shell=True)
        cmd = '''rm {}.bai'''.format(bam_file_tumor)
        subprocess.check_output(cmd, shell=True)

        # remove workspace to save space
        cmd = '''rm -rf {}/strelka/workspace/ '''.format(outpath)
        subprocess.check_output(cmd, shell=True)

if __name__ == '__main__':

    # run the script --> python tcga_calling.py --blood_sample XXX  --tumour_sample XXX --cpus XXX
    run()