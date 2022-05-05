import csv
import json
from tqdm import tqdm
import os
from collections import defaultdict
import csv
import gzip
import contextlib
import gzip
import pandas as pd 
import sys

def _open(file):
    if file.endswith('.gz'):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')
    
class DatasetError(Exception):
    pass

class TSVReader:
    def __init__(self, file):
        self._filename = file
        self._f = None

    def __enter__(self):
        self._f = _open(self._filename)
        reader = csv.DictReader(self._f, delimiter='\t')
        return reader

    def __exit__(self, exc_type, exc_value, traceback):
        self._f.close()

    def __iter__(self):
        with self as reader:
            for row in reader:
                yield row


@contextlib.contextmanager
def out_open(file=None, mode='w'):
    """
    Open file in write mode or send to STDOUT
    if file is None.
    Args:
        file: file path
        mode:
    Returns:
        File descriptor.
    """
    if file is None:
        fh = sys.stdout
    elif file.endswith('.gz'):
        fh = gzip.open(file, mode)
    else:
        fh = open(file, mode)
    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()

                
VALID_CONSEQUENCES = {
        "transcript_ablation", "splice_donor_variant", "splice_acceptor_variant", "stop_gained", "frameshift_variant",
        "stop_lost", "initiator_codon_variant", "transcript_amplification", "inframe_insertion", "inframe_deletion",
        "missense_variant", "splice_region_variant", "incomplete_terminal_codon_variant", "stop_retained_variant",
        "synonymous_variant", "coding_sequence_variant"
    }
TRANSCRIPTS = set()
GENES = set()
transcripts_file =  '../data/misc/ensembl_canonical_transcripts.tsv'

with open(transcripts_file) as fd:
    for r in csv.reader(fd, delimiter='\t'):
        TRANSCRIPTS.add(r[1])
        GENES.add(r[0])
        
def valid_consequence(consequences):
    csq_set = {c for c in consequences.split(',')}
    return len(csq_set & VALID_CONSEQUENCES) > 0

def filter_(file, stats):
    
    counts = defaultdict(int)
    chromosomes, genes, consequence = defaultdict(int), defaultdict(int), defaultdict(int)
    genes_skipped, genes_processed = set(), set()
    data = defaultdict(list)
    
    for v in tqdm(TSVReader(file)):
        counts['before'] += 1
        if not valid_consequence(v['Consequence']):
            counts['skip_consequence'] += 1
            continue
        if v['Feature'] not in TRANSCRIPTS:
            if v['Gene'] in GENES:
                # Selected gene without matching transcript id
                genes_skipped.add(v['Gene'])
            continue
        genes_processed.add(v['Gene'])
        # Remove multiple consequences
        # TODO shouldn't we get the first valid consequence?
        v['Consequence'] = v['Consequence'].split(',')[0]
        chromosome = v['Location'].split(":")[0]
        chromosomes[chromosome] += 1
        counts['after'] += 1
        consequence[v['Consequence']] += 1
        genes[v['SYMBOL']] += 1
        # Store all the lines matching a single position in a given patient
        # to get rid of positions that lie in more than one gene
        sample = v['#Uploaded_variation'].split('__')[1]
        data[(v['Location'], sample)].append(v)
    if counts['after'] == 0:
        print('There is no VEP output')
    synonymous_variants = consequence.get('synonymous_variant', 0)
    if synonymous_variants < 5:
        raise DatasetError(f'There are only {synonymous_variants} synonymous variants')
    if synonymous_variants == 0:
        ratio_missense = None
    else:
        ratio_missense = (consequence.get('missense_variant', 0) /
                          synonymous_variants)
    if ratio_missense is None:
        print("There are no synonymous variants")
    elif ratio_missense > 10:
        print("The ratio of missense / synonymous variants is too high")
        
    multiple_matches = {}
    for k, v in data.items():
        # If a (position , sample) (k) lie in more than one gene,
        # the first one in alphabetical order is selected and the
        # other are discarded
        if len(v) > 1:
            multiple_matches[k] = sorted([x['SYMBOL'] for x in v])
        yield sorted(v, key=lambda x: x['SYMBOL'])[0]
        
    stats['consequence'] = consequence
    stats['chromosomes'] = chromosomes
    stats['genes'] = genes
    stats['count'] = counts
    stats['ratio_missense'] = ratio_missense
    genes_orphan = genes_skipped - genes_processed
    
    if len(genes_orphan) > 0:
        stats['orphan_genes'] = list(genes_orphan)
        stats['warning_orphan_genes'] = f"There are {len(genes_orphan)} orphan genes"
    if len(chromosomes) < 23:
        stats['warning_few_chromosomes_with_mutations'] = f"There are only {len(chromosomes)} chromosomes with mutations"
    if len(multiple_matches) > 0:
        
        msg = ['{} in {} mapped to {}. Discarded {}'.format(k[0], k[1], v[0], ','.join(set(v[1:]))) for k, v in multiple_matches.items()]
        stats['warning_mutations_match_various_genes'] = msg
        
        
def process_VEP(infile, outfile):
    stats = {}
    try:
        with out_open(outfile, 'wt') as fd:
            writer = None
            for row in tqdm(filter_(infile, stats)):
                if writer is None:
                    fields = [f for f in row.keys()]
                    writer = csv.DictWriter(fd, fieldnames=fields, delimiter='\t')
                    writer.writeheader()
                writer.writerow(row)
    except DatasetError as e:
        # remove the file
        os.unlink(outfile)
        raise
    else:
        # write stats to file
        with open(outfile + '.stats.json', 'w') as fd:
            json.dump(stats, fd, indent=4)


finput = '/Users/picho/projects/collaborations/bbglab/CH/clonal_hematopoiesis_code/data_VAF/TCGA/vaf_filter/TCGA_FULL/TCGA_full.vep.sorted.tsv.out.gz'
fiout = '/Users/picho/projects/collaborations/bbglab/CH/clonal_hematopoiesis_code/data_VAF/TCGA/vaf_filter/TCGA_FULL/TCGA_full.vep.sorted.tsv.filtered.gz'

#process_VEP(finput, fiout)

finput = '/Users/picho/projects/collaborations/bbglab/CH/clonal_hematopoiesis_code/data_VAF/HMF/vaf_filter/HMF_FULL/HMF_full.vep.sorted.tsv.out.gz'
fiout = '/Users/picho/projects/collaborations/bbglab/CH/clonal_hematopoiesis_code/data_VAF/HMF/vaf_filter/HMF_FULL/HMF_full.vep.sorted.tsv.filtered.gz'
#process_VEP(finput, fiout)


df = pd.read_csv("/Users/picho/projects/collaborations/bbglab/CH/clonal_hematopoiesis_code/data_VAF/TCGA/vaf_filter/TCGA_FULL/TCGA_full.vep.sorted.tsv.filtered.gz", sep ='\t')
df['SAMPLE'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[1])
df['REF'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[2])
df['ALT'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[3])
df['POS'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[4])
df['VAF'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[5])
df['VARIANT_CLASS'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[6])
df['PHASING_PREDICT'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[7])
df['ID'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[8])
df['repeated_counts'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[9])

df.to_csv('/Users/picho/projects/collaborations/bbglab/CH/clonal_hematopoiesis_code/data_VAF/TCGA/vaf_filter/TCGA_FULL/TCGA_full.vep.sorted.tsv.filtered.toUSE.gz', 
         sep ='\t', index = False, header = True, compression = 'gzip')


df = pd.read_csv("/Users/picho/projects/collaborations/bbglab/CH/clonal_hematopoiesis_code/data_VAF/HMF/vaf_filter/HMF_FULL/HMF_full.vep.sorted.tsv.filtered.gz", sep ='\t')
df['SAMPLE'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[1])
df['REF'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[2])
df['ALT'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[3])
df['POS'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[4])
df['VAF'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[5])
df['VARIANT_CLASS'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[6])
df['PHASING_PREDICT'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[7])
df['ID'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[8])
df['FOUND_GERMLINE'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[9])
df['FOUND_SOMATIC'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[10])
df['MUTECT'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[11])
df['PHASING'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[12])
df['repeated_counts'] = df['#Uploaded_variation'].apply(lambda x : x.split('__')[13])
df.to_csv('/Users/picho/projects/collaborations/bbglab/CH/clonal_hematopoiesis_code/data_VAF/HMF/vaf_filter/HMF_FULL/HMF_full.vep.sorted.tsv.filtered.toUSE.gz', 
         sep ='\t', index = False, header = True, compression = 'gzip')