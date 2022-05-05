from bgparsers import readers
import gzip
import subprocess

# this is adapted from Loris code
def process_to_TCGA_vep(f):
    outf = f.replace('tsv.gz', 'vep.tsv.gz')
    with gzip.open(outf, 'wt') as outfile:
        i = 0
        for m in readers.variants(f, 
                                  extra=['VAF', 'VARIANT_CLASS', 'PHASING_PREDICTED',
                                         'ID', 'count_repeat', 
                                         ]): 
            
            start = end = int(m["POSITION"])
            if m["REF"] == "-":  # is an insertion
                start = end + 1
            elif m["ALT"] == "-":  # is a deletion
                end = start + len(m["REF"]) - 1
            else:  # snv/mnv
                end = start + len(m["ALT"]) - 1
            fields = [
                m['CHROMOSOME'],
                f"{start}",
                f"{end}",
                f"{m['REF']}/{m['ALT']}",
                m['STRAND'],
                f"I{i:010d}__{m['SAMPLE']}__{m['REF']}__{m['ALT']}__{m['POSITION']}__{m['VAF']}__{m['VARIANT_CLASS']}__{m['PHASING_PREDICTED']}__{m['ID']}__{m['count_repeat']}", 
            ]
            i+=1
            out_f = '\t'.join(fields)
            outfile.write(out_f+'\n')    

def process_to_HMF_vep(f):
    outf = f.replace('tsv.gz', 'vep.tsv.gz')
    with gzip.open(outf, 'wt') as outfile:
        i = 0
        for m in readers.variants(f, extra=['VAF', 'PHASING', 'VARIANT_CLASS', 'PHASING_PREDICTED', 'ID', 
                                "FOUND_GERMLINE","FOUND_SOMATIC", 'MUTECT', 'count_repeat', "minorCN"]): 

            start = end = int(m["POSITION"])
            if m["REF"] == "-":  # is an insertion
                start = end + 1
            elif m["ALT"] == "-":  # is a deletion
                end = start + len(m["REF"]) - 1
            else:  # snv/mnv
                end = start + len(m["ALT"]) - 1
            fields = [
                m['CHROMOSOME'],
                f"{start}",
                f"{end}",
                f"{m['REF']}/{m['ALT']}",
                m['STRAND'],
                f"I{i:010d}__{m['SAMPLE']}__{m['REF']}__{m['ALT']}__{m['POSITION']}__{m['VAF']}__{m['VARIANT_CLASS']}__{m['PHASING_PREDICTED']}__{m['ID']}__{m['FOUND_GERMLINE']}__{m['FOUND_SOMATIC']}__{m['MUTECT']}__{m['PHASING']}__{m['count_repeat']}__{m['minorCN']}", 
            ]

            i+=1
            out_f = '\t'.join(fields)
            outfile.write(out_f+'\n')
    
# as always, change your paths accordingly
path_tcga = '/Users/picho/projects/collaborations/bbglab/CH/clonal_hematopoiesis_code/data_VAF/TCGA/vaf_filter/'
f = '{}/TCGA_FULL/TCGA_full.tsv.gz'.format(path_tcga)
process_to_TCGA_vep(f)
cmd = "gunzip -c {0}/TCGA_FULL/TCGA_full.vep.tsv.gz | sort -k1V -k2n -k3n\
| gzip > {0}/TCGA_FULL/TCGA_full.vep.sorted.tsv.gz".format(path_tcga)
subprocess.check_output(cmd, shell =  True)

path_hmf = '/Users/picho/projects/collaborations/bbglab/CH/clonal_hematopoiesis_code/data_VAF/HMF/vaf_filter/'
f = '{}/HMF_FULL/HMF_full.tsv.gz'.format(path_hmf)
process_to_HMF_vep(f)
cmd = "gunzip -c {0}/HMF_FULL/HMF_full.vep.tsv.gz | sort -k1V -k2n -k3n\
| gzip > {0}/HMF_FULL/HMF_full.vep.sorted.tsv.gz".format(path_hmf)
subprocess.check_output(cmd, shell = True)



