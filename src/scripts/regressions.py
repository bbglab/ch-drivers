import statsmodels.api as sm
import pandas as pd
import numpy as np
from collections import defaultdict
import pandas as pd

def get_age_patients(df):
    
    try:
        bio = df['biopsyDate'].split('-')[0]
        birth = df['birthYear']
        age = int(bio)-int(birth)
    except:
        age = 'unknown'
    
    return age

def return_age_treatment_dict():
 
    pat = pd.read_csv(path_metadata, sep='\t', encoding='latin-1')

    pat['primaryTumorLocation'] = pat['primaryTumorLocation'].replace('Bone/soft tissue', 'Bone/Soft tissue')

    # fix primary location
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation'].apply(
        lambda x: str(x).replace(' ', '-').replace('/', '-')
    )
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('Head-and-Neck', 'Head-and-neck')
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('nan', 'Unknown')
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('CUP', 'Unknown')
    no_pretreat = pat[(pat['hasSystemicPreTreatment']=='No')&(pat['hasRadiotherapyPreTreatment']=='No')]
    no_pretreat_sample = no_pretreat['sampleId'].tolist()
    pretreat_sample = pat[(pat['hasSystemicPreTreatment']=='Yes')|(pat['hasRadiotherapyPreTreatment']=='Yes')]['sampleId'].tolist()

    pat['age'] = pat.apply(get_age_patients, axis = 1)
    dic_aging = dict(zip(pat['sampleId'], pat['age']))
    df_treat = pd.read_csv(path_treatments, sep ='\t', encoding='latin-1')

    return no_pretreat_sample, pretreat_sample, dic_aging

def do_logistic_reg(testing):

    all_cols = testing.columns
    
    testing.dropna(inplace = True)
    vals = testing[all_cols[:-1]]
    targ = testing[all_cols[-1]]
    vals = sm.add_constant(vals)
    logit_model=sm.Logit(targ, vals, )
    result=logit_model.fit(disp=0, method='bfgs')
    tup = result.params, result.conf_int()[0], result.conf_int()[1], result.pvalues

    return result, tup


def clinical_annot():

    df_treat = pd.read_csv(path_treatments, sep ='\t', encoding='latin-1')

    pat = pd.read_csv(path_metadata, sep='\t', encoding='latin-1')

    list_radiotherapy = pat[pat['hasRadiotherapyPreTreatment']=='Yes']['#patientId'].tolist()

    wanted = ['Alkaloid', 'Pyrimidine (ant)agonist', 'Topoisomerase inhibitor', 
            'Anthracycline', 'Alkylating','Platinum', 'Radionuclide', 'Taxane']
    
    pat['NEW_GENDER'] = pat['gender'].apply(lambda x : 1 if x == 'female' else 0 if x =='male' else 'not_found')

    togen = pat[pat['NEW_GENDER']!='not_found']
    dic_gender = dict(zip(togen['sampleId'], togen['NEW_GENDER']))

    wanted_t = ['Platinum', 'Pyrimidine (ant)agonist', 'Anthracycline', 'Topoisomerase inhibitor', 
            'Alkylating', 'Radionuclide',  'Alkaloid','Taxane',
            ]

    composed = [
                'Pyrimidine (ant)agonist, anthracycline, alkylating', 
                'Platinum, taxane', 'Radiofrequency ablation', 'Platinum, Pyrimidine (ant)agonist', 
                'Folate antagonist, Platinum', 'Taxane, anthracycline, alkylating',
                'Taxane, anthracycline, alkylating', 
                'Folinic acid, pyrimidine (ant)agonist, platinum', 'Anthracycline, alkylating', 
                'Alkylating, vinca alkaloid, alkylating, glucocorticoid',
                'Alkylating, antifolate, pyrimidine (ant)agonist', 
                'Folinic acid, pyrimidine (ant)agonist, topoisomerase inhibitor, platinum', 
                'Anti-CD20, alkylating, anthracycline, vinca alkaloid, glucocorticoid', 
                'Pyrimidine (ant)agonist, anthracycline, alkylating, taxane', 
                'Alkylating, alkylating, vinca alkaloid',
                'Anthracycline, antitumor antibiotic, vinca alkaloid, alkylating',
                'Vinca Alkaloid, antitumor antibiotic, alkylating',
                'Pyrimidine (ant)agonist, other',
                'Pyrimidine (ant)agonist, folinic acid',
                'Folate antagonist, Platinum', 
                'Folinic acid, pyrimidine (ant)agonist, topoisomerase inhibitor',
                'Vinca Alkaloid, alkaloid, glucocorticoid, anthracycline',
                'Glucocorticoid, pyrimidine (ant)agonist, platinum, alkylating, anthracycline, alkaloid',
                'Anti-mitotic vinca alkaloid',
                'Alkylating, alkaloid, pyrimidine (ant)agonist, alkylating', 
                'Taxane, Bisphosphonate','Microtubule inhibitor',
    ]

    dic_correction = {'Platinum':['Platinum', 
                        'Pyrimidine (ant)agonist, anthracycline, alkylating', 
                        'Platinum, taxane', 
                        'Platinum, Pyrimidine (ant)agonist', 
                        'Folate antagonist, Platinum', 
                        'Folinic acid, pyrimidine (ant)agonist, platinum', 
                        'Folate antagonist, Platinum', 
                        'Pyrimidine (ant)agonist, platinum',
                ],
                'Pyrimidine (ant)agonist':['Pyrimidine (ant)agonist',
                        'Pyrimidine (ant)agonist, platinum',
                        'Pyrimidine (ant)agonist, anthracycline, alkylating', 
                        'Platinum, Pyrimidine (ant)agonist', 
                        'Folinic acid, pyrimidine (ant)agonist, platinum', 
                        'Alkylating, antifolate, pyrimidine (ant)agonist', 
                        'Folinic acid, pyrimidine (ant)agonist, topoisomerase inhibitor, platinum', 
                        'Pyrimidine (ant)agonist, anthracycline, alkylating, taxane', 
                        'Pyrimidine (ant)agonist, other', 
                        'Pyrimidine (ant)agonist, folinic acid', 
                        'Folinic acid, pyrimidine (ant)agonist, topoisomerase inhibitor', 
                        'Glucocorticoid, pyrimidine (ant)agonist, platinum, alkylating, anthracycline, alkaloid', 
                        'Alkylating, alkaloid, pyrimidine (ant)agonist, alkylating'],

                'Anthracycline': ['Anthracycline',
                    'Pyrimidine (ant)agonist, anthracycline, alkylating', 
                    'Taxane, anthracycline, alkylating',
                    'Taxane, anthracycline, alkylating',
                    'Anthracycline, alkylating',
                    'Anti-CD20, alkylating, anthracycline, vinca alkaloid, glucocorticoid', 
                            'Pyrimidine (ant)agonist, anthracycline, alkylating, taxane', 
                    'Anthracycline, antitumor antibiotic, vinca alkaloid, alkylating',
                    'Vinca Alkaloid, alkaloid, glucocorticoid, anthracycline',
                    'Glucocorticoid, pyrimidine (ant)agonist, platinum, alkylating, anthracycline, alkaloid',

                    ],
                'Topoisomerase inhibitor':['Topoisomerase inhibitor', 

                        'Folinic acid, pyrimidine (ant)agonist, topoisomerase inhibitor, platinum', 
                        'Folinic acid, pyrimidine (ant)agonist, topoisomerase inhibitor',

                    ],

                'Alkylating':['Alkylating', 
                    'Pyrimidine (ant)agonist, anthracycline, alkylating', 
                    'Taxane, anthracycline, alkylating',
                    'Taxane, anthracycline, alkylating', 
                    'Folinic acid, pyrimidine (ant)agonist, platinum', 'Anthracycline, alkylating', 
                    'Alkylating, vinca alkaloid, alkylating, glucocorticoid',
                    'Alkylating, antifolate, pyrimidine (ant)agonist', 
                    'Anti-CD20, alkylating, anthracycline, vinca alkaloid, glucocorticoid', 
                    'Pyrimidine (ant)agonist, anthracycline, alkylating, taxane', 
                    'Alkylating, alkylating, vinca alkaloid',
                    'Anthracycline, antitumor antibiotic, vinca alkaloid, alkylating',
                    'Vinca Alkaloid, antitumor antibiotic, alkylating',
                    'Glucocorticoid, pyrimidine (ant)agonist, platinum, alkylating, anthracycline, alkaloid',
                    'Alkylating, alkaloid, pyrimidine (ant)agonist, alkylating',

                ], 

                'Radionuclide':['Radionuclide',
                    'Radiofrequency ablation',
                ], 

                'Taxane':['Taxane',
                    'Microtubule inhibitor', 
                    'Platinum, taxane',
                    'Taxane, anthracycline, alkylating',
                    'Pyrimidine (ant)agonist, anthracycline, alkylating, taxane', 
                    'Taxane, Bisphosphonate',

                ],

                'Alkaloid':['Alkaloid',
                    'Alkylating, vinca alkaloid, alkylating, glucocorticoid',
                    'Anti-CD20, alkylating, anthracycline, vinca alkaloid, glucocorticoid', 
                    'Alkylating, alkylating, vinca alkaloid',
                    'Anthracycline, antitumor antibiotic, vinca alkaloid, alkylating',
                    'Glucocorticoid, pyrimidine (ant)agonist, platinum, alkylating, anthracycline, alkaloid',
                    'Anti-mitotic vinca alkaloid',
                    ]
    }

    samples_with_treatment = {}

    allcytotreat = []
    for k, v in dic_correction.items():
        allcytotreat.extend(v)
        
    other_treatments = set()
    for sample, data in df_treat.groupby(by='#patientId'):
        for mechanism, l in data.groupby(by='mechanism'):
            if mechanism not in allcytotreat:
                other_treatments.add(sample)
                
    dic_treat = defaultdict(set)
    all_cytotoxic_s = set()
    for treat, l in dic_correction.items():
        for d in l:
            l_treat = df_treat[df_treat['mechanism']==d]['#patientId'].tolist()
            for s in l_treat:
                dic_treat[treat].add(s)
                all_cytotoxic_s.add(s)
    
    return all_cytotoxic_s, other_treatments, dic_treat, dic_gender

def gene_level_regression(testing, subs_myeloid_full_muts):

    l = ['ATM', 'ASXL1', 'TET2', 'DNMT3A', 'CHEK2', 'TP53', 'PPM1D']
    dic_CH_gene = {}
    allowed_samples = testing.index
    mut_s = {}
    for gene, data in subs_myeloid_full_muts.groupby(by='SYMBOL'):
        try:
            if gene in l:
                samples_affected = data['SAMPLE'].unique().tolist()
                wanted_s = [s for s in samples_affected if s in allowed_samples]
                mut_s[gene] = len(wanted_s)
                gene_pos = testing.loc[wanted_s]
                gene_neg = testing.loc[[s for s in testing.index if s not in samples_affected]]
                gene_pos['CH'] = 1
                gene_neg['CH'] = 0

                toreg = pd.concat([gene_pos, gene_neg])
                res, tup = do_logistic_reg(toreg)
                dic_CH_gene[gene] = tup
        except:
            continue

    return dic_CH_gene, mut_s

def gene_level_combined(testing, subs_myeloid_full_muts):
    
    dic_CH_gene = {}
    allowed_samples = testing.index

    l = ['CHEK2', 'TP53', 'PPM1D']
    data = subs_myeloid_full_muts[(subs_myeloid_full_muts['SYMBOL'].isin(l))]
    gene = 'Combination'
    if len(data)>0:

        samples_affected = data['SAMPLE'].unique().tolist()
        wanted_s = [s for s in samples_affected if s in allowed_samples]
        gene_pos = testing.loc[wanted_s]
        gene_neg = testing.loc[[s for s in testing.index if s not in samples_affected]]
        gene_pos['CH'] = 1
        gene_neg['CH'] = 0
        toreg = pd.concat([gene_pos, gene_neg])

        toreg = toreg.loc[:, (toreg != 0).any(axis=0)]

        res, tup = do_logistic_reg(toreg)
        dic_CH_gene[gene] = tup
    
    return dic_CH_gene

def main_regression(file, list_drivers):

    no_pretreat_sample, pretreat_sample, dic_aging = return_age_treatment_dict()

    df = pd.read_csv(file, sep ='\t')

    df['TREAT'] = df['SAMPLE'].apply(lambda x : 'YES' if x.split('_')[0] in pretreat_sample 
                                    else 'NO' if x.split('_')[0] in no_pretreat_sample else 'UNKNOWN')

    df['NEW_S'] = df['SAMPLE'].apply(lambda x : x.split('_')[0])
    df['AGE'] = df['NEW_S'].map(dic_aging)

    forbidden_samples = df[df['AGE'] =='unknown']['SAMPLE'].tolist()
    df = df[df['AGE'] !='unknown']
    dic_aging = dict(zip(df['SAMPLE'], df['AGE']))

    muts_in_drivers = df[df['SYMBOL'].isin(list_drivers)]
    muts_in_drivers = muts_in_drivers[~muts_in_drivers['Consequence'].isin(['synonymous_variant', 
                                                                        'inframe_deletion', 
                                                                        'inframe_insertion'])]
    subs_myeloid_full_muts = muts_in_drivers
    subs_myeloid_full = muts_in_drivers['SAMPLE'].unique()

    set_cytotoxic, other_treatments, dic_treat, dic_gender = clinical_annot()

    pat = pd.read_csv(path_metadata, sep='\t', encoding='latin-1')

    # samples with treatment information
    samples_without_treatment = pat[pat['hasSystemicPreTreatment'].isin(['No'])]['sampleId'].tolist()
    samples_radiotherapy = pat[(pat['hasRadiotherapyPreTreatment']=='Yes')]['sampleId'].tolist()
    samples_citotoxic_treatment = pat[(pat['hasSystemicPreTreatment']=='Yes')& (pat['#patientId'].isin(set_cytotoxic)) ]['sampleId'].tolist()
    samples_NOcitotoxic_treatment = pat[(pat['hasSystemicPreTreatment']=='Yes')& (pat['#patientId'].isin(other_treatments)) ]['sampleId'].tolist()

    samples_with_given_treat = set(samples_citotoxic_treatment + samples_NOcitotoxic_treatment + samples_radiotherapy)
    samples_to_check = set(samples_without_treatment + samples_citotoxic_treatment + samples_NOcitotoxic_treatment + samples_radiotherapy)

    alls = df['SAMPLE'].unique()

    unique_samples_treatment_info = [s for s in alls if s.split('_')[0] in samples_to_check]

    median_age = int(np.median([x for x in list(dic_aging.values()) if str(x)!='unknown']) )

    dic_to_reg = defaultdict(dict)

    # select data for treated samples
    for sample in unique_samples_treatment_info:
        
        age  = int(str(dic_aging.get(sample, 0)).replace('unknown', str(median_age)))
        dic_to_reg[sample]['Age'] = age/15
        dic_to_reg[sample]['TreatmentCITO'] = 0
        dic_to_reg[sample]['TreatmentNOCITO'] = 0
        dic_to_reg[sample]['Radiotherapy'] = 0

        dic_to_reg[sample]['Gender'] = dic_gender[sample.split('_')[0]]

        if sample.split('_')[0] in samples_citotoxic_treatment:
            dic_to_reg[sample]['TreatmentCITO'] = 1

        if sample.split('_')[0] in samples_NOcitotoxic_treatment:
            dic_to_reg[sample]['TreatmentNOCITO'] = 1
        
        if sample.split('_')[0] in samples_radiotherapy:
            dic_to_reg[sample]['Radiotherapy'] = 1
            
        dic_to_reg[sample]['CH'] = 0

        if sample in subs_myeloid_full : 
            dic_to_reg[sample]['CH'] = 1
        
    testing = pd.DataFrame(dic_to_reg).T
    all_cols = testing.columns
    res1, tup1 = do_logistic_reg(testing)

    #------------------- 
    # Second regression
    #-------------------

    dic_to_reg = defaultdict(dict)
    treat_set = set()
    good_sample_treated = set()

    # loop over all treated samples with cytotixic!
    for sample in unique_samples_treatment_info:  
        
        sample_t = sample.split('_')[0]
        sample_w = sample.split('T_')[0]
        
        #if (sample_t in samples_citotoxic_treatment):
        for treat in dic_treat:
            if sample_w in dic_treat[treat]:
                dic_to_reg[sample][treat] = 1
                good_sample_treated.add(sample)
            else:
                dic_to_reg[sample][treat] = 0

        age  = int(str(dic_aging.get(sample, 0)).replace('unknown', str(median_age)))
        dic_to_reg[sample]['Age'] = age/15
        dic_to_reg[sample]['Gender'] = dic_gender[sample.split('_')[0]]

        dic_to_reg[sample]['OtherTreatment'] = 0

        if sample_t in samples_NOcitotoxic_treatment:
            dic_to_reg[sample]['OtherTreatment'] = 1

        dic_to_reg[sample]['Radiotherapy'] = 0
        if sample_t in samples_radiotherapy:
            dic_to_reg[sample]['Radiotherapy'] = 1       

        dic_to_reg[sample]['CH'] = 0
        if sample in subs_myeloid_full : 
            dic_to_reg[sample]['CH'] = 1

    testing_full = pd.DataFrame(dic_to_reg).T

    testing = testing_full
    all_cols = testing.columns
    testing = testing[testing['Age']>0]
    res, tup = do_logistic_reg(testing)

    res_tt = defaultdict(dict)

    const, low_i, high_i, pval = tup 
    for c in const.keys():
        if c!='const':
            res_tt[c]= np.exp(const[c]), np.exp(low_i[c]), np.exp(high_i[c]), pval[c]
            
    tt_df = pd.DataFrame(res_tt, index = ['Val', 'LC', 'HC', 'pvalue']).T
    dic_CH_gene, mut_s = gene_level_regression(testing, subs_myeloid_full_muts)
    dic_CH_gene2 = gene_level_combined(testing, subs_myeloid_full_muts)

    return tup1, dic_CH_gene, mut_s, dic_CH_gene2

path_metadata = '/Users/picho/projects/collaborations/bbglab/CH/cluster/metadata.tsv'
path_treatments = '/Users/picho/projects/collaborations/bbglab/CH/cluster/pre_biopsy_drugs_by_patient.tsv'