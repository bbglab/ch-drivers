This requires you to previously run INTOGEN for HMF, TCGA and MSKCC CH datasets (the last one won't be the full INTOGEN pipeline, please look at the Methods section of our paper). The scripts should be run in the following order:

` python driver_CH_list.py
`

` python impact_discovery.py
`

`python create_discovery.py
`

`python web_tables.py
`
