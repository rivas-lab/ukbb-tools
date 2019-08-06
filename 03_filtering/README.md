# Filtering

We do not employ sample- and variant-level QC in the traditional sense in this pipeline. Sample-level QC is often done by the BioBank itself before the data is sent out, and variant-level QC is done after GWAS is run (in case any signal can be recovered despite bad data).

The notebooks in this folder do two distinct tasks:

1) [Go over the overall marker quality of different BioBank arrays](https://github.com/rivas-lab/ukbb-tools/blob/master/03_filtering/Marker_QC.ipynb) (marker_QC.ipynb)
2) [Redact individuals from the dataset and define populations](https://github.com/rivas-lab/ukbb-tools/blob/master/03_filtering/sample_qc_v3.ipynb) (rsample_qc_v3.ipynb)
