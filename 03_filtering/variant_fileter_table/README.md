# variant filter table

This directory contains variant-level QC for the genotype dataset.

Please also check [`17_annotation`](/17_annotation) directory for the relevant information.

## [`snp_qc_patch.ipynb`](snp_qc_patch.ipynb)

- The [snp QC file provided by UK Biobank](http://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=1955) has marker-level QC information.
- It seems like the file has some errors in header line: `PC19_loading`, `PC29_loading`, and `PC39_loading` are all mislabed as `PC9_loading`
- This notebook applies patches to the issues above and merge it with the pvar file.

Data: `/oak/stanford/groups/mrivas/ukbb24983/snp/ukb_snp_qc.pvar.zst`
