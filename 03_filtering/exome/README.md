# Exome dataset

## data location

- Please see our wiki:
  - https://github.com/rivas-lab/wiki/wiki/UKBB-dataset#exome-dataset

## Exome spb dataset for 50k sample (2020)

We used [`ukb_exm_spb_2020.download.sh`](ukb_exm_spb_2020.download.sh) to download the exome SPB dataset. The log file is [`ukb_exm_spb_2020.download.20200611-093925.log`](ukb_exm_spb_2020.download.20200611-093925.log).

### flipcheck

We confirmed that there is no flip.

```{bash}
cd /oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb_2020

bash ~/repos/rivas-lab/ukbb-tools/09_liftOver/flipcheck.sh --assembly hg38 ukb_exm_spb_2020.pvar.zst | awk 'toupper($4) != toupper($6)'
#CHROM  POS     ID      REF     ALT     FASTA_REF
```
