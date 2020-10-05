# LD map for the `cal + hla + imp` dataset

Yosuke Tanigawa

## Version history

- 2019 Sept: we used the updated variant filtering criteria for the imputed variants (we don't use HWE filter anymore).
  - Starting this version, we use the "array_imp_combined" dataset instead of using "array_imp_no_cnv" dataset simply because the creation of the new PLINK file is not required. We can simply specify the set of variants we'd like to use (i.e. array + HLA + imputed) and pass it to plink with `--extract`.
- 2020 March: rerunning the analysis with `--ld-window 100000000`
- 2019 Dec.: the initial attempt. It turend out that there was a bug in the pipe (we forgot to specify `--ld-window` arg) and the resulting LD map is incomplete.

## File location

We have three directories corresponding to the three versions.

- `$OAK` space `/oak/stanford/groups/mrivas/ukbb24983/array_imp_combined_no_cnv/ldmap/`
  - `ldmap_20200928`
  - `ldmap_20200304`
  - `ldmap_20191216`
- [Google Drive](https://drive.google.com/drive/folders/1Jq62lwncX37zaAArLA1n8D1-PxHnkhpn)

Please see [here](https://github.com/rivas-lab/ukbb-tools/tree/master/14_LD_map) for types of the output files.

## Copy of the data in GBE (based on Dec. 2019)

- `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser/static/ldmap`

## Public release

WB LD map (version March 2020) was released to public: https://bit.ly/rivas-lab_covid19_UKB_LD_public_release
