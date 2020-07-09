# Exome variant annotation quick fix

2020/7/8

This is a partial fix to https://github.com/rivas-lab/ukbb-tools/issues/29.

We joined the pvar file, allele frequency (computed on WB), and a subset of variant annotation.
The resulting file has the following annotation files.

- `CHROM`: chromosome
- `POS`: position
- `ID`: variant ID
- `REF`: the reference allele
- `ALT`: the alternate allele
- `ALT_FREQ_white_british`: the allele frequency of the alternate allele in the white british population
- `rsID`: the annotated variant ID. This column can be empty or NA.
- `Gene`: the annotated Ensembl gene ID. This column can be empty or NA.
- `Gene_symbol`: the annotated gene symbol. This column can be empty or NA.
- `Consequence`: the predicted consequence of the variant.

The resulting file is written to `/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb_2020/ukb_exm_spb_2020.pvar.zst` and the original file was moved to `/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb_2020/ukb_exm_spb_2020.pvar.original.zst`
