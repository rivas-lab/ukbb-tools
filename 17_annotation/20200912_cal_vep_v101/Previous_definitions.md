# The older rules to aggregate Consequence fields to Csq groups

## Biomarker project etc.

```{R}
variant_anno_f %>%
fread(select=c('#CHROM', 'POS', 'REF', 'ALT', 'ID', 'Consequence', 'ld_indep', 'Gene_symbol'), colClasses='character') %>%
rename('CHROM'='#CHROM') -> annot.arr

annot.arr$Csq[
    !(annot.arr$Consequence %in% c("frameshift_variant","splice_donor_variant","stop_gained","stop_lost","start_lost","splice_acceptor_variant","splice_region_variant","missense_variant","inframe_insertion","inframe_deletion"))
] = "nc"
annot.arr$Csq[
    annot.arr$Consequence %in% c("splice_region_variant","missense_variant","inframe_insertion","inframe_deletion")
] = "PAVs"
annot.arr$Csq[
    annot.arr$Consequence %in% c("frameshift_variant","splice_donor_variant","stop_gained","stop_lost","start_lost","splice_acceptor_variant")
] = "PTVs"

```

## Guhan's MRP criteria

```{Python}
    ptv = [
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "start_lost",
        "stop_lost",
    ]
    pav = [
        "protein_altering_variant",
        "inframe_deletion",
        "inframe_insertion",
        "splice_region_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "missense_variant",
    ]
    pcv = [
        "synonymous_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "coding_sequence_variant",
        "incomplete_terminal_codon_variant",
        "TF_binding_site_variant",
    ]
    intron = [
        "regulatory_region_variant",
        "intron_variant",
        "intergenic_variant",
        "downstream_gene_variant",
        "mature_miRNA_variant",
        "non_coding_transcript_exon_variant",
        "upstream_gene_variant",
        "NA",
        "NMD_transcript_variant",
    ]
```