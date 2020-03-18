# liftOver and flip fix

The scripts in this directory provide liftOver & flip fix functions.
Those scripts are available as part of `ukbb-tools` module. Please see the usage below.
You can, of course, use the latest scripts by directly calling the bash scripts in this directory.

## load module

```{bash}
ml load ukbb-tools
```

## liftOver

### `$liftOver_wrapper_sh` script

When you load the module, you can call `$liftOver_wrapper_sh`. Please note that you need to type `$` because the path to the script is stored in this environmental variable.

```{bash}
$liftOver_wrapper_sh
liftOver_wrapper.sh: 3 positional arguments are required
liftOver_wrapper.sh (version 0.1.0)
Apply UCSC liftOver. If the input file contains OR or BETA column, the script will also automatically apply the flipfix.
  Author: Yosuke Tanigawa (ytanigaw@stanford.edu)

Usage: liftOver_wrapper.sh [options] in_file out_mapped out_unmapped
  in_file       The input file. Please see the notes blow about our assumptions on the input file format.
  out_mapped    The output file of the mapped elements.
  out_unmapped  The output file of the unmapped elements.

Options:
  --threads  (-t) Number of CPU cores
  --src_genome    The genome build for the input file
  --dst_genome    The genome build for the output file

Notes:
  - We assume the input file has a header line (first line) that starts with .
  - The input file needs to have the following columns: CHROM, POS, REF, and ID.
  - If effect size column (BETA or OR) is provided, we will automatically apply flipfix.
  - You may include additional columns.

Default configurations (please use the options above to modify them):
  to_bed_field_sep=!
  bed_chr_prefix=chr
  threads=4
  src_genome="hg19"
  dst_genome="hg38"
  ```

### Example usage (liftOver)

```{bash}
ml load ukbb-tools

data_dir=/oak/stanford/groups/mrivas/projects/hearing/finngen
p=H8_CONSENHEARINGLOSS

$liftOver_wrapper_sh \
--src_genome hg38 --dst_genome hg19 \
${data_dir}/summary_stats_plink/${p}.gz \
${data_dir}/summary_stats_plink_hg19/${p}.gz \
${data_dir}/summary_stats_plink_hg19/${p}.unmapped.txt.gz
```

## flip fix

```{bash}
$flipfix_sh
flipfix.sh: 1 positional arguments are required
flipfix.sh (version 1.1.1)
Apply flipfix.
  Author: Yosuke Tanigawa (ytanigaw@stanford.edu)

Usage: flipfix.sh [options] in_file
  in_file       The input file. Please see the notes blow about our assumptions on the input file format.

Options:
  --assembly    The genome build for the input file
  --ref_fa      The reference genome sequence. It it's specified as AUTO, it will be automatically grab one from /scratch/groups/mrivas/public_data/genomes (the default behavior).

Notes:
  - We assume the input file has a header line (first line) that starts with .
  - The input file needs to have the following columns: CHROM, POS, ID, REF, ALT, A1, and one of the following: OR or BETA.
  - You may include additional columns.

Default configurations (please use the options above to modify them):
  to_bed_field_sep=!
  bed_chr_prefix=chr
  ref_fa="AUTO"
  assembly="hg19"
  ```

### Example usage (flipfix)

```{bash}
ml load ukbb-tools

[ytanigaw@sh-109-53 /oak/stanford/groups/mrivas/private_data/summary_stats/BBJ_Masa_20191020]$ $flipfix_sh plink_format/BBJ.Hearing_Loss.gz | bgzip -l9 --threads 4 > plink_format_hg19/BBJ.Hearing_Loss.gz
```
