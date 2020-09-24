# GWAS plot scripts
## Yosuke Tanigawa (2020/9/23)

We have two scripts and a helper function:

- [`gwas_plot_qq.R`](gwas_plot_qq.R): this script generates GWAS qq plot
- [`gwas_plot_manhattan.R`](gwas_plot_manhattan.R): this script generates Manhattan plot
- [`gwas_plot_misc.R`](gwas_plot_misc.R): this file contains helpfer functions for the two scripts

## Usage

```{bash}
Rscript gwas_plot_qq.R <summary.statistics.tsv.gz> qq.png

Rscript gwas_plot_manhattan.R <summary.statistics.tsv.gz> <variant annotation file> manhattan.png
```

Please see the [`17_annotation`](/17_annotation#location-of-relevant-annotation-files-on-oak) for the location of the variant annotation file.

