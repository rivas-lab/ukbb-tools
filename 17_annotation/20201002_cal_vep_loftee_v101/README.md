# The genotyped variant (`cal` dataset) annotation using VEP v101 with loftee plugin

## Data location

- `/oak/stanford/groups/mrivas/ukbb24983/cal/annotation_2021002`
- `/scratch/groups/mrivas/ukbb24983/cal/annotation_20201002`

## Annotation with `VEP` with loftee plugin

We [installed VEP with loftee plugin via Docker/Singularity image](https://github.com/rivas-lab/sherlock-modules/tree/master/vep).

Here, [`vep_example_single_thread.sh`](vep_example_single_thread.sh) illustrates the usage of the software.
We subsequently call the vcf to table conversion script, `/oak/stanford/groups/mrivas/software/loftee/src/tableize_vcf.py`.

## change log

- 2020/10/2: VEP version 101 with loftee. In September run, we did not use loftee plugin, but now it's incorporated into the pipeline.
- 2020/9/12: VEP version 101 with the latest data.
