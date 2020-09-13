## Future development goals

In this version, we used Docker/Singularity version of VEP.

Compared to [what we had before](https://github.com/rivas-lab/ukbb-tools/blob/c95a924ae30bb33d8f188aea061a68accb10b126/17_annotation/annotate_bims.sh#L53), this approach does not incorporate the results from `LoF` plugin.

The `loftee` plugin has some dependencies on Perl modules. A proper way to handle this is to write a Docker file specifying those dependencies.

I (Yosuke) tried to pass the `PERL5LIB` variable, only to find that there are additional dependencies.


```{bash}
SINGULARITYENV_PERL5LIB="/oak/stanford/groups/mrivas/users/ytanigaw/repos/konradjk/loftee:/opt/vep/src/ensembl-vep:/opt/vep/src/ensembl-vep/modules" singularity run /oak/stanford/groups/mrivas/software/vep/v101/ensembl-vep_release_101.0.sif \

--plugin LoF,loftee_path:/oak/stanford/groups/mrivas/users/ytanigaw/repos/konradjk/loftee,human_ancestor_fa:${public_d}/loftee_human_ancestor_20170411/human_ancestor.fa.gz,conservation_file:${public_d}/loftee/phylocsf_gerp.sql \
```

Also, `--force_overwrite` was used in the original version, but we decided to drop them.
