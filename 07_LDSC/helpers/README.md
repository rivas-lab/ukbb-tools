# LDSC wrapper functions

The input file needs to be in the LDSC format. For that, we have a conversion script called [`ukb_ldsc_munge_wrapper.sh`](ukb_ldsc_munge_wrapper.sh).

Once input files are converted, then you can call [`ukb_ldsc_rg_helper.sh`](ukb_ldsc_rg_helper.sh) for the analysis.

- Note: both script calls Python script.
  - Yosuke: he has a conda environment called `ldsc_ytanigaw` (`$ conda activate ldsc_ytanigaw`).
  - The LDSC repo is cloned (and patched?) in `/oak/stanford/groups/mrivas/users/ytanigaw/repos/bulik/ldsc`.
  - For others, it should be possible to install it via `$ conda env create --file /oak/stanford/groups/mrivas/users/ytanigaw/repos/bulik/ldsc/environment.ytanigaw.yml`
