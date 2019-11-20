# Bulk data download

Some of the data types are categorized as "Bulk" and does not come with the table file.

The instruction for bulk download can be found here:
http://biobank.ctsu.ox.ac.uk/showcase/docs/ukbfetch_instruct.html

## Data location

Please check: https://github.com/rivas-lab/ukbb24983wiki/tree/master/bulk

## Method (How to download the data)

### Step 0. Identify the UKB24983 table that contains the field of your interest

```
$ ml load ukbb-query
$ ukbb-query_find_table_by_field_id.sh 20204
http://bit.ly/UKB24983_tables
ukb37855        8735    f.20204.2.0     f       20204   2       0
ukb35059        116     f.20204.2.0     f       20204   2       0
ukb28983        78      f.20204.2.0     f       20204   2       0
ukb25826        4282    f.20204.2.0     f       20204   2       0
ukb21731        4281    f.20204.2.0     f       20204   2       0
ukb10137        4281    f.20204.2.0     f       20204   2       0
```

From this results and the information on http://bit.ly/UKB24983_tables, we found
ukb37855 is the most recent table that has the field.

### Step 1. Apply `ukbconv` to get list of individuals with bulk file

Using `ukbconv` program (in `ukbb-showcase-utils` module), you can generate the list of bulk files.

Example:
```
$ ml load ukbb-showcase-utils
$ cd /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/2005693/37855/download
$ ukbconv ukb37855.enc_ukb bulk -s20204
$ mv ukb37855.bulk ukb37855.20204.bulk # we highly recommend to rename the output files (before overwritten by another run).
```

[The section 5.1 of the official document](http://biobank.ctsu.ox.ac.uk/showcase/docs/ukbfetch_instruct.html) explains how to get the bulk list.

### Step 2. Place `.key` file

Ask PI and get the key file and place the key file as `.ukbkey` in the download directory.

### Step 3. bulk download wrapper script

We have a custom wrapper script for batch download. 
https://github.com/rivas-lab/ukbb-tools/blob/master/08_bulk_DL/ukbfetch_bulk_wrapper.sh

#### example: how to submit a download job

```
[ytanigaw@sh-109-53 /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/2005693/37855/bulk/ukb2005693.37855.20204]$ cat dl20204.sbatch
#!/bin/bash
#SBATCH --job-name=dl20204
#SBATCH --output=dl20204.%A.out
#SBATCH  --error=dl20204.%A.err
#SBATCH --nodes=1
#SBATCH --cores=4
#SBATCH --mem=30000
#SBATCH --time=7-00:00:00
#SBATCH -p mrivas,normal

set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

out_d="/scratch/groups/mrivas/ukbb/24983/bulk/ukb2005693.37855.20204"
src="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/08_bulk_DL/ukbfetch_bulk_wrapper.sh"

bash ${src} ${out_d}/ukb37855.20204.bulk ${out_d} ${out_d}/k24983.key ${cores}

[ytanigaw@sh-109-53 /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/2005693/37855/bulk/ukb2005693.37855.20204]$ sbatch dl20204.sbatch
Submitted batch job 53180525
```

#### usage of the wrapper script

```
[ytanigaw@sh-102-07 /scratch/groups/mrivas/ukbb/24983/bulk/ukb24615.21015]$ bash /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/08_bulk_DL/ukbfetch_bulk_wrapper.sh ukb24615.21015.bulk .
```

#### Why we wrote a custom wrapper script?
In [the official document](http://biobank.ctsu.ox.ac.uk/showcase/docs/ukbfetch_instruct.html), there is a section called "5.1 Automation Example". They say `ukbfetch` cannot take batch file with more than 1000 lines. They illustrate some examples to use `-s` and `-n` flags to mitigate this limitation, but this functionality is broken (as of 2019/2/27). I (Yosuke) wrote a wrapper script that can download bulk file with more than 1000 files in parallel.
