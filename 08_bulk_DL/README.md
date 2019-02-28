# Bulk data download

Some of the data types are categorized as "Bulk" and does not come with the table file.

The instruction for bulk download can be found here:
http://biobank.ctsu.ox.ac.uk/showcase/docs/ukbfetch_instruct.html


## `ukbconv` to get list of individuals with bulk file

The section 5.1 explains how to get the bulk list.


Example:
```
$ ukbconv  ukb789.enc_ukb  bulk  -s145
```


## bulk download wrapper script

In the official documentation, there is a section called "5.1 Automation Example". 
They say `ukbfetch` cannot take batch file with more than 1000 lines. They illustrate some examples to use `-s` and `-n` flags to mitigate this limitation, but this functionality is broken (as of 2019/2/27). I (Yosuke) wrote a wrapper script that can download bulk file with more than 1000 files in parallel.


### usage of the wrapper script

```
[ytanigaw@sh-102-07 /scratch/groups/mrivas/ukbb/24983/bulk/ukb24615.21015]$ bash /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/08_bulk_DL/ukbfetch_bulk_wrapper.sh ukb24615.21015.bulk .
```
