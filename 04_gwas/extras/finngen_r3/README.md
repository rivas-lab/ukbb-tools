# FinnGen R3 public release

- https://www.finngen.fi/en/access_results

## files

finngen_r3_variants.hg19.unmapped.gz
finngen_r3_variants.master.tsv.gz
finngen_r3_variants.master.tsv.gz.tbi

## scripts

- [`1_download.sh`](1_download.sh): download the sumstats
- [`2_liftOver.sh`](2_liftOver.sh): apply liftOver


## instruction

```
---------- Forwarded message ---------
From: <humgen-servicedesk@helsinki.fi>
Subject: Registration for FinnGen GWAS Summary Statistics download


Dear researcher,

Thanks for your interest in FinnGen data. 
Below you can find the information on how to download the data.

Released FinnGen GWAS summary statistics can be downloaded from Google cloud storage free of charge.
______________________________________________________________   

INSTRUCTIONS FOR WEB BROWSER-BASED ACCESS:
1) Open web browser (Google Chrome is recommended)
2) Navigate: 
https://console.cloud.google.com/storage/browser/finngen-public-data-r3/summary_stats/
https://console.cloud.google.com/storage/browser/finngen-public-data-r3/finemapping/
or
https://console.cloud.google.com/storage/browser/finngen-public-data-r2/summary_stats/

3) Login with your google account
4) Select the files to be downloaded
5) Use ... at the right-hand side to start downloading  

______________________________________________________________

INSTRUCTIONS FOR COMMAND-LINE ACCESS:
Using wget utility  https://www.gnu.org/software/wget
Example:
wget https://storage.googleapis.com/finngen-public-data-r3/summary_stats/finngen_r3_AB1_ARTHROPOD.gz

Using curl utility https://curl.haxx.se/docs/
Example:
curl https://storage.googleapis.com/finngen-public-data-r3/summary_stats/finngen_r3_AB1_ARTHROPOD.gz -o finngen_r3_AB1_ARTHROPOD.gz

______________________________________________________________

INSTRUCTIONS FOR GOOGLE CLOUD-BASED ACCESS:
To install Google Cloud SDK follow directions https://cloud.google.com/sdk/install
1) List the files
gsutil ls gs://finngen-public-data-r3/summary_stats/

2) Copy the files
gsutil cp gs://finngen-public-data-r3/summary_stats/finngen_r3_AB1_ARTHROPOD.gz /path/to/your/incoming_folder/

______________________________________________________________

FURTHER INFORMATION:
The Manifest file with the link to all the downloadable summary statistics is available at: https://storage.googleapis.com/finngen-public-data-r3/summary_stats/r3_manifest.tsv

Linkage disequilibrium (LD) estimations data based on Finnish SISU panel v3 is available at:
https://console.cloud.google.com/storage/browser/finngen-public-data-ld/imputation_panel_v1/

Guidelines how to use ".bcor" LD data files is available:
https://finngen.gitbook.io/documentation/methods/genotype-imputation/ld-estimation#example-usage

More information about the data QC, PheWAS methodology can be obtained at: https://finngen.gitbook.io/documentation/



Explore the results at: http://r3.finngen.fi/ or http://r2.finngen.fi/ 

______________________________________________________________
Please consider visiting the study website(https://www.finngen.fi/en) and follow FinnGen on twitter: @FinnGen_FI(https://twitter.com/finngen_fi?lang=en)

If you want to host the FinnGen summary statistics on your website, please get in contact with us at: humgen-servicedesk@helsinki.fi
```
