# imaging dataset download
## 2020/9/22 Yosuke Tanigawa

- Liver images - gradient echo - DICOM (UKB Field: [20203](http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20203))
- Pancreatic fat - DICOM (UKB Field: [20202](http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20202))
- Pancreas Images - gradient echo - DICOM (UKB Field: [20260](http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20260))

Thoese data fields are included in `ukb41413` (releae date: 2020/4/1). We used the table file from this release to identify the list of individuals for bulk download.

## Data location

- data directory: `/scratch/groups/mrivas/ukbb24983/phenotypedata/2005693/41413/bulk`
- We have 3 sub-directories, each corresponds to one UKB Field ID.
  - ukb2005693.41413.20202
  - ukb2005693.41413.20203
  - ukb2005693.41413.20260

Please ignore (but pleast don't delete) `*.bulk` and `*.lis` files - they are the index files used in the bulk download script.

## scripts

- [`1_ukbconv.sh`](1_ukbconv.sh): this script generates the list of individuals with the specified bulk field
- [`2_bulk_dl.sh`](2_bulk_dl.sh): this SLURM job script download the bulk data
