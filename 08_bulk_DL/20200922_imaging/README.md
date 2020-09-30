# imaging dataset download
## 2020/9/22 Yosuke Tanigawa

| UKB Field ID                                                     | Description                             | # data files | # individuals |
|------------------------------------------------------------------|-----------------------------------------|--------------|---------------|
| [20202](http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20202) | Pancreatic fat - DICOM                  | 44108 *      | 41754         |
| [20203](http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20203) | Liver images - gradient echo - DICOM    | 10108        | 10108         |
| [20204](http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20204) | Liver Imaging - T1 ShMoLLI - DICOM      | 43267        | 45685         |
| [20254](http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20254) | Liver images - IDEAL protocol - DICOM   | 34839        | 37267         |
| [20260](http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20260) | Pancreas Images - gradient echo - DICOM | 37751        | 35415         |

* We were not able to download one file. Please see the issue ticket [#38](https://github.com/rivas-lab/ukbb-tools/issues/38).

Thoese data fields are included in `ukb41413` (releae date: 2020/4/1). We used the table file from this release to identify the list of individuals for bulk download.

## Data location

- data directory: `/scratch/groups/mrivas/ukbb24983/phenotypedata/2005693/41413/bulk`
- For each UKB Field ID, we have a separate sub-directory: `ukb2005693.41413.<UKB Field ID>`

Please ignore (but pleast don't delete) `*.bulk` and `*.lis` files - they are the index files used in the bulk download script.

## scripts

- [`1_ukbconv.sh`](1_ukbconv.sh): this script generates the list of individuals with the specified bulk field
- [`2_bulk_dl.sh`](2_bulk_dl.sh): this SLURM job script download the bulk data
