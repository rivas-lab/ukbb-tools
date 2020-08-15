# assigining phenotype groups

Yosuke Tanigawa, 2020/8/12

## phenotype grouping

We provide phenotype grouping information on the PheWAS page on GBE. The notebook/scripts in this directory were used to geerate the phenotype grouping.

We have several groups of phenotypes that are defined in a specific way:

### Cancer phenotypes ("cancer" phenotype group in the PheWAS page)

### High confidence phenotype definitions ("Disease_outcome" phenotype group in the PheWAS page)

### Family history phenotype definitions ("Family_history" phenotype group in the PheWAS page)

### Biomarkers ("Biomarkers" phenotype group in the PheWAS page)

We applied technical covariate correction for 35 blood and urine biomarker phenotypes. Those phenotypes are listed under "Biomarkers" group and the detailed description on the phenotyping procedure is explained in the following manuscript:

N. Sinnott-Armstrong*, Y. Tanigawa*, et al, Genetics of 38 blood and urine biomarkers in the UK Biobank. bioRxiv, 660506 (2019). doi:10.1101/660506

### Other phenotypes from UK Biobank.

For the UK Biobank dataset hosted on Global Biobank Engine, we used the following fields to define the phenotypes.

Please check [`phenotype_description.html`](phenotype_description.html) for more description, which is draft texts for [the FAQ section on Global Biobank Engine](https://gbe.stanford.edu/faq).

The following notebooks were used for the curation of the phenotype grouping info:

- [`1_curate_UKB_category.ipynb`](1_curate_UKB_category.ipynb): in this notebook, we fetched [the Data Dictrionary Showcase table](http://biobank.ctsu.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.csv) from UKB, parsed the "path" field (that represents the tree structure in UKB field browser), applied some filtering using regular expression, and saved the results into a file.
  - [`UKB_fields_with_category.20200812.tsv`](UKB_fields_with_category.20200812.tsv): this is the results of the curated groupings of the UKB fields
- [`2_assign_GBE_category.ipynb`](2_assign_GBE_category.ipynb): using the curated groupings of the UKB fields, and also considering the special phenotype groups in GBE (`HC` disease definition, `cancer`, `FH` family history information, and `Biomarkers` phenotype), we provide the groupings on the `GBE_ID`s. The key difference of this procedure from the one in the previous step is that [`1_curate_UKB_category.ipynb`](1_curate_UKB_category.ipynb) provides the groupings on the UKB fields whereas [`2_assign_GBE_category.ipynb`](2_assign_GBE_category.ipynb) provides the groupings on the `GBE_ID`s.
  - [`GBE_category.20200812.tsv`](GBE_category.20200812.tsv): the resulting groupins for `GBE_ID`s as well as case counts across populations. Please see below for the column descriptor.
- [`3_icdinfo.ipynb`](3_icdinfo.ipynb): this notebook was used to generate `icdinfo` files for (7 + 1) population groups. The results are saved in [`icdinfo`](icdinfo) directory.

#### set of regular expression used in the curation of UKB field grouping

Some phenotype grouping provided in the UKB field browser is too specific.
We removed bottom-level groupings so that we have a moderate number of categories for UKB fields.

The followings are the set of custom regular expression as of 8/13/2020. Please check [`1_curate_UKB_category.ipynb`](1_curate_UKB_category.ipynb) for the latest set of RegEx.

```{R}
# here we have a custom Regex to clean-up the tree structure
    Path = str_replace(Path, 'Touchscreen > ([^>]+) > .*', 'Touchscreen > \\1'),
    Path = str_replace(Path, 'Online follow-up > ([^>]+) > .*', 'Online follow-up > \\1'),
    Path = str_replace(Path, 'Brain MRI > ([^>]+) > .*', 'Brain MRI > \\1'),
    Path = str_replace(Path, 'Health-related outcomes > ([^>]+) > .*', 'Health-related outcomes > \\1'),
    Path = str_replace(Path, 'Additional exposures > ([^>]+) > .*', 'Additional exposures > \\1'),
    Path = str_replace(Path, 'Genomics > ([^>]+) > .*', 'Genomics > \\1'),
    Path = str_replace(Path, 'Assay results > ([^>]+) > .*', 'Assay results > \\1'),
    Path = str_replace(Path, 'Physical measures > ([^>]+) > .*', 'Physical measures > \\1'),
    Path = str_replace(Path, '^Population characteristics >.*', 'Population characteristics'),
    Path = str_replace(Path, 'Sample inventory >.*', 'Sample inventory'),
    Path = str_replace(Path, 'Biological sampling >.*', 'Biological sampling'),
    Path = str_replace(Path, 'Recruitment >.*', 'Recruitment'),
    Path = str_replace(Path, 'Procedural metrics >.*', 'Procedural metrics'),
    Path = str_replace(Path, 'Cognitive function >.*', 'Cognitive function'),
    Path = str_replace(Path, 'Cognitive function online', 'Cognitive function'),
    Path = str_replace(Path, 'Abdominal MRI >.*', 'Abdominal MRI'),
    Path = str_replace(Path, 'Heart MRI >.*', 'Heart MRI'),
    Path = str_replace(Path, 'DXA assessment >.*', 'DXA assessment')
```

Additional notes on RegEx. There are two types of RegEx.

The first class is this kind: `Path = str_replace(Path, '^Population characteristics >.*', 'Population characteristics')`.
Using this expression, for example, we keep [`Population characteristics`](http://biobank.ctsu.ox.ac.uk/crystal/browse.cgi?id=1&cd=category) as a category, but remove all of their descendants.

The second class is this kind: `Path = str_replace(Path, 'Touchscreen > ([^>]+) > .*', 'Touchscreen > \\1')`.
Using this expression, for example, we keep all the direct children of [`Touchscreen`](http://biobank.ctsu.ox.ac.uk/crystal/browse.cgi?id=100025&cd=category), but remove all of their descendants.

#### column descriptor for [`GBE_category.20200812.tsv`](GBE_category.20200812.tsv)

- GBE_category: the phenotype category listed in the PheWAS page on GBE
- GBE_ID: the GBE_ID of the phenotype
- N: the total number of case individuals (binary traits) or individuals with non-NA measurements (quantitative traits) within the full 500k set in UKB
- N_GBE: the `N` value in white British population
- N_NBW: the `N` value in non-British white population
- N_AFR: the `N` value in African population
- N_SAS: the `N` value in South Asian population
- N_EAS: the `N` value in East Asian population
- N_SMR: the `N` value in semi-related population
- N_OTH: the `N` value in the "others" population
- N_META: the total `N` value used in the UKB-wide meta-analysis
- GBE_NAME: the full phenotype name listed on GBE
- GBE_short_name: the short phenotype name listed on GBE
- GBE_short_name_len: the length of the short phenotype name
