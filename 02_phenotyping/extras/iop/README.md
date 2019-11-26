## IOP (Intraocular pressure) phenotype definition

We will define the median (of left and right) for using the following fields


```
#GBE_ID GBE_NAME
INI5254 Intra-ocular_pressure,_corneal-compensated_(right)
INI5255 Intra-ocular_pressure,_Goldmann-correlated_(right)
INI5262 Intra-ocular_pressure,_corneal-compensated_(left)
INI5263 Intra-ocular_pressure,_Goldmann-correlated_(left)
```

To use the most recent data, we check in which table we have those fields.

```
$ ukbb-query_find_table_by_field_id.sh 5254
http://bit.ly/UKB24983_tables
ukb37855        4041    f.5254.0.0      f       5254    0       0
ukb25826        2070    f.5254.0.0      f       5254    0       0
ukb21731        2070    f.5254.0.0      f       5254    0       0
ukb10137        2070    f.5254.0.0      f       5254    0       0
```

It seems like we have them in `2005693/37855`.

```
$ cat /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/2005693/37855/download/ukb37855.tab.columns | egrep 'f.5254|f.5255|f.5262|f.5263'
f.5254.0.0      f       5254    0       0
f.5254.1.0      f       5254    1       0
f.5255.0.0      f       5255    0       0
f.5255.1.0      f       5255    1       0
f.5262.0.0      f       5262    0       0
f.5262.1.0      f       5262    1       0
f.5263.0.0      f       5263    0       0
f.5263.1.0      f       5263    1       0
```

We defined the followings

- `INI2005254`: Intra-ocular_pressure,_corneal-compensated_(median)
- `INI2005255`: Intra-ocular_pressure,_Goldmann-correlated_(median)
