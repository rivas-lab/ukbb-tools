#!/bin/bash

cd ../tables/Data_Dictionary_Showcases/
wget http://biobank.ctsu.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.csv && mv Data_Dictionary_Showcase.csv Data_Dictionary_Showcase.$(date '+%Y%m%d').csv
cd ..
ln -sf Data_Dictionary_Showcases/Data_Dictionary_Showcase.$(date '+%Y%m%d').csv Data_Dictionary_Showcase.csv
