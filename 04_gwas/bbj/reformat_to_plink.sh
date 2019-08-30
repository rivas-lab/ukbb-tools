#!/bin/bash

# Usage: sbatch -p mrivas -t 12:00:00 --mem=150000 reformat_to_plink.sh 

#QT

##CHROM POS ID REF ALT A1 TEST OBS_CT BETA SE T_STAT P

#biomarkers (majority)
#SNP	CHR	POS	REF	ALT	Frq	Rsq	BETA	SE	P	LOG10P	N 
for file in $(ls /oak/stanford/groups/mrivas/bbj/raw_sumstats/qt); do
    echo $file
    output=$(dirname $(dirname $(dirname $file)))/plink_format/qt/$(basename ${file%.*})
    zcat /oak/stanford/groups/mrivas/bbj/raw_sumstats/qt/$file | awk -F'\t' 'BEGIN{OFS='\t'} {if (NR > 1) {print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$5"\tADD\t"$6"\t"$8"\t"$9"\t"$8/$9"\t"$10} else { print "#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP"}}' >$output
    gzip -f $output
done

#smoking QT
#SNP	CHR	POS	A1	A2	A1Frq	Rsq	BETA	SE	P
for file in $(ls /oak/stanford/groups/mrivas/bbj/raw_sumstats/qt/*smok*); do
    echo $file
    output=$(dirname $(dirname $(dirname $file)))/plink_format/qt/$(basename ${file%.*})
    zcat $file | awk -F'\t' 'BEGIN{OFS='\t'} {if (NR > 1) {print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$4"\tADD\t"1-$6"\t"$8"\t"$9"\t"$8/$9"\t"$10} else { print "#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP"}}' >$output
    gzip -f $output
done

for file in $(ls /oak/stanford/groups/mrivas/bbj/raw_sumstats/qt/*cig*); do
    echo $file
    output=$(dirname $(dirname $(dirname $file)))/plink_format/qt/$(basename ${file%.*})
    zcat $file | awk -F'\t' 'BEGIN{OFS='\t'} {if (NR > 1) {print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$4"\tADD\t"1-$6"\t"$8"\t"$9"\t"$8/$9"\t"$10} else { print "#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP"}}' >$output
    gzip -f $output
done

#BIN

#CHROM POS ID REF ALT A1 FIRTH? TEST OBS_CT OR SE Z_STAT P

#smoking BIN
#SNP	CHR	POS	A1	A2	A1Frq	Rsq	BETA	SE	P
for file in $(ls /oak/stanford/groups/mrivas/bbj/raw_sumstats/bin/*smok*); do
    echo $file
    output=$(dirname $(dirname $(dirname $file)))/plink_format/bin/$(basename ${file%.*})
    zcat $file | awk -F'\t' 'BEGIN{OFS='\t'} {if (NR > 1) {print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$4"\tN\tADD\t"1-$6"\t"exp($8)"\t"$9"\t"$8/$9"\t"$10} else { print "#CHROM\tPOS\tID\tREF\tALT\tA1\tFIRTH?\tTEST\tOBS_CT\tOR\tSE\tZ_STAT\tP"}}' >$output
    gzip -f $output
done

#AF BIN
#SNP	CHR	POS	A1	A2	FREQ1	CASE_FREQ1	CTRL_FREQ1	RSQR	EFFECT1	OR	STDERR	WALDCHISQ	PVALUE
for file in $(ls /oak/stanford/groups/mrivas/bbj/raw_sumstats/bin/*AF*); do
    echo $file
    output=$(dirname $(dirname $(dirname $file)))/plink_format/bin/$(basename ${file%.*})
    zcat $file | awk -F'\t' 'BEGIN{OFS='\t'} {if (NR > 1) {print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$4"\tN\tADD\t"1-$6"\t"$11"\t"$12"\t"$10/$12"\t"$14} else { print "#CHROM\tPOS\tID\tREF\tALT\tA1\tFIRTH?\tTEST\tOBS_CT\tOR\tSE\tZ_STAT\tP"}}' >$output
    gzip -f $output
done

#POAG BIN
#SNP	CHR	POS	REF	ALT	Case_frq	Ctrl_frq	Rsq	BETA	OR	SE	P
for file in $(ls /oak/stanford/groups/mrivas/bbj/raw_sumstats/bin/*POAG*); do
    echo $file
    output=$(dirname $(dirname $(dirname $file)))/plink_format/bin/$(basename ${file%.*})
    zcat $file | awk -F'\t' 'BEGIN{OFS='\t'} {if (NR > 1) {print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$5"\tN\tADD\t"$7"\t"$10"\t"$11"\t"$9/$11"\t"$12} else { print "#CHROM\tPOS\tID\tREF\tALT\tA1\tFIRTH?\tTEST\tOBS_CT\tOR\tSE\tZ_STAT\tP"}}' >$output
    gzip -f $output
done

#RA BIN
#SNP	CHR	POS	A1	A2	CASE_FREQ1	CTRL_FREQ1	BETA	SE	PVALUE	N_CASE	N_CTRL
for file in $(ls /oak/stanford/groups/mrivas/bbj/raw_sumstats/bin/*RA*); do
    echo $file
    output=$(dirname $(dirname $(dirname $file)))/plink_format/bin/$(basename ${file%.*})
    zcat $file | awk -F'\t' 'BEGIN{OFS='\t'} {if (NR > 1) {print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$4"\tN\tADD\t"$7"\t"exp($8)"\t"$9"\t"$8/$9"\t"$10} else { print "#CHROM\tPOS\tID\tREF\tALT\tA1\tFIRTH?\tTEST\tOBS_CT\tOR\tSE\tZ_STAT\tP"}}' >$output
    gzip -f $output
done

#colorectal BIN
##MARKER	CHR	POS	EffectAllele	NonEffectAllele	FREQ1	CASE_FREQ1	CTRL_FREQ1	RSQR	EFFECT1	OR	STDERR	WALDCHISQ	PVALUE	LRCHISQ	LRPVAL	contig	contigloc	num_of_gene	gene	relativeloc
for file in $(ls /oak/stanford/groups/mrivas/bbj/raw_sumstats/bin/*colorectal*); do
    echo $file
    output=$(dirname $(dirname $(dirname $file)))/plink_format/bin/$(basename ${file%.*})
    zcat $file | awk -F'\t' 'BEGIN{OFS='\t'} {if (NR > 1 && $13 != "NA") {print $2"\t"$3"\t"$1"\t"$5"\t"$4"\t"$4"\tN\tADD\t"1-$6"\t"$11"\t"$12"\t"$10/$12"\t"$14} else if (NR > 1 && $13 == "NA") {print $2"\t"$3"\t"$17"\t"$5"\t"$4"\t"$4"\tY\tADD\t"1-$6"\t"$11"\t"$12"\t"$10/$12"\t"$16} else { print "#CHROM\tPOS\tID\tREF\tALT\tA1\tFIRTH?\tTEST\tOBS_CT\tOR\tSE\tZ_STAT\tP"}}' >$output
    gzip -f $output
done

#T2D
for file in $(ls /oak/stanford/groups/mrivas/bbj/raw_sumstats/bin/*T2D*); do
    echo $file
    output=$(dirname $(dirname $(dirname $file)))/plink_format/bin/$(basename ${file%.*})
    zcat $file | awk -F'\t' 'BEGIN{OFS='\t'} {if (NR > 1) {print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$5"\tN\tADD\t"1-$6"\t"exp($7)"\t"$8"\t"$7/$8"\t"$9} else { print "#CHROM\tPOS\tID\tREF\tALT\tA1\tFIRTH?\tTEST\tOBS_CT\tOR\tSE\tZ_STAT\tP"}}' >$output
    gzip -f $output
done
