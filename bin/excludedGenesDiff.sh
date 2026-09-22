#!/bin/sh
#
# Purpose:
#	Generate list of tsvFile, tsvFileCount
#	This ${DIFFRAW_INPUTDIR}/tsvFileCount will be used by preprocessDiff.py
#	to find the list of genes that are to be excluded from processing.
#
#	Only genes that appear in all tsv files are to be processed.
#

cd `dirname $0`

CONFIG_LOAD=../rnaseqload.config

#
# Make sure the common configuration file exists and source it.
#
if [ -f ${CONFIG_LOAD} ]
then
    . ${CONFIG_LOAD}
else
    echo "Missing configuration file: ${CONFIG_LOAD}"
    exit 1
fi

cd ${DIFFRAW_INPUTDIR}
rm -rf tsvFile tsvCount tsvFileCount tsvGenesFile1 tsvGenesFile2 tsvGenesExcluded
touch tsvFile tsvCount tsvFileCount tsvGenesFile1 tsvGenesFile2 tsvGenesExcluded

# tsvFile = name of file
# tsvCount = count of genes
# tsvFileCount = tsvCount, tsfFile sorted by tsvCount
for i in *raw-counts.tsv
do
ls $i >> tsvFile
cut -f1 $i | sort | uniq | wc -l >> tsvCount
done
paste tsvCount tsvFile | sort > tsvFileCount

# ensembl ids from first row
# this is the good set of genes
head -n 1 tsvFileCount | cut -f2 | xargs cat | cut -f1 > tsvGenesFile1

# ensembl ids from last row which may contain differences
# this contains genes that need to be excluded from processing by preprocessDiff.py
# because these genes are only found in some of the tsv files
tail -n 1 tsvFileCount | cut -f2 | xargs cat | cut -f1 > tsvGenesFile2

# diff the first and last rows to find any genes that will need to be excluded
sort tsvGenesFile1 tsvGenesFile2 | uniq -u  > tsvGenesExcluded

# tsvGenesExcluded will be read by preprocessDiff.py
