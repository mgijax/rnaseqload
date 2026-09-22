#!/bin/sh
#
# Purpose:
#	Generate list of tsvFile, tsvFileCount
#	This ${DIFFRAW_INPUTDIR}/tsvFileCount will be used by preprocessDiff.py
#	to find the list of genes that are to be excluded from processing.
#
#	Only genes that appear in all tsv files are to be processed.
#
# example:
#	The tsv with the least amount of genes is the baseline
#
#	55574	E-MTAB-8840-raw-counts.tsv 
#	55574	E-MTAB-8964-raw-counts.tsv 
#	56749	E-MTAB-9000-raw-counts.tsv 
#
#	E-MTAB-8840, E-MTAB-8964 has the lowest tsvCount, so that's the baseline count
#	E-MTAB-9000 count is *not* the baseline count
#
#	preprocessDiff.py/loadExcludedGenes
#		compare baseline ensembl ids with non-baseline ensembl ids
#		add the diff to the exclucedGenes list
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
rm -rf tsvFile tsvCount tsvFileCount
touch tsvFile tsvCount tsvFileCount
for i in *raw-counts.tsv
do
ls $i >> tsvFile
cut -f1 $i | sort | uniq | wc -l >> tsvCount
done
paste tsvCount tsvFile | sort > tsvFileCount

