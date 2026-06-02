#!/bin/sh
#
# Purpose: Wrapper for Diff RNASeq load
#
# Step 1: process withdrawn markers
#
# Step 2: delete existing Diff RNASeqSet data
#
# only run if new Diff member is added
# Step 3: run differential MGI_Set, MGI_SetMember
#
# only run if Step 3 is run OR if need to refresh raw input files
# Step 4: run differential download (raw_input_differential)
# Step 5: run differential pre processing (input_differential)
#
# Step 6: run differential: RNASeqSet, RNASeq_SetMember, RNASeqCombined
#
# MGI_Set = Diff RNASeq Experiments
# MGI_SetMember
# GXD_HTSample_RNASeqSet
# GXD_HTSample_RNASeqSetMember
# GXD_HTSample_RNASeqSetCombined
#

if [ "${MGICONFIG}" = "" ]
then
        MGICONFIG=/usr/local/mgi/live/mgiconfig
	export MGICONFIG
	source ${MGICONFIG}/master.config.sh
fi

. ${RNASEQLOAD}/rnaseqload.config

rm -rf ${DIFFLOG}
touch ${DIFFLOG}

# make sure this is uncommented before tagging
#date >> ${DIFFLOG} 2>&1
#echo "Step 1: process withdrawn markers" >> ${DIFFLOG} 2>&1
#${RNASEQLOAD}/bin/processWithdrawnMarkers.sh >> ${DIFFLOG} 2>&1

date >> ${DIFFLOG} 2>&1
echo "Step 2: delete existing Diff RNASeqSet data" >> ${DIFFLOG} 2>&1
${RNASEQLOAD}/bin/rnaseqDiffDelete.sh >> ${DIFFLOG} 2>&1

# only run if adding a new Diff member
date >> ${DIFFLOG} 2>&1
echo "Step 3: run differential MGI_Set, MGI_SetMember" >> ${DIFFLOG} 2>&1
${RNASEQLOAD}/bin/run_setDiff.sh >> ${DIFFLOG} 2>&1

# only needs to run if Step 4 is run
date >> ${DIFFLOG} 2>&1
echo "Step 4: run differential download (raw_input_differential)" >> ${DIFFLOG} 2>&1
${RNASEQLOAD}/bin/run_downloadDiffFiles.sh >> ${DIFFLOG} 2>&1

# only needs to run if Step 5 is run
date >> ${DIFFLOG} 2>&1
echo "Step 5: run differential pre processing (input_differential)" >> ${DIFFLOG} 2>&1
rm -rf ${DIFFINPUTDIR}/*
${PYTHON} ${RNASEQLOAD}/bin/preprocessDiff.py >> ${DIFFLOG} 2>&1

#date >> ${DIFFLOG} 2>&1
#echo "Step 6: run differential: RNASeqSet, RNASeq_SetMember, RNASeqCombined" >> ${DIFFLOG} 2>&1
#${MGD_DBSCHEMADIR}/index/GXD_HTSample_RNASeqCombined_drop.object >> ${DIFFLOG} 2>&1
#${MGD_DBSCHEMADIR}/key/GXD_HTSample_RNASeqCombined_drop.object >> ${DIFFLOG} 2>&1
#${PYTHON} ${RNASEQLOAD}/bin/rnaseqDiff.py >> ${DIFFLOG} 2>&1
#${MGD_DBSCHEMADIR}/key/GXD_HTSample_RNASeqCombined_create.object >> ${DIFFLOG} 2>&1
#${MGD_DBSCHEMADIR}/index/GXD_HTSample_RNASeqCombined_create.object >> ${DIFFLOG} 2>&1

date >> ${DIFFLOG} 2>&1
