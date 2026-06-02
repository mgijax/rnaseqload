#!/bin/sh
#
# Purpose: Wrapper for Baseline RNASeq load
#
# Step 1: process withdrawn markers
# Step 2: delete existing Baseline RNASeqSet data
# Step 3: run baseline MGI_Set, MGI_SetMember
# Step 4: run baseline download (raw_input_baseline)
# Step 5: run baseline pre processing (input_baseline)
# Step 6: run baseline: RNASeqSet, RNASeq_SetMember, RNASeqCombined
#
# MGI_Set = Baseline RNASeq Experiments
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

rm -rf ${BASELINELOG}
touch ${BASELINELOG}

# make sure this is uncommented before tagging
#date >> ${BASELINELOG} 2>&1
#echo "Step 1: process withdrawn markers" >> ${BASELINELOG} 2>&1
#${RNASEQLOAD}/bin/processWithdrawnMarkers.sh >> ${BASELINELOG} 2>&1

date >> ${BASELINELOG} 2>&1
echo "Step 2: delete existing Baseline RNASeqSet data" >> ${BASELINELOG} 2>&1
${RNASEQLOAD}/bin/rnaseqBaselineDelete.sh >> ${BASELINELOG} 2>&1

date >> ${BASELINELOG} 2>&1
echo "Step 3: run baseline MGI_Set, MGI_SetMember" >> ${BASELINELOG} 2>&1
${RNASEQLOAD}/bin/run_setBaseline.sh >> ${BASELINELOG} 2>&1

date >> ${BASELINELOG} 2>&1
echo "Step 4: run baseline download (raw_input_baseline)" >> ${BASELINELOG} 2>&1
${RNASEQLOAD}/bin/run_downloadBaselineFiles.sh >> ${BASELINELOG} 2>&1

date >> ${BASELINELOG} 2>&1
echo "Step 5: run baseline pre processing (input_baseline)" >> ${BASELINELOG} 2>&1
rm -rf ${BASELINEINPUTDIR}/*
${PYTHON} ${RNASEQLOAD}/bin/preprocessBaseline.py >> ${BASELINELOG} 2>&1

date >> ${BASELINELOG} 2>&1
echo "Step 6: run baseline: RNASeqSet, RNASeq_SetMember, RNASeqCombined" >> ${BASELINELOG} 2>&1
${MGD_DBSCHEMADIR}/index/GXD_HTSample_RNASeqCombined_drop.object >> ${BASELINELOG} 2>&1
${MGD_DBSCHEMADIR}/key/GXD_HTSample_RNASeqCombined_drop.object >> ${BASELINELOG} 2>&1
${PYTHON} ${RNASEQLOAD}/bin/rnaseqBaseline.py >> ${BASELINELOG} 2>&1
${MGD_DBSCHEMADIR}/key/GXD_HTSample_RNASeqCombined_create.object >> ${BASELINELOG} 2>&1
${MGD_DBSCHEMADIR}/index/GXD_HTSample_RNASeqCombined_create.object >> ${BASELINELOG} 2>&1

date >> ${BASELINELOG} 2>&1
