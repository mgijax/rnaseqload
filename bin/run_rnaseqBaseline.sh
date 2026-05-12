#!/bin/sh
#
# Purpose:
#	Run Baselinne Processing
#

cd `dirname $0`
LOG=`pwd`/rnaseqBaseline.log
rm -rf ${LOG}

CONFIG_LOAD=../rnaseqload.config

#
# verify & source the configuration file
#

if [ ! -r ${CONFIG_LOAD} ]
then
    echo "Cannot read configuration file: ${CONFIG_LOAD}"
    exit 1
fi

. ${CONFIG_LOAD}

rm -rf ${BASELINELOG}

#echo "Step 1: Baseline MGI_Set only, all RNASeqSet, RNASeqSetMember, RNASeqSet_Cache" 
#${RNASEQLOAD}/bin/run_setbaseline.sh >> ${BASELINELOG} 2>&1

#echo "Step 2: Baseline Pre Processing" 
#${PYTHON} ${RNASEQLOAD}/bin/preprocessBaseline.py >> ${BASELINELOG} 2>&1

#echo "Step 3: Process Withdrawn Markers"
#${RNASEQLOAD}/bin/processWithdrawnMarkers.sh >> ${BASELINELOG} 2>&1

echo "Step 4: Baseline Processing : RNASeq, RNASeqCombined" 
${PYTHON} ${RNASEQLOAD}/bin/rnaseqBaseline.py >> ${BASELINELOG} 2>&1

date | tee -a ${BASELINELOG}
