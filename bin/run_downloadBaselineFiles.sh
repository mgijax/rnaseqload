#!/bin/sh
#
# Purpose:
#	Download Raw Baseline files
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

rm -rf ${BASELINELOG_DOWNLOAD}

LASTRUN_FILE=${BASELINERAW_INPUTDIR}/lastrun
if [ -f ${LASTRUN_FILE} ]
then
        echo "${LASTRUN_FILE} exists - skipping run_downloadBaselineFiles.sh" | tee -a ${BASELINELOG_DOWNLOAD}
        exit 0
fi

date | tee -a ${BASELINELOG_DOWNLOAD}

echo "Downloading input files" 
rm -rf ${BASELINERAW_INPUTDIR}/*
${PYTHON} ${RNASEQLOAD}/bin/downloadBaselineFiles.py >> ${BASELINELOG_DOWNLOAD} 2>&1

touch ${LASTRUN_FILE}

date | tee -a ${BASELINELOG_DOWNLOAD}
