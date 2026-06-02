#!/bin/sh
#
# Purpose:
#	Download Baseline AES & EAE files
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

if [ -f ${BASELINELOG_DOWNLOAD} ]
    rm -rf ${BASELINELOG_DOWNLOAD}
then
    touch ${BASELINELOG_DOWNLOAD}
fi

LASTRUN_FILE=${BASELINERAW_INPUTDIR}/lastrun
if [ -f ${LASTRUN_FILE} ]
then
        echo "${LASTRUN_FILE} exists - skipping load" | tee -a ${BASELINELOG_DOWNLOAD}
        exit 0
fi

date | tee -a ${BASELINELOG_DOWNLOAD}

echo "Downloading input files" 
rm -rf ${BASELINELOG_DOWNLOAD}
rm -rf ${BASELINERAW_INPUTDIR}/*
${PYTHON} ${RNASEQLOAD}/bin/downloadBaselineFiles.py >> ${BASELINELOG_DOWNLOAD} 2>&1
STAT=$?

#
# Touch the "lastrun" file to note when the load was run.
#
if [ ${STAT} = 0 ]
then
    touch ${LASTRUN_FILE}
fi

date | tee -a ${BASELINELOG_DOWNLOAD}
