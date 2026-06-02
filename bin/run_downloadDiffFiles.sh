#!/bin/sh
#
# Purpose:
#	Download Diff AES & EAE files
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

if [ -f ${DIFFLOG_DOWNLOAD} ]
    rm -rf ${DIFFLOG_DOWNLOAD}
then
    touch ${DIFFLOG_DOWNLOAD}
fi

LASTRUN_FILE=${DIFFRAW_INPUTDIR}/lastrun
if [ -f ${LASTRUN_FILE} ]
then
        echo "${LASTRUN_FILE} exists - skipping load" | tee -a ${DIFFLOG_DOWNLOAD}
        exit 0
fi

date | tee -a ${DIFFLOG_DOWNLOAD}

echo "Downloading input files" 
rm -rf ${DIFFLOG_DOWNLOAD}
rm -rf ${DIFFRAW_INPUTDIR}/*
${PYTHON} ${RNASEQLOAD}/bin/downloadDiffFiles.py >> ${DIFFLOG_DOWNLOAD} 2>&1
STAT=$?

#
# Touch the "lastrun" file to note when the load was run.
#
if [ ${STAT} = 0 ]
then
    touch ${LASTRUN_FILE}
fi

date | tee -a ${DIFFLOG_DOWNLOAD}
