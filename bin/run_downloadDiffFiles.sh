#!/bin/sh
#
# Purpose:
#	Download Raw Differential files
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

rm -rf ${DIFFLOG_DOWNLOAD}

#LASTRUN_FILE=${DIFFRAW_INPUTDIR}/lastrun
#if [ -f ${LASTRUN_FILE} ]
#then
#        echo "${LASTRUN_FILE} exists - skipping run_downloadDiffFiles.sh" | tee -a ${DIFFLOG_DOWNLOAD}
#        exit 0
#fi

date | tee -a ${DIFFLOG_DOWNLOAD}

echo "Downloading input files" 
rm -rf ${DIFFRAW_INPUTDIR}/*
${PYTHON} ${RNASEQLOAD}/bin/downloadDiffFiles.py >> ${DIFFLOG_DOWNLOAD} 2>&1

touch ${LASTRUN_FILE}

date | tee -a ${DIFFLOG_DOWNLOAD}
