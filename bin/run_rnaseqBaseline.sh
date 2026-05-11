#!/bin/sh
#
# Purpose:
#	Run Baselinne Processing
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

echo "Baseline Processing" 
rm -rf ${BASELINELOG}
${PYTHON} ${RNASEQLOAD}/bin/rnaseqBaseline.py >> ${BASELINELOG} 2>&1
STAT=$?

date | tee -a ${BASELINELOG_DOWNLOAD}
