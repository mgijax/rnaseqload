#!/usr/bin/bash
#
#  run_downloadFiles.sh
###########################################################################
#
#  Purpose:
#
  Usage=run_downloadFiles.sh
#
#  Env Vars:
#
#      See the configuration file
#
#  Inputs:
#
#      - Common configuration file -
#               /usr/local/mgi/live/mgiconfig/master.config.sh
#      - Load configuration file - rnaseqload.config
#
#
#  Outputs:
#      - Array Express Files in INPUT_DIR
#      - Expression Atlas Files in the INPUT_DIR
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  Fatal error occurred
#      2:  Non-fatal error occurred
#
#  Assumes:  Nothing
#
# History:
#
# sc	06/25/2019 - created
#

cd `dirname $0`
LOG=`pwd`/run_downloadFiles.log
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

if [ -f ${LOG_DOWNLOAD} ]
    rm ${LOG_DOWNLOAD}
then
    touch ${LOG_DOWNLOAD}
fi

# check the lastrun file, exit if input file not updated
LASTRUN_FILE=${INPUTDIR}/lastrun
if [ -f ${LASTRUN_FILE} ]
then
    if test ${LASTRUN_FILE} -nt ${INPUT_FILE_DEFAULT}
    then

        echo "Input file has not been updated - skipping file download" | tee -a ${LOG_DOWNLOAD}
        exit 0
    fi
fi

# file updated; rm the download_ok file if it exists
echo "DOWNLOAD_OK: ${DOWNLOAD_OK}"
if [ -f ${DOWNLOAD_OK} ]
then
	rm ${DOWNLOAD_OK}
fi

#####################################
#
# Main
#
#####################################

date 

echo "Downloading input files" 
${RNASEQLOAD}/bin/downloadFiles.py #>> ${LOG_DOWNLOAD} 2>&1
STAT=$?

#
# Touch the "download_ok" file if all downloads successful
#
if [ ${STAT} = 0 ]
then
    touch ${DOWNLOAD_OK}
fi

