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

#####################################
#
# Main
#
#####################################

date  >>  ${LOG_DOWNLOAD}

echo "Downloading input files" 
${RNASEQLOAD}/bin/downloadFiles.py >> ${LOG_DOWNLOAD} 2>&1
STAT=$?

date  >>  ${LOG_DOWNLOAD}
