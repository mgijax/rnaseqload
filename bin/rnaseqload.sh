#!/usr/bin/bash
#
#  rnaseqload.sh
###########################################################################
#
#  Purpose:
#
  Usage=rnaseqload.sh
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
#      - input file - see python script header
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
# sc	11/30/2015 - created
#

cd `dirname $0`
LOG=`pwd`/rnaseqload.log
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

#
#  Source the DLA library functions.
#

if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}" | tee -a ${LOG}
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined." | tee -a ${LOG}
    exit 1
fi

#####################################
#
# Main
#
#####################################
# remove logs (if not assocload logs will be appended)
#cleanDir ${LOGDIR}

#
# createArchive including OUTPUTDIR, startLog, getConfigEnv
# sets "JOBKEY"

preload ${OUTPUTDIR}

# remove files from output directory
cleanDir ${OUTPUTDIR}

echo "Downloading input files" >> ${LOG_DIAG}
${RNASEQLOAD}/bin/downloadFiles.py >> ${LOG_DIAG}
STAT=$?
checkStatus ${STAT} "downloadFiles.py"

echo "Running rnaseqload.py" >> ${LOG_DIAG}
${RNASEQLOAD}/bin/rnaseqload.py >> ${LOG_DIAG}
STAT=$?
checkStatus ${STAT} "rnaseqload.py"

#
# run postload cleanup and email logs
#
shutDown
