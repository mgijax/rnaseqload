#!/bin/sh
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

# remove all but Connie's file from the INPUT directory
#rm *.eae.txt
#rm *.aes.txt

echo "Downloading input files" >> ${LOG_DIAG}
${RNASEQLOAD}/bin/downloadFiles.py >> ${LOG_DIAG}
STAT=$?
checkStatus ${STAT} "downloadFiles.py"
#wget --no-check-certificate --tries=10 -nd -o /data/loads/sc/mgi/rnaseqload/logs/E-MTAB-7279.aes.log -O /data/loads/sc/mgi/rnaseqload/input/E-MTAB-7279.aes.txt https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7279/E-MTAB-7279.sdrf.txt
#wget --no-check-certificate --tries=10 -nd -o E-MTAB-7279.aes.log -O E-MTAB-7279.aes.txt https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7279/E-MTAB-7279.sdrf.txt
#
# run postload cleanup and email logs
#
shutDown
