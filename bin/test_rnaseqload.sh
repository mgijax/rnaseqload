#!/usr/bin/bash
#
#  rnaseqload.sh
###########################################################################
#
#  Purpose: Wrapper that determines if the load needs to run, then runs
#           the load. This script is responsible for truncating the RNA Seq
#	    tables and dropping/recreating indexes
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
# sc	05/23/2019 - created
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

#
# createArchive, startLog, getConfigEnv
# sets "JOBKEY"

preload 

# remove files from output directory and the input directory where 
# the intermediate files are created.
cleanDir ${OUTPUTDIR} ${INPUTDIR}

date | tee -a ${LOG_DIAG}
echo "Running test_rnaseqload.py" | tee -a ${LOG_DIAG}
${RNASEQLOAD}/bin/test_rnaseqload.py >> ${LOG_DIAG} 2>&1
STAT=$?
checkStatus ${STAT} "test_rnaseqload.py"

shutDown
