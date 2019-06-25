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

#
# createArchive, startLog, getConfigEnv
# sets "JOBKEY"

preload 

# remove files from output directory
cleanDir ${OUTPUTDIR}

# This cascades to GXD_HTSample_RNASeqSetMember
date | tee -a ${LOG_DIAG}
echo "Truncate GXD_HTSample_RNASeqSet table"  | tee -a ${LOG_DIAG}
${MGD_DBSCHEMADIR}/table/GXD_HTSample_RNASeqSet_truncate.object >> ${LOG_DIAG} 2>&1

date | tee -a ${LOG_DIAG}
echo "Loading Biological Replicates into Set tables" | tee -a ${LOG_DIAG}
${RNASEQLOAD}/bin/loadBioReps.py
STAT=$?
checkStatus ${STAT} "loadBioReps.py"

#date | tee -a ${LOG_DIAG}
#echo "Downloading input files" | tee -a ${LOG_DIAG}
#${RNASEQLOAD}/bin/downloadFiles.py >> ${LOG_DIAG}
#STAT=$?
#checkStatus ${STAT} "downloadFiles.py"

date | tee -a ${LOG_DIAG}
echo "Truncate GXD_HTSample_RNASeq table"  | tee -a ${LOG_DIAG}
${MGD_DBSCHEMADIR}/table/GXD_HTSample_RNASeq_truncate.object >> ${LOG_DIAG} 2>&1

date | tee -a ${LOG_DIAG}
echo "Truncate GXD_HTSample_RNASeqCombined table"  | tee -a ${LOG_DIAG}
${MGD_DBSCHEMADIR}/table/GXD_HTSample_RNASeqCombined_truncate.object >> ${LOG_DIAG} 2>&1

date | tee -a ${LOG_DIAG}
echo "Drop Indexes GXD_HTSample_RNASeq table"  | tee -a ${LOG_DIAG}
${MGD_DBSCHEMADIR}/index/GXD_HTSample_RNASeq_drop.object >> ${LOG_DIAG} 2>&1

date | tee -a ${LOG_DIAG}
echo "Drop Indexes GXD_HTSample_RNASeqCombined table"  | tee -a ${LOG_DIAG}
${MGD_DBSCHEMADIR}/index/GXD_HTSample_RNASeqCombined_drop.object >> ${LOG_DIAG} 2>&1

date | tee -a ${LOG_DIAG}
echo "Running rnaseqload.py" | tee -a ${LOG_DIAG}
${RNASEQLOAD}/bin/rnaseqload.py >> ${LOG_DIAG} 2>&1
STAT=$?
checkStatus ${STAT} "rnaseqload.py"

date | tee -a ${LOG_DIAG}
echo "Create Indexes GXD_HTSample_RNASeq table"  | tee -a ${LOG_DIAG}
${MGD_DBSCHEMADIR}/index/GXD_HTSample_RNASeq_create.object >> ${LOG_DIAG}  2>&1

date | tee -a ${LOG_DIAG}
echo "Create Indexes GXD_HTSample_RNASeqCombined table"  | tee -a ${LOG_DIAG}
${MGD_DBSCHEMADIR}/index/GXD_HTSample_RNASeqCombined_create.object >> ${LOG_DIAG}  2>&1

date | tee -a ${LOG_DIAG}
echo "Grant Database Permissions" | tee -a ${LOG_DIAG}
${PG_DBUTILS}/bin/grantPublicPerms.csh ${PG_DBSERVER} ${PG_DBNAME} >> ${LOG_DIAG} 2>&1

#
# run postload cleanup and email logs
#
shutDown
