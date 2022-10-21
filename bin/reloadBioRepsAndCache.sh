#!/usr/bin/bash
#
# reloadBioRepsAndCache.sh
###########################################################################
#
#  Purpose: Wrapper that reloads the biological replicates (set and set members)
#           and the gxd_htsample_rnaseqset_cache
#
  Usage=reloadBioRepsAndCache.sh
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
# sc	10/21/2022 - created
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

# override LOG_DIAG, and MAIL_LOADNAME
LOG_DIAG=${LOGDIR}/reloadBioRepsAndCache.diag.log
MAIL_LOADNAME="GXD RNA-Seq Load Reload Bio Replicates and Cache"
export LOG_DIAG MAIL_LOADNAME

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

# This cascades to GXD_HTSample_RNASeqSetMember
date | tee -a ${LOG_DIAG}
echo "Truncate GXD_HTSample_RNASeqSet table"  | tee -a ${LOG_DIAG}
${MGD_DBSCHEMADIR}/table/GXD_HTSample_RNASeqSet_truncate.object >> ${LOG_DIAG} 2>&1

date | tee -a ${LOG_DIAG}
echo "Loading Biological Replicates into Set tables" | tee -a ${LOG_DIAG}
${PYTHON} ${RNASEQLOAD}/bin/loadBioReps.py
STAT=$?
checkStatus ${STAT} "loadBioReps.py"

date | tee -a ${LOG_DIAG}
echo "Running loadSeqSetCache.sh" | tee -a ${LOG_DIAG}
${RNASEQLOAD}/bin/loadSeqSetCache.sh >> ${LOG_DIAG} 2>&1
STAT=$?
checkStatus ${STAT} "loadSeqSetCache.sh"

echo "Grant Database Permissions" | tee -a ${LOG_DIAG}
${PG_DBUTILS}/bin/grantPublicPerms.csh ${PG_DBSERVER} ${PG_DBNAME} mgd >> ${LOG_DIAG} 2>&1

shutDown
