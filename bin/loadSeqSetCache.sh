#!/usr/bin/bash
#
# loadSeqSetCache.sh 
###########################################################################
#
#  Purpose: Wrapper that runs the loadSeqSetCache python script
#           and truncates the cache table
#
#  Usage=loadSeqSetCache.sh
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
#  Outputs:
#
#  Exit Codes:
#
#  Assumes:  Nothing
#
# History:
#
# sc	09/09/2022 - created
#

cd `dirname $0`
LOG=`pwd`/loadSeqSetCache.log
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

#####################################
#
# Main
#
#####################################

date >> ${LOG_DIAG}
echo "Create temp tables for the input data" >> ${LOG_DIAG}
cat - <<EOSQL | psql -h${MGD_DBSERVER} -d${MGD_DBNAME} -U mgd_dbo -e >>  ${LOG_DIAG}

create table rnaseqsetcombined (
    _rnaseqcombined_key int not null,
    _rnaseqset_key      int not null
);

create index idx1 on rnaseqsetcombined(_rnaseqset_key);
create index idx2 on rnaseqsetcombined (_rnaseqcombined_key);

grant all on rnaseqsetcombined to public;
grant all on rnaseqsetcombined to mgd_dbo;

insert into rnaseqsetcombined
select distinct rs._rnaseqcombined_key,
            sm._rnaseqset_key
        from gxd_htsample_rnaseq rs, gxd_htsample_rnaseqsetmember sm
        where rs._sample_key = sm._sample_key;

EOSQL


date | tee -a ${LOG_DIAG}
echo "Truncate GXD_HTSample_RNASeqSet_Cache table"  | tee -a ${LOG_DIAG}
${MGD_DBSCHEMADIR}/table/GXD_HTSample_RNASeqSet_Cache_truncate.object >> ${LOG_DIAG} 2>&1

date | tee -a ${LOG_DIAG}
echo "Loading Seq Set Cache table" | tee -a ${LOG_DIAG}
${PYTHON} ${RNASEQLOAD}/bin/loadSeqSetCache.py
STAT=$?
echo "STAT: ${STAT}"

date | tee -a ${LOG_DIAG}
echo "Grant Database Permissions" | tee -a ${LOG_DIAG}
${PG_DBUTILS}/bin/grantPublicPerms.csh ${PG_DBSERVER} ${PG_DBNAME} mgd >> ${LOG_DIAG} 2>&1

#
# Drop the temp table
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo "Drop the temp table" >> ${LOG_DIAG}
cat - <<EOSQL | psql -h${MGD_DBSERVER} -d${MGD_DBNAME} -U mgd_dbo -e  >> ${LOG_DIAG}

drop table rnaseqsetcombined;

EOSQL

