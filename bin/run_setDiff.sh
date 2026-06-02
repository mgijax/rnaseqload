#!/bin/csh
#
#  Purpose:
#	quick script to run the baseline setload
#

if ( ${?MGICONFIG} == 0 ) then
        setenv MGICONFIG /usr/local/mgi/live/mgiconfig
endif

source ${MGICONFIG}/master.config.csh
source ${RNASEQLOAD}/baselinesetload.config

cd ${LOADDIR}

setenv LOG $0.log
rm -rf $LOG
touch $LOG

setenv LASTRUN_FILE ${SETDATADIR}/lastrun.diff
if ( -e ${LASTRUN_FILE} ) then
        echo "LASTRUN_FILE exists - skipping load" | tee -a ${LOG}
        exit 0
endif

echo "Run baseline setload"  | tee -a ${LOG}
${SETLOAD}/setload.csh ${RNASEQLOAD}/baselinesetload.config | tee -a ${LOG}

cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 | tee -a $LOG
select s.* from MGI_Set s where s.name = 'Baseline RNASeq Load Experiments';
select a.accid, m.* from MGI_Set s, MGI_SetMember m , ACC_Accession a
where s.name = 'Baseline RNASeq Load Experiments'
and s._set_key = m._set_key
and s._mgitype_key = a._mgitype_key
and m._object_key = a._object_key
and a._logicaldb_key = 189
and a.preferred = 1
;
EOSQL

#
# Touch the "lastrun" file to note when the load was run.
#
touch ${LASTRUN_FILE}
