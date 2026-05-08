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

echo "Run baseline setload" 
${SETLOAD}/setload.csh ${RNASEQLOAD}/baselinesetload.config

cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 | tee -a $LOG
select s.* from MGI_Set s where s.name = 'Baseline RNASeq Load Experiment';
select a.accid, m.* from MGI_Set s, MGI_SetMember m , ACC_Accession a
where s.name = 'Baseline RNASeq Load Experiment'
and s._set_key = m._set_key
and s._mgitype_key = a._mgitype_key
and m._object_key = a._object_key
and a._logicaldb_key = 189
and a.preferred = 1
;
EOSQL


