#
#  Purpose:
#	refresh:
#		mgi_set
#		mgi_setmember
#		gxd_htsample_rnaseqset
#		gxd_htsample_rnaseqsetmember
#		gxd_htsample_rnaseqset_cache
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

echo "Run baseline setload" | tee -a $LOG
${SETLOAD}/setload.csh ${RNASEQLOAD}/baselinesetload.config | tee -a $LOG

echo "Reload RNASet/SetMember & RNASeqSet_Cache All MGI_Sets" | tee -a $LOG
${RNASEQLOAD}/bin/reloadBioRepsAndCache.sh | tee -a $LOG

cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 | tee -a $LOG

select s.* from MGI_Set s where s.name = 'Baseline RNASeq Load Experiment';

select a.accid, m.* 
into temp table baseline
from MGI_Set s, MGI_SetMember m , ACC_Accession a
where s.name = 'Baseline RNASeq Load Experiment'
and s._set_key = m._set_key
and s._mgitype_key = a._mgitype_key
and m._object_key = a._object_key
and a._logicaldb_key = 189
and a.preferred = 1
;

select * from baseline;

select b.*, s.*, m.*
from baseline b, GXD_HTSample_RNASeqSet s, GXD_HTSample_RNASeqSetMember m
where b._object_key = s._experiment_key
and s._rnaseqset_key = m._rnaseqset_key
;

--select * from GXD_HTSample_RNASeqSet_Cache limit 1;

EOSQL

