#!/bin/csh
#
#  Purpose:
#	refresh:
#		mgi_set = Baseline RNASeq Experiment
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

rm -rf ${BASELINELOG}
touch ${BASELINELOG}

#
# delete the Baseline RNASeqSet, RNASeqSetMember
#
cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 | tee -a ${BASELINELOG}

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

delete from GXD_HTSample_RNASeqSet
using baseline
where baseline._object_key = GXD_HTSample_RNASeqSet._experiment_key
;

select b.*, s.*, m.*
from baseline b, GXD_HTSample_RNASeqSet s, GXD_HTSample_RNASeqSetMember m
where b._object_key = s._experiment_key
and s._rnaseqset_key = m._rnaseqset_key
;

EOSQL

#date | tee -a ${BASELINELOG}
#echo "Run baseline setload" | tee -a ${BASELINELOG}
#cd ${LOADDIR}
#setenv SETLOG_DIAG $0.log
#rm -rf ${SETLOG_DIAG}
#touch ${SETLOG_DIAG}
#${SETLOAD}/setload.csh ${RNASEQLOAD}/baselinesetload.config | tee -a ${SETLOG_DIAG}

date | tee -a ${BASELINELOG}
echo "Run baseline: RNASeqSet, RNASeq_SetMember, RNASeq, RNASeqCombined" | tee -a ${BASELINELOG}
$PYTHON ${RNASEQLOAD}/bin/rnaseqBaseline.py | tee -a ${BASELINELOG}

#date | tee -a ${BASELINELOG}
#echo "Run loadSeqSetCache.sh" | tee -a ${BASELINELOG}
#${RNASEQLOAD}/bin/loadSeqSetCache.sh | tee -a ${BASELINELOG}

