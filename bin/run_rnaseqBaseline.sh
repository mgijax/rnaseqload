#!/bin/sh
#
#  Purpose:
#	refresh:
#		mgi_set = Baseline RNASeq Experiment
#		mgi_setmember
#		gxd_htsample_rnaseqset
#		gxd_htsample_rnaseqsetmember
#		gxd_htsample_rnaseqset_cache
#

if [ "${MGICONFIG}" = "" ]
then
        MGICONFIG=/usr/local/mgi/live/mgiconfig
	export MGICONFIG
	source ${MGICONFIG}/master.config.sh
fi

. ${RNASEQLOAD}/rnaseqload.config

rm -rf ${BASELINELOG}
touch ${BASELINELOG}

#
# delete the Baseline RNASeqSet, RNASeqSetMember
#
cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 >> ${BASELINELOG} 2>&1

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

--delete from GXD_HTSample_RNASeqSet
--using baseline
--where baseline._object_key = GXD_HTSample_RNASeqSet._experiment_key
--;

select b.*, s.*, m.*
from baseline b, GXD_HTSample_RNASeqSet s, GXD_HTSample_RNASeqSetMember m
where b._object_key = s._experiment_key
and s._rnaseqset_key = m._rnaseqset_key
;

EOSQL

#date >> ${BASELINELOG} 2>&1
#echo "Step 1: run baseline download (raw_input_baseline)" >> ${BASELINELOG} 2>&1
#${RNASEQLOAD}/bin/run_downloadBaselineFiles.sh >> ${BASELINELOG} 2>&1

#date >> ${BASELINELOG} 2>&1
#echo "Step 2: run baseline MGI_Set, MGI_SetMember" >> ${BASELINELOG} 2>&1
#${RNASEQLOAD}/bin/run_setbaseline.sh >> ${BASELINELOG} 2>&1

date >> ${BASELINELOG} 2>&1
echo "Step 3: run baseline pre processing (input_baseline)" 
${PYTHON} ${RNASEQLOAD}/bin/preprocessBaseline.py >> ${BASELINELOG} 2>&1

#date >> ${BASELINELOG} 2>&1
#echo "Step 4: process withdrawn markers"
#${RNASEQLOAD}/bin/processWithdrawnMarkers.sh >> ${BASELINELOG} 2>&1

#date >> ${BASELINELOG} 2>&1
#echo "Step 4: run baseline: RNASeqSet, RNASeq_SetMember, RNASeq, RNASeqCombined" >> ${BASELINELOG} 2>&1
#${PYTHON} ${RNASEQLOAD}/bin/rnaseqBaseline.py >> ${BASELINELOG} 2>&1

#date >> ${BASELINELOG} 2>&1
#echo "Step 5: run loadSeqSetCache.sh" >> ${BASELINELOG} 2>&1
#${RNASEQLOAD}/bin/loadSeqSetCache.sh >> ${BASELINELOG} 2>&1

date >> ${BASELINELOG} 2>&1
