#!/bin/sh
#
#  Purpose:
# 	delete the Baseline RNASeqSet, RNASeqSetMember
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

cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 >> ${BASELINELOG} 2>&1

-- start: delete _CreatedBy_key = 1613 for the baseline experiments
-- this only needs to be run once & can then be removed

select s.*
into temp toDelete
from GXD_HTSample_RNASeqSet s, MGI_SetMember sm
where s._createdBy_key = 1613
and s._experiment_key = sm._object_key
and sm._set_key = 1061
;

create index idx1 on toDelete (_rnaseqset_key);

delete from GXD_HTSample_RNASeqCombined s
using toDelete d
where d._rnaseqset_key = s._rnaseqset_key
;

delete from GXD_HTSample_RNASeqSet s
using toDelete d
where d._rnaseqset_key = s._rnaseqset_key
;

-- end: delete _CreatedBy_key = 1613 for the baseline experiments

delete from GXD_HTSample_RNASeqCombined where _createdBy_key = 1673;
delete from GXD_HTSample_RNASeqSet where _createdBy_key = 1673;

EOSQL

