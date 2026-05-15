#!/bin/csh -f

if ( ${?MGICONFIG} == 0 ) then
        setenv MGICONFIG /usr/local/mgi/live/mgiconfig
endif

source ${MGICONFIG}/master.config.csh

cd `dirname $0`

setenv LOG $0.log
rnamem -rf $LOG
touch $LOG
 
date | tee -a $LOG
 
# first time only; remove 1613 that are now in 1673
cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 >> ${BASELINELOG} 2>&1

-- remove old combo associated with markers that are now in baseline
select old._rnaseqcombined_key, oldset._experiment_key
into temp table toDelete
from GXD_HTSample_RNASeqCombined old, GXD_HTSample_RNASeqSet oldset
where old._CreatedBy_key = 1613
and old._rnaseqset_key = oldset._rnaseqset_key
and exists (select 1 from GXD_HTSample_RNASeqCombined new,  GXD_HTSample_RNASeqSet newset
	where oldset._experiment_key = newset._experiment_key
        and new._createdby_key = 1673
	and new._rnaseqset_key = newset._rnaseqset_key
        )
;

create index idx1 on toDelete(_rnaseqcombined_key);

delete from GXD_HTSample_RNASeqCombined
using toDelete
where toDelete._rnaseqcombined_key = GXD_HTSample_RNASeqCombined._rnaseqcombined_key
;

EOSQL

date |tee -a $LOG

