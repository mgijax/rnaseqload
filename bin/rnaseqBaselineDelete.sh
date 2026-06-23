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

delete from GXD_HTSample_RNASeqCombined where _createdBy_key = 1673;
delete from GXD_HTSample_RNASeqSet where _createdBy_key = 1673;

EOSQL

