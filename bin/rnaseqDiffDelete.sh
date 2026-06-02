#!/bin/sh
#
#  Purpose:
# 	delete the Differential RNASeqSet, RNASeqSetMember
#

if [ "${MGICONFIG}" = "" ]
then
        MGICONFIG=/usr/local/mgi/live/mgiconfig
	export MGICONFIG
	source ${MGICONFIG}/master.config.sh
fi

. ${RNASEQLOAD}/rnaseqload.config

rm -rf ${DIFFERENTIALLOG}
touch ${DIFFERENTIALLOG}

cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 >> ${DIFFERENTIALLOG} 2>&1

delete from GXD_HTSample_RNASeq where _CreatedBy_key = 1613;
delete from GXD_HTSample_RNASeqCombined where _CreatedBy_key = 1613;
delete from GXD_HTSample_RNASeqSet where _CreatedBy_key = 1613;

EOSQL

