#!/bin/csh -f

#
# Process Withdrawn Markers
#
# if there is no “new” marker, then delete the GXD_HTSample_RNASeq & Combined rows
# else, change “old” marker key to “new marker key
#

if ( ${?MGICONFIG} == 0 ) then
        setenv MGICONFIG /usr/local/mgi/live/mgiconfig
endif

source ${MGICONFIG}/master.config.csh

cd `dirname $0`

setenv LOG ${DATALOADSOUTPUT}/mgi/rnaseqload/logs/processWithdrawnMarkers.sh.log
rm -rf $LOG
touch $LOG

date |tee -a $LOG
 
cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 |tee -a $LOG

-- find rnaseq rows that contains withdrawn markers (_marker_status_key = 2)
select distinct m._marker_key, m._marker_status_key, m.symbol, m.name
into temp table rna
from mrk_marker m, gxd_htsample_rnaseq rna
where rna._marker_key = m._marker_key
and m._marker_status_key = 2
union
select distinct m._marker_key, m._marker_status_key, m.symbol, m.name
from mrk_marker m, gxd_htsample_rnaseqcombined rna
where rna._marker_key = m._marker_key
and m._marker_status_key = 2
;
create index idx1 on rna(_marker_key);
select count(distinct _marker_key) from rna;

-- get "new" marker key
select distinct r.*, m._marker_key as newMarkerKey, m.symbol as newSymbol
into temp table newrna
from rna r, mrk_marker m, mrk_history h
where r._marker_key = h._history_key
and h._marker_key = m._marker_key
and r.name != 'withdrawn'
;
create index idx2 on newrna(_marker_key);
select * from newrna;
select count(rna.*) from gxd_htsample_rnaseq rna, newrna n where rna._marker_key = n._marker_key;
select distinct n.*, m._marker_key, m.symbol, m.name
from gxd_htsample_rnaseq rna, newrna n, mrk_marker m
where rna._marker_key = n._marker_key
and n.newMarkerKey = m._marker_key
;
update gxd_htsample_rnaseq rna
set _marker_key = newrna.newMarkerkey
from newrna
where rna._marker_key = newrna._marker_key
;
update gxd_htsample_rnaseqcombined rna
set _marker_key = newrna.newMarkerkey
from newrna
where rna._marker_key = newrna._marker_key
;
select rna.*, m._marker_key, m.symbol, m.name
from gxd_htsample_rnaseq rna, newrna n, mrk_marker m
where rna._marker_key = n._marker_key
and n.newMarkerKey = m._marker_key
;

select distinct r.*, m._marker_key as newMarkerKey, m.symbol as newSymbol
into temp table nonewrna
from rna r, mrk_marker m, mrk_history h
where r._marker_key = h._history_key
and h._marker_key = m._marker_key
and r.name = 'withdrawn'
;
create index idx3 on nonewrna(_marker_key);
select * from nonewrna;
select distinct n.*, m._marker_key, m.symbol, m.name
from gxd_htsample_rnaseq rna, nonewrna n, mrk_marker m
where rna._marker_key = n._marker_key
and n.newMarkerKey = m._marker_key
;
delete from gxd_htsample_rnaseq
using nonewrna
where gxd_htsample_rnaseq._marker_key = nonewrna._marker_key
;
delete from gxd_htsample_rnaseqcombined
using nonewrna
where gxd_htsample_rnaseqcombined._marker_key = nonewrna._marker_key
;

-- if all went well, then counts should be 0
select count(m.*)
from mrk_marker m, gxd_htsample_rnaseq rna
where rna._marker_key = m._marker_key
and m._marker_status_key = 2
;

select count(m.*)
from mrk_marker m, gxd_htsample_rnaseqcombined rna
where rna._marker_key = m._marker_key
and m._marker_status_key = 2
;

EOSQL

date |tee -a $LOG

