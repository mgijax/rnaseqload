#rnaseqload

SELECT
    'DROP TABLE IF EXISTS ' || quote_ident(schemaname) || '.' ||
    quote_ident(tablename) || ';'
FROM pg_tables
WHERE schemaname = 'mgd'
  AND tablename LIKE 'tm_temp%'
ORDER BY tablename;