export HSC_SSP_CAS_PASSWORD={your stars account password}
for f in /your/path/to/sqls/*.sql; do
  base=$(basename "$f" .sql)
  echo "Running $f..."
  python3 hscSspQuery3.py --user kongyouwen --release-version dr4 --format fits "$f"
done
