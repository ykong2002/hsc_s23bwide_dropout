for f in sql_files_s23b/*.sql; do
  base=$(basename "$f" .sql)
  echo "Running $f..."
  python3 hscSspQuery3.py --user kongyouwen --release-version dr4 --format fits "$f"
done


python3 run_sqls.py --user kongyouwen --release-version dr4 --format fits --sql-dir sql_files/

export HSC_SSP_CAS_PASSWORD='Kyw200206045221@'

export HSC_SSP_CAS_PASSWORD='Kyw200206045221@'

for f in sql_files/2*.sql; do
  base=$(basename "$f" .sql)
  echo "Running $f..."
  python3 hscSspQuery3.py --user kongyouwen --release-version dr4 --format fits "$f"
done

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export HSC_SSP_CAS_PASSWORD='Kyw200206045221@'
for f in /home/kongyw/hsc_wide_dropout_2510/vlass/sqls/*.sql; do
  base=$(basename "$f" .sql)
  echo "Running $f..."
  python3 hscSspQuery3.py --user kongyouwen --release-version dr4 --format fits "$f"
done
