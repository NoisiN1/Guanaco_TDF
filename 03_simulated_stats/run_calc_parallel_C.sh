#!/bin/bash
# Uso: ./run_calc_parallel_C.sh <1|2>
scen=$1
if [ -z "$scen" ]; then echo "Uso: $0 <1|2>"; exit 1; fi

cd /home/apena/Splatche/datasetC
mkdir -p chunks_scen${scen}

# Chunks de 100: [1-100], [101-200], ..., [9901-10000]  → 100 chunks totales
seq 1 100 9901 | parallel -j 90 \
  "Rscript calc_stats_C.R $scen {} \$(({} + 99)) > logs/calc_scen${scen}_chunk_{}.log 2>&1"

# Combinar chunks en el archivo final sim<scen>_abc.txt
echo -e "simid\tK\tm\tH\tFST" > sim${scen}_abc_C.txt
for f in chunks_scen${scen}/chunk_*.txt; do
  tail -n +2 "$f" >> sim${scen}_abc_C.txt
done

echo "Total simulaciones procesadas:"
wc -l sim${scen}_abc_C.txt
