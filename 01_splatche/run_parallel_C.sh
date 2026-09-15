#!/bin/bash
# Uso: ./run_parallel_C.sh <scenario_number>
scen=$1
if [ -z "$scen" ]; then
  echo "Uso: $0 <1|2>"
  exit 1
fi

homedir=$(pwd)

seq 1 10000 | parallel -j 96 \
  "cd $homedir/scenario${scen}_C/subscen{} && bash $homedir/runSplat_C.sh settings.txt {} $scen"
