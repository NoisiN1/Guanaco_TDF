#!/bin/bash
mkdir -p scenario1_C
for i in $(seq 1 10000)
do
  echo "$i"
  mkdir -p ./scenario1_C/subscen$i
  k=$(sed -n "${i}p" param.txt | cut -d" " -f3)
  m=$(sed -n "${i}p" param.txt | cut -d" " -f2)
  cp -r ./Splatche_scen1_deme_TDF_C/* ./scenario1_C/subscen${i}/.
  sed -i "s/PARAM/$k/g" ./scenario1_C/subscen${i}/datasets/veg2K_pop1_time_3.txt
  sed -i "s/MigrationRate=0.01/MigrationRate=$m/g" ./scenario1_C/subscen${i}/settings.txt
done
