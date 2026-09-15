#!/bin/bash
mkdir -p scenario2_C
for i in $(seq 1 10000)
do
  echo "$i"
  mkdir -p ./scenario2_C/subscen$i
  k=$(sed -n "${i}p" param.txt | cut -d" " -f4)
  m=$(sed -n "${i}p" param.txt | cut -d" " -f2)
  cp -r ./Splatche_scen2_hum_TDF_C/* ./scenario2_C/subscen${i}/.
  sed -i "s/PARAM/$k/g" ./scenario2_C/subscen${i}/datasets/veg2K_pop1_time_5.txt
  sed -i "s/PARAM/$k/g" ./scenario2_C/subscen${i}/datasets/veg2K_pop1_time_6.txt
  sed -i "s/MigrationRate=0.01/MigrationRate=$m/g" ./scenario2_C/subscen${i}/settings.txt
done
