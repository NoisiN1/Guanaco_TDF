#!/bin/bash
sub=$1; scen=$2
out="/home/apena/Splatche/datasetC/calib_${scen}.txt"
home=$(pwd)
tmp="/tmp/cal_$$_${scen}_${sub}"
mkdir -p $tmp; rsync -a $home/* $tmp/ 2>/dev/null
cd $tmp; chmod +x ./SPLATCHE3-Linux-64b
./SPLATCHE3-Linux-64b settings.txt > /dev/null 2>&1
G=datasets/GeneticsOutput
tot=0; ven=0
for f in $G/*Immigrants_*_P1.nm; do
  [ -f "$f" ] || continue
  cel=$(basename $f | sed 's/.*Immigrants_\([0-9]*\)_P1.nm/\1/')
  [ "$cel" = "558" ] && continue
  read a b <<< $(awk 'NR>1{s+=$2; if($1>=560 && $1<=760) v+=$2} END{print s+0, v+0}' "$f")
  tot=$((tot+a)); ven=$((ven+b))
done
m=$(grep -oP 'MigrationRate=\K[0-9.]+' settings.txt)
echo -e "${sub}\t${m}\t${tot}\t${ven}" >> $out
cd $home; rm -rf $tmp
