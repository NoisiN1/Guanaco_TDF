#!/bin/bash
file=$1; sub=$2; scen=$3
outdir="/home/apena/Splatche/datasetC/calib_scen${scen}"
homedir=$(pwd)
tmp="calib_$$_${scen}_${sub}"
mkdir -p /tmp/$tmp "$outdir"
rsync -a $homedir/* /tmp/$tmp/.
cd /tmp/$tmp
chmod +x ./SPLATCHE3-Linux-64b
./SPLATCHE3-Linux-64b $file > /dev/null 2>&1
cp datasets/Arrival_cell_output.txt "$outdir/arrival_${sub}.txt" 2>/dev/null
cp datasets/GeneticsOutput/settings_GeneSamples_C_1.arp "$outdir/out_${sub}.arp" 2>/dev/null
find . -name "*ensity*" -o -name "*ccupation*" 2>/dev/null | grep -vi bmp | \
  while read f; do cp "$f" "$outdir/$(basename $f)_${sub}" 2>/dev/null; done
cd $homedir; rm -rf /tmp/$tmp
