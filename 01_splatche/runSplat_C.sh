#!/bin/bash
file=$1
subscen=$2
scen=$3
outdir="/home/apena/Splatche/datasetC/outscenario${scen}_C"
homedir=$(pwd)
NOW=$(date '+%F_%H_%M_%S')
tmpdir="tmp_${USER}_${NOW}_${scen}_subscen${subscen}"
mkdir /tmp/$tmpdir
rsync -a $homedir/* /tmp/$tmpdir/.
cd /tmp/$tmpdir
chmod +x ./SPLATCHE3-Linux-64b
./SPLATCHE3-Linux-64b $file > /dev/null 2>&1
rsync -a datasets/GeneticsOutput/settings_GeneSamples_C_1.arp $outdir/out${scen}_${subscen}.arp
cd $homedir
rm -r /tmp/$tmpdir
