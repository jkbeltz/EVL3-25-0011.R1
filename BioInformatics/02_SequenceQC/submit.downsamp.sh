#!/bin/bash

while read line ; do 
    set $line
    sbatch ./02I_DownsampBams.orch2020.sbatch $1 $2
done <./Orch2020.DownsampParams.txt 
