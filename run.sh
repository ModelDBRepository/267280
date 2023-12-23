#!/bin/bash 
nodes=1
ppn=10
let nmpi=$nodes*$ppn
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=${ppn}]"
#BSUB -R "select[type=any]"
##BSUB -R "select[cpus=24]"
#BSUB -n ${nmpi}
#BSUB -q ppc_24h 
##BSUB -q p8_6h 
##BSUB -q ppc_6h 
#BSUB -W 1440 
##BSUB -W 360  
#---------------------------------
mpirun --tag-output --bind-to core -np ${nmpi} ./piriformENDO
EOF
#---------------------------------------
bsub  <batch.job
