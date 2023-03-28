#!/bin/bash
source ~/.bashrc
module add c3ddb/glibc/2.14
module add c3ddb/singularity/3.5.2
#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH
singularity exec -B $PWD /scratch/users/anniz44/bin/pro/CRISPRCasFinder/CrisprCasFinder.simg perl /scratch/users/anniz44/bin/pro/CRISPRCasFinder/CRISPRCasFinder.pl -so /scratch/users/anniz44/bin/pro/CRISPRCasFinder/sel392v2.so -cf /scratch/users/anniz44/bin/pro/CRISPRCasFinder/CasFinder-2.0.3 -drpt /scratch/users/anniz44/bin/pro/CRISPRCasFinder/supplementary_files/repeatDirection.tsv -rpts /scratch/users/anniz44/bin/pro/CRISPRCasFinder/supplementary_files/Repeat_List.csv -cas -def G -out /scratch/users/anniz44/bin/pro/CRISPRCasFinder/am_TuSa_g33_test -in /scratch/users/anniz44/bin/pro/CRISPRCasFinder/am_TuSa_g33.fasta
