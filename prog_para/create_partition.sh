#!/bin/bash

NB_PARTITIONS=4
NB_VERTEX=5768

rm -rf metis scotch
mkdir metis scotch
for i in `seq 2 $NB_PARTITIONS`;
do
    mkdir metis/partition$i
    mkdir scotch/partition$i
done 


#metis

for i in `seq 1 $NB_PARTITIONS`;
do
    cd TP_METIS_SCOTCH
    if [ $i -eq 1 ]
    then
	for j in `seq 1 $NB_VERTEX`;
	do
	    echo 0 >> dualformetis.dat.epart.1
	done
    else
	metis-4.0.3/partdmesh dualformetis.dat $i
    fi
    mv dualformetis.dat.epart.$i ../TP_PROG_EF
    cd ../TP_PROG_EF
    cat meshprogc.data dualformetis.dat.epart.$i > mesh_for_progc.data
    gcc -o Preprocess.exe Preprocess.c
    ./Preprocess.exe $i mesh_for_progc.data
    mpicc -std=c99 -lm -o fem.exe FemPar.c
    cd ..
    rsync --delete -r TP_PROG_EF/* metis/partition$i
    rm TP_PROG_EF/Data0*
    rm TP_PROG_EF/dualformetis.dat.epart.*
done 


