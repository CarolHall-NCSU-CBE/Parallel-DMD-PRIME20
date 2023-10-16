#!/bin/tcsh
#BSUB -n 8
#BSUB -J Parallel_userfriendly
#BSUB -W 24:00
#BSUB -R "span[hosts=1]"
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB standard_ib

module load PrgEnv-intel

rm -f results/*
rm -f outputs/*
rm -f parameters/*
rm -f inputs/*
rm -f checks/*

# Generate initial configuration for new simulation:
/share/hall2/vtnguye4/Finalcleanversion101223/src/initconfig

#Annealing:	
foreach i (`seq 1 9`)
mpirun /share/hall2/vtnguye4/Finalcleanversion101223/src/DMDPRIME20 < inputs/annealtemp_$i > outputs/out_annealtemp_$i
end

#DMD simulations
foreach i (`seq 1 20`)
mpirun /share/hall2/vtnguye4/Finalcleanversion101223/src/DMDPRIME20 < inputs/simtemp > outputs/out_simtemp_$i
end
