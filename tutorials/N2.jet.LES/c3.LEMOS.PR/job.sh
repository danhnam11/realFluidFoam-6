#!/bin/bash

# this is comment line

#parallel environment request

#$ -pe mpi_64 64

# our Job name 
#$ -N N2.PRc3Smar

#$ -q cpl3.q
#$ -S /bin/bash

##$ -l h=!(cpl2-01)
##$ -l h=(cpl2-03-ib)

#$ -cwd


export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="hfi1_0:1"


source $HOME/OpenFOAM/OpenFOAM-6.x/etc/bashrc WM_COMPILER_TYPE=ThirdParty WM_COMPILER=Gcc48 WM_LABEL_SIZE=64 WM_MPLIB=OPENMPI FOAMY_HEX_MESH=yes

export MPI_EXEC=~/OpenFOAM/ThirdParty-6.x/platforms/linux64Gcc48/openmpi-2.1.1/bin/mpirun

solver=realFluidFoam

$MPI_EXEC -np $NSLOTS $solver -parallel > run_out.log

