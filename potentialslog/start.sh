#!/bin/bash

#SBATCH -D /lustre/nyx/hades/user/vnikolae/SimData/SMASH/TMP/
#SBATCH -J smash_gen
#SBATCH -p main
#SBATCH --time=05:00:00
#SBATCH -a 1-1
#SBATCH -o /lustre/nyx/hades/user/vnikolae/SimData/SMASH/SGE_OUT/slurm_%A_%a.out
#SBATCH -e /lustre/nyx/hades/user/vnikolae/SimData/SMASH/SGE_OUT/slurm_%A_%a.err

#Energy in cms
ecm=5.0
#Number of events
nev=100

export COMMIT=SMASH_ROOT_11.5gev_test

export START_POSITION=$PWD

export MAIN_DIR=/lustre/nyx/hades/user/$USER/SimData/SMASH

export INPUTFILE=$MAIN_DIR/inputfiles/config.yaml

export DATE=${SLURM_ARRAY_JOB_ID} # `date '+%Y%m%d_%H%M%S'`
# export DATE=`date '+%Y%m%d_%H%M%S'`

export OUT=$MAIN_DIR/OUT/$COMMIT/$DATE
export OUT_LOG=$OUT/log
export OUT_FILE=$OUT/files

mkdir -p $OUT
mkdir -p $OUT_LOG
mkdir -p $OUT_FILE

_log() {

local format='+%Y/%m/%d-%H:%M:%S'
echo [`date $format`] "$@"

}


ROOTSYS_cvmfs=/cvmfs/it.gsi.de/root/v6-06-06/

export ROOTSYS=${ROOTSYS_cvmfs}
export PATH=${PATH}:${ROOTSYS}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib/root:${ROOTSYS}/lib
export Root=${ROOTSYS}/bin/root

# export McDST_DIR=/lustre/nyx/hades/user/parfenov/Soft/McDst
# source $McDST_DIR/mcdst_environment.sh

export SMASH_DIR=/lustre/nyx/hades/user/parfenov/Soft/Models/SMASH/smash/build
export SMASH_BIN=${SMASH_DIR}/smash

export TMPALL=$MAIN_DIR/TMP
export TMPDIR=$TMPALL/TMP_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR_OUT=$TMPDIR/OUT

mkdir -p $TMPDIR_OUT

cd $TMPDIR

cp $INPUTFILE ${TMPDIR}/config.yaml

sed -e "s|energyincms|$ecm|" -i ./config.yaml
sed -e "s|numberofevents|$nev|" -i ./config.yaml
sed -e "s|randomrandom|`shuf -i 1-1000000 -n 1`|" -i ./config.yaml


_log ${TMPDIR}

_log ${INPUTFILE}

_log `ls ${TMPDIR}`

_log ${ROOTSYS_cvmfs}

_log ${ROOTSYS}

_log ${LD_LIBRARY_PATH}

_log ${PATH}

#cat inputfile >> $OUT/JOB_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log
cat config.yaml >> $OUT_LOG/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

_log -------

_log "Running SMASH..."
$SMASH_BIN -i ${TMPDIR}/config.yaml -o ${TMPDIR_OUT}/ &>>$OUT_LOG/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

_log -------

cd $TMPDIR
_log `ls $TMPDIR`

_log -------

_log Moving output files from $TMPDIR_OUT to ${OUT_FILE}...
mv ${TMPDIR_OUT}/Particles.root $OUT_FILE/particles_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root
mv ${TMPDIR_OUT}/Collisions.root $OUT_FILE/collisions_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root

_log Cleaning temporary directory...
rm -rf $TMPDIR

cd $START_POSITION
_log "Job is done!"

