#!/bin/bash
#SBATCH --job-name=mc_A_core
#SBATCH --account=f202500010hpcvlabuminhox
#SBATCH --partition=normal-x86
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=48:00:00
#SBATCH --array=19,21,23,55,56,59,60,65,66%1
#SBATCH --output=logs/mc_A_core_%A_%a.out
#SBATCH --error=logs/mc_A_core_%A_%a.err

set -euo pipefail

module purge
module load R/4.4.2-gfbf-2024a

export R_LIBS_USER=/projects/F202500010HPCVLABUMINHO/cecilia/Rlibs/4.4.2

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd /projects/F202500010HPCVLABUMINHO/cecilia/ese_heckman

mkdir -p logs results

export RUN_MODE=scenario
export SCENARIO_INDEX="${SLURM_ARRAY_TASK_ID}"

export MASTER_SEED=20260611
export ROOT_DIR=/projects/F202500010HPCVLABUMINHO/cecilia/ese_heckman/results/mc_family_esn_est_A_array_final_seed20260611

export RUN_GRID_STUDY=TRUE
export RUN_NSWEEP_STUDY=TRUE

export REPS_GRID=100
export REPS_NSWEEP=100
export N_GRID=1500
export N_SWEEP=500,1000,1500,2000,2500,3000

export GRID_RHOS=0.0,0.3,0.7
export GRID_ETA_MAGS=0.3,0.6,1.0
export GRID_SEL_LEVELS=COMMON,MODERATE,RARE
export GRID_TAU_BAR=0.4
export GRID_EST_NU=6

export INCLUDE_OPPOSITE_SIGN_DIRECTION=FALSE
export START_STRATEGY=data_multistart
export RESUME_EXISTING=TRUE
export REPL_TIMEOUT=43200
export MAX_WORKERS=128

echo "===================================================================="
echo "A-core corrected recovery study"
echo "job id       : ${SLURM_ARRAY_JOB_ID}"
echo "task id      : ${SLURM_ARRAY_TASK_ID}"
echo "scenario idx : ${SCENARIO_INDEX}"
echo "node         : $(hostname)"
echo "cpus         : ${SLURM_CPUS_PER_TASK}"
echo "root         : ${ROOT_DIR}"
echo "start time   : $(date)"
echo "===================================================================="

Rscript --version
Rscript scripts/mc_family_esn_est_study_array_aligned.R

echo "===================================================================="
echo "completed A-core scenario idx: ${SCENARIO_INDEX}"
echo "end time                    : $(date)"
echo "===================================================================="
