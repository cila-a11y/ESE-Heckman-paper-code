#!/bin/bash
#SBATCH --job-name=mc_A_array
#SBATCH --account=f202500010hpcvlabuminhox
#SBATCH --partition=normal-x86
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=48:00:00
#SBATCH --array=1-66%4
#SBATCH --output=logs/mc_A_array_%A_%a.out
#SBATCH --error=logs/mc_A_array_%A_%a.err

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
export RUN_TAG=array_final100
export MASTER_SEED=20260611
export ROOT_DIR=/projects/F202500010HPCVLABUMINHO/cecilia/ese_heckman/results/mc_family_esn_est_A_array_final_seed20260611
export RESUME_EXISTING=TRUE
export AUTO_INSTALL=FALSE

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

export CDF_MAXPTS=100000
export CDF_ABSEPS=1e-7
export CDF_RELEPS=1e-7
export OPT_MAXIT=2000
export REPL_TIMEOUT=43200

printf '============================================================\n'
printf 'A array final | job %s | task %s | node %s\n' "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "$(hostname)"
printf 'root: %s\n' "${ROOT_DIR}"
printf '============================================================\n'
Rscript --version
# array precreate shared output directories
mkdir -p "${ROOT_DIR}/results"
mkdir -p "${ROOT_DIR}/report/tables"
mkdir -p "${ROOT_DIR}/report/figures"
mkdir -p "${ROOT_DIR}/report/metadata"
mkdir -p "${ROOT_DIR}/report/objects"
chmod -R u+rwX,g+rwX "${ROOT_DIR}"
Rscript scripts/mc_family_esn_est_study_array_aligned.R
