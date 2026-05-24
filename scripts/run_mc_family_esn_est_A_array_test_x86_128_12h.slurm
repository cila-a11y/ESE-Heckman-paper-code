#!/bin/bash
#SBATCH --job-name=mc_A_test
#SBATCH --account=f202500010hpcvlabuminhox
#SBATCH --partition=normal-x86
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=12:00:00
#SBATCH --array=1-4%2
#SBATCH --output=logs/mc_A_test_%A_%a.out
#SBATCH --error=logs/mc_A_test_%A_%a.err

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
export RUN_TAG=array_test
export MASTER_SEED=20260621
export ROOT_DIR=/projects/F202500010HPCVLABUMINHO/cecilia/ese_heckman/results/mc_family_esn_est_A_array_test_seed20260621
export RESUME_EXISTING=TRUE
export AUTO_INSTALL=FALSE

export REPS_GRID=2
export REPS_NSWEEP=2
export N_GRID=300
export N_SWEEP=500,1000
export CDF_MAXPTS=5000
export CDF_ABSEPS=1e-5
export CDF_RELEPS=1e-5
export OPT_MAXIT=300
export REPL_TIMEOUT=3600

printf '============================================================\n'
printf 'A array test | job %s | task %s | node %s\n' "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "$(hostname)"
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
