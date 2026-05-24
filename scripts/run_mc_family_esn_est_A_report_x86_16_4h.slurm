#!/bin/bash
#SBATCH --job-name=mc_A_report
#SBATCH --account=f202500010hpcvlabuminhox
#SBATCH --partition=normal-x86
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=logs/mc_A_report_%j.out
#SBATCH --error=logs/mc_A_report_%j.err

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

export RUN_MODE=report
export RUN_TAG=array_final100
export MASTER_SEED=20260611
export ROOT_DIR=/projects/F202500010HPCVLABUMINHO/cecilia/ese_heckman/results/mc_family_esn_est_A_array_final_seed20260611
export RESUME_EXISTING=TRUE
export AUTO_INSTALL=FALSE

printf '============================================================\n'
printf 'A report | job %s | node %s\n' "${SLURM_JOB_ID}" "$(hostname)"
printf 'root: %s\n' "${ROOT_DIR}"
printf '============================================================\n'
Rscript --version
Rscript scripts/mc_family_esn_est_study_array_aligned.R
