#!/bin/bash
#SBATCH --job-name=bench_B3
#SBATCH --account=f202500010hpcvlabuminhox
#SBATCH --partition=normal-x86
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=48:00:00
#SBATCH --array=1-3%1
#SBATCH --output=logs/benchmark_B3_%A_%a.out
#SBATCH --error=logs/benchmark_B3_%A_%a.err

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

SCENARIOS=(
  nearly_gaussian
  symmetric_heavy_tails
  contaminated_normal
)

SCENARIO="${SCENARIOS[$((SLURM_ARRAY_TASK_ID-1))]}"

export PROFILE=final
export REPS_PER_SCEN=50
export N_PER_DATASET=1000
export MASTER_SEED=20260601
export ROOT_DIR=/projects/F202500010HPCVLABUMINHO/cecilia/ese_heckman/results/reviewer_heckmanem_comparison_B_final_seed20260601
export RESUME_EXISTING=TRUE
export AUTO_INSTALL=FALSE
export MAX_WORKERS=50
export SCENARIO_FILTER="${SCENARIO}"

echo "===================================================================="
echo "benchmark B rerun3 scenario array"
echo "job id       : ${SLURM_ARRAY_JOB_ID}"
echo "task id      : ${SLURM_ARRAY_TASK_ID}"
echo "scenario     : ${SCENARIO_FILTER}"
echo "node         : $(hostname)"
echo "cpus         : ${SLURM_CPUS_PER_TASK}"
echo "root         : ${ROOT_DIR}"
echo "start time   : $(date)"
echo "===================================================================="

Rscript --version
Rscript scripts/reviewer_heckmanem_comparison_efficient_corrected.R

echo "===================================================================="
echo "completed scenario: ${SCENARIO_FILTER}"
echo "end time           : $(date)"
echo "===================================================================="
