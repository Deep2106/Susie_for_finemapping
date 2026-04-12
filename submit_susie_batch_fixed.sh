#!/bin/bash
# =============================================================================
# submit_susie_batch_fixed.sh
# Submits SuSiE finemapping SLURM jobs for eGenes only (~1,000-5,000 genes).
#
# Pipeline order:
#   1. Rscript filter_egenes.R                             → gene.tsv
#   2. Rscript calc_zscore_fixed.R <cis_all_INFO.RDS>      → per-gene .fminput.z
#   3. bash submit_susie_batch_fixed.sh                    → SLURM jobs
#   4. Rscript combine_pip_fixed.R                         → merged outputs
# =============================================================================

# ---------------------------------------------------------------------------
# Settings — replace placeholders before running
# ---------------------------------------------------------------------------
R_SCRIPT="run_susie_batch_fixed.R"
GENE_INFO="gene.tsv"       # eGenes only — from filter_egenes.R
                                     # NOT the full 17K gene list
WORK_DIR="/path/to/CIS_TRANS_EQTLS_SEPT_2025/COVS_default"
LOG_DIR="logs"

BATCH_SIZE=200    # eGenes per job
                  # 200 genes × 180 sec/gene worst case = 10h → safe in 48h
MEM="15G"
CORES=2
WALLTIME="48:00:00"
PARTITION="compute"

# ---------------------------------------------------------------------------
# Conda — using base environment
# R and all required packages (susieR, seqminer, data.table,
# MASS, Matrix, filelock) must be installed in base.
# If you ever move to a named environment, replace the activation
# block below with:
#   source "$(conda info --base)/etc/profile.d/conda.sh"
#   conda activate YOUR_ENV_NAME
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Validate prerequisites
# ---------------------------------------------------------------------------
GENE_FILE="${WORK_DIR}/${GENE_INFO}"

if [ ! -f "${GENE_FILE}" ]; then
  echo "ERROR: ${GENE_FILE} not found."
  echo "  Run filter_egenes.R first to generate gene_remaining.tsv"
  exit 1
fi

if [ ! -f "${WORK_DIR}/${R_SCRIPT}" ]; then
  echo "ERROR: ${WORK_DIR}/${R_SCRIPT} not found."
  exit 1
fi

# ---------------------------------------------------------------------------
# Derive TOTAL_GENES from gene_remaining.tsv — never hardcode
# ---------------------------------------------------------------------------
TOTAL_GENES=$(( $(wc -l < "${GENE_FILE}") - 1 ))   # subtract header
N_JOBS=$(( (TOTAL_GENES + BATCH_SIZE - 1) / BATCH_SIZE ))
MAX_H=$(( BATCH_SIZE * 180 / 3600 ))

echo "========================================"
echo " SuSiE finemapping submission"
echo "========================================"
echo " eGenes to process : ${TOTAL_GENES}"
echo " Batch size        : ${BATCH_SIZE}"
echo " Jobs to submit    : ${N_JOBS}"
echo " Max time/job      : ~${MAX_H}h (at 180s/gene)"
echo " Wall time limit   : ${WALLTIME}"
echo " Memory/job        : ${MEM}"
echo " Conda env         : base (default)"
echo "========================================"
echo ""

if [ "$TOTAL_GENES" -gt 10000 ]; then
  echo "WARNING: TOTAL_GENES=${TOTAL_GENES} is unexpectedly large."
  echo "  Confirm gene_remaining.tsv was produced by filter_egenes.R"
  echo "  and is NOT the full 17K gene coordinate file."
  echo ""
fi

mkdir -p "${WORK_DIR}/${LOG_DIR}"

# ---------------------------------------------------------------------------
# Submission loop
# ---------------------------------------------------------------------------
JOB_COUNT=0

for (( START=1; START<=TOTAL_GENES; START+=BATCH_SIZE )); do
  END=$(( START + BATCH_SIZE - 1 ))
  [ "$END" -gt "$TOTAL_GENES" ] && END=$TOTAL_GENES

  echo "Submitting genes ${START}–${END}"

  sbatch <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -p ${PARTITION}
#SBATCH -J SuSiE_${START}_${END}
#SBATCH --mem=${MEM}
#SBATCH -c ${CORES}
#SBATCH -t ${WALLTIME}
#SBATCH -o ${WORK_DIR}/${LOG_DIR}/SuSiE_${START}_${END}_%j.out
#SBATCH -e ${WORK_DIR}/${LOG_DIR}/SuSiE_${START}_${END}_%j.err

# Using base conda environment — no activation needed
# All required R packages must be installed in base:
#   susieR, seqminer, data.table, dplyr, MASS, Matrix, filelock
echo "Conda env  : \${CONDA_DEFAULT_ENV:-base}"
echo "R version  : \$(R --version | head -1)"

# Limit BLAS/LAPACK threads to allocated cores
# SuSiE is single-threaded; these affect eigen() and nearPD() only
export OMP_NUM_THREADS=${CORES}
export OPENBLAS_NUM_THREADS=${CORES}
export MKL_NUM_THREADS=${CORES}
export NUMEXPR_NUM_THREADS=${CORES}
export VECLIB_MAXIMUM_THREADS=${CORES}

cd ${WORK_DIR}

echo "Host     : \$(hostname)"
echo "Started  : \$(date)"
echo "Genes    : ${START} to ${END}"

Rscript ${R_SCRIPT} ${START} ${END}
EXIT_CODE=\$?

echo "Finished : \$(date)"
echo "Exit code: \${EXIT_CODE}"
exit \${EXIT_CODE}
EOF

  JOB_COUNT=$(( JOB_COUNT + 1 ))
done

echo ""
echo "========================================"
echo " Submitted ${JOB_COUNT} jobs"
echo "========================================"
echo ""
echo "Monitor:"
echo "  squeue -u \$USER"
echo "  watch -n 30 squeue -u \$USER"
echo ""
echo "Check for failures after completion:"
echo "  grep -h 'ERROR\|NOT_CONVERGED' ${WORK_DIR}/${LOG_DIR}/gene_runtime_summary_susie.txt"
echo "  grep -h 'ERROR' ${WORK_DIR}/${LOG_DIR}/error_log_susie.txt"
echo ""
echo "After all jobs finish — Step 4:"
echo "  Rscript combine_pip_fixed.R"
