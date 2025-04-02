#!/bin/bash

# === SETTINGS ===
total_genes=20000 ## change as per your data
batch_size=1000 ## change as per your data
r_script="run_susie_batch.R"
log_dir="logs"
mkdir -p "$log_dir"

# === BATCH SUBMISSION LOOP ===
for ((start=1; start<=total_genes; start+=batch_size)); do
  end=$((start + batch_size - 1))
  if [ "$end" -gt "$total_genes" ]; then end=$total_genes; fi

  echo "Submitting SuSiE job for gene range: ${start}-${end}"

  sbatch <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -p compute ## change name as per your HPC queue system
#SBATCH -J Prc_${start}_${end}
#SBATCH --mem=12G
#SBATCH -c 2 ## change as per your HPC resources.
#SBATCH -t 96:00:00
#SBATCH -o "${log_dir}/SuSiE_${start}_${end}_%j.out"
#SBATCH -e "${log_dir}/SuSiE_${start}_${end}_%j.err"

export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
export NUMEXPR_NUM_THREADS=2
export VECLIB_MAXIMUM_THREADS=2
export R_BLAS_USE_FLOCK=TRUE
#######################################################
cd ~/EQTL/TMM_log2 ; ### change as per your folder strutures
#######################################################
echo "Starting SuSiE batch: genes ${start} to ${end}"
Rscript "$r_script" "$start" "$end"
echo "Finished SuSiE batch: genes ${start} to ${end}"
EOF

done

