# 🔬 Susie for Finemapping with SLURM Automation

Welcome to **Susie_for_finemapping** - a powerful and automated solution for batch-wise fine-mapping using SuSiE, built to work seamlessly with SLURM for high-performance computing environments.

---

## 📌 Overview

This script dynamically generates SLURM job scripts to **process genes in batches** for fine-mapping using the **SuSiE (Sum of Single Effects)** model.

---

# 🔬 Susie for Finemapping with SLURM Automation

Welcome to **Susie_for_finemapping** - a powerful and automated solution for batch-wise fine-mapping using SuSiE, built to work seamlessly with SLURM for high-performance computing environments.

---

## 🛠️ Prerequisites

- **R version 4.2 or later**  
  Check your installed R version with:
  ```r
  R.version.string


## 🧠 What It Does

- 📥 Takes **Z-scores** from eQTL analysis pipelines:
  - `MatrixEQTL`
  - `FastQTL`
  - `TensorQTL`

- ⚙️ Incorporates **covariates** to generate a **covariate-adjusted LD matrix**

- 🔬 Performs **SuSiE fine-mapping**, returning:
  - Posterior inclusion probabilities (PIPs)
  - The full SuSiE object
  - (Optional) **Credible sets** for downstream analysis

---

## 🛠️ Features

- ✅ SLURM job generation for large-scale computation
- ✅ Supports multiple eQTL formats
- ✅ Modular and customizable
- ✅ Suitable for HPC/cluster environments

---

## 📂 Output

- **PIP values** per SNP
- **SuSiE object** (for reproducibility and downstream analysis)
- **Credible sets** (optional, for high-confidence variant identification)

---

## 🚀 Coming Soon / To-Do

- [ ] Add support for other LD matrix sources
- [ ] Parallelization improvements
- [ ] Containerization (Docker/Singularity)

---

## 🧬 Citation / Attribution

If you use this script in your research, please cite or credit the repository.

---

Feel free to contribute, suggest improvements, or ask questions. Happy fine-mapping! 🙌




