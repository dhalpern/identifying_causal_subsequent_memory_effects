#!/bin/bash
#
#SBATCH --job-name=glmnet
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=20GB
#SBATCH --mail-type=END
#SBATCH --mail-user=david.halpern@nyu.edu
#SBATCH --output=slurm_glmnet_%A_%a.out

module purge
module load r/gcc/4.0.4

SRCDIR=$HOME/identifying_causal_subsequent_memory_effects/models/

Rscript --vanilla \
$SRCDIR/run_models.R $SLURM_ARRAY_TASK_ID irt_test_loc

