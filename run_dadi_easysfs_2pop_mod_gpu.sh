#!/bin/bash
#SBATCH --job-name dadi_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=wrosenth@uwyo.edu
#SBATCH --account=ysctrout
#SBATCH --partition=teton-gpu
#SBATCH --gres=gpu:2
#SBATCH --output=dadi_uncerts_0.3_mod_nomaf_%A.out

echo "SLURM_JOB_ID:" $SLURM_JOB_ID
echo "SLURM_JOB_NAME:" $SLURM_JOB_NAME
echo "SLURM_JOB_PARTITION" $SLURM_JOB_PARTITION
echo "SLURM_JOB_NUM_NODES:" $SLURM_JOB_NUM_NODES
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "SLURM_JOB_CPUS_PER_NODE:" $SLURM_JOB_CPUS_PER_NODE
echo "SLURM_TASKS_PER_NODE:" $SLURM_TASKS_PER_NODE
echo "SLURM_CPUS_PER_TASK:" $SLURM_CPUS_PER_TASK
echo "SLURM_CPUS_ON_NODE:" $SLURM_CPUS_ON_NODE
echo "OMP_NUM_THREADS:" $OMP_NUM_THREADS

module load gcc/11.2.0
module load miniconda3/23.1.0
module load cuda/11.8.0

srun nvidia-smi -L

# Path to installed conda environment
conda activate entropy

start=$(date +'%D %T')
echo "Start:" $start

time srun python3 hypoth_bootstrap_2pop_mod_easysfs.py

end=$(date +'%D %T')
echo "End:" $end

start_secs=$(date --date="$start" '+%s')
end_secs=$(date --date="$end"   '+%s')
duration=$((end_secs - start_secs))

echo "Duration:" $duration"sec"
conda deactivate
echo "Done."