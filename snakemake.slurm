#!/bin/bash -login
#SBATCH -p bmm                 # partition, or queue, to assign to
#SBATCH -J charcoal            # name for job
#SBATCH -N 1                   # one "node", or computer
#SBATCH -n 1                   # one task for this node
#SBATCH -c 16                  # cores per task
#SBATCH -t 1-0                 # ask for a full day
#SBATCH --mem=300000             # memory (300,000 mb = 300gb)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=titus@idyll.org

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate your desired conda environment
conda activate sgc

# fail on weird errors
set -o nounset
set -o errexit
set -x

# go to the directory you ran 'sbatch' in, OR just hardcode it...
#cd $SLURM_SUBMIT_DIR
cd ~/charcoal

# run the snakemake!
snakemake -p -j 16 --use-conda --configfile=conf/ibd2.conf --resources mem=${SLURM_MEM_PER_NODE:-150000} -n -- all reports all_make_tree_viz

# print out various information about the job
env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
