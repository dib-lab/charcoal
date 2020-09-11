#!/bin/bash -login
#SBATCH -p med2                # partition, or queue, to assign to
#SBATCH -J charcoal            # name for job
#SBATCH -N 1                   # one "node", or computer
#SBATCH -n 1                   # one task for this node
#SBATCH -c 32                  # cores per task
#SBATCH -t 2:00:00                 # ask for an hour
#SBATCH --mem=60000             # memory (30,000 mb = 30gb)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=titus@idyll.org

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate your desired conda environment
conda activate charcoal

# fail on weird errors
set -o nounset
set -o errexit
set -x

# go to the directory you ran 'sbatch' in, OR just hardcode it...
#cd $SLURM_SUBMIT_DIR
cd ~/charcoal

# run the snakemake!
#charcoal run conf/gtdb-random-dna.conf -p -j 32 -k --no-use-conda
#charcoal run conf/gtdb-contam-dna.conf -p -j 32 -k --no-use-conda
#charcoal run eval/tara-delmont-relaxed.conf -p -j 32 -k --no-use-conda --unlock
#charcoal run eval/tara-delmont-relaxed.conf -p -j 32 -k --no-use-conda
#charcoal run eval/tara-delmont-strict.conf -p -j 32 -k --no-use-conda --unlock
#charcoal run eval/tara-delmont-strict.conf -p -j 32 -k --no-use-conda
python -m charcoal run conf/ibd2.conf -j 32 -k -p --unlock
python -m charcoal run conf/gtdb-contam-dna.conf -j 32 -k -p --unlock
python -m charcoal run conf/tara-delmont.conf -j 32 -k -p --unlock

python -m charcoal run conf/ibd2.conf -j 32 -k -p
python -m charcoal run conf/gtdb-contam-dna.conf -j 32 -k -p
python -m charcoal run conf/tara-delmont.conf -j 32 -k -p
#charcoal run conf/gtdb-random-dna.conf -j 32 -k -p

# print out various information about the job
env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
