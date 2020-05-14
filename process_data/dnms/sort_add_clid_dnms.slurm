#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jedidiah@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=2000
#SBATCH --time 1:00:00
#SBATCH --job-name=add_clid
#SBATCH --partition=bipolar
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH --output="/net/snowwhite/home/jedidiah/slurm/slurmJob-%J.out" --error="/net/snowwhite/home/jedidiah/slurm/slurmJob-%J.err"
# s#SBATCH --export=PATH=/net/snowwhite/home/jedidiah/anaconda3/envs/anno/bin:$PATH
# source /net/snowwhite/home/jedidiah/anaconda3/bin/activate /net/snowwhite/home/jedidiah/anaconda3/envs/anno

# set paths
POP=afr
PYTHONPATH=/net/snowwhite/home/jedidiah/anaconda3/envs/anno/bin/python
SCRIPTDIR=/net/snowwhite/home/jedidiah/scripts

# run with params
srun ${PYTHONPATH} ${SCRIPTDIR}/sort_add_clid_dnms.py --chr ${SLURM_ARRAY_TASK_ID}
