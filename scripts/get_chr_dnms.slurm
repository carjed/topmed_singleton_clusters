#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jedidiah@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=2000
#SBATCH --time 8:00:00
#SBATCH --job-name=get_topmed_freeze3_singletons
#SBATCH --partition=bipolar
#SBATCH --array=1-300
#SBATCH --requeue
#SBATCH --output="/net/snowwhite/home/jedidiah/slurm/slurmJob-%J.out" --error="/net/snowwhite/home/jedidiah/slurm/slurmJob-%J.err"
# s#SBATCH --export=PATH=/net/snowwhite/home/jedidiah/anaconda3/envs/anno/bin:$PATH
# source /net/snowwhite/home/jedidiah/anaconda3/bin/activate /net/snowwhite/home/jedidiah/anaconda3/envs/anno

# set paths
PYTHONPATH=/net/snowwhite/home/jedidiah/anaconda3/envs/anno/bin/python
SAMPLEFILE=/net/snowwhite/home/jedidiah/topmed_freeze5_trio_children_nr.txt
SCRIPTDIR=/net/snowwhite/home/jedidiah/scripts

# Parse input parameters
WINDOWFILE=/net/snowwhite/home/jedidiah/sandbox_data/GRCh37.10mb.sorted.bed
LINE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $WINDOWFILE)

CHR=$(echo $LINE | cut -d ' ' -f 1 | sed s/chr//)
BEGIN=$(echo $LINE | cut -d ' ' -f 2)
END=$(echo $LINE | cut -d ' ' -f 3)

# confirm outputs
echo $CHR
echo $BEGIN
echo $END

# run with params
srun ${PYTHONPATH} ${SCRIPTDIR}/get_chr_dnms.py --samplefile ${SAMPLEFILE} --chr ${CHR} --begin ${BEGIN} --end ${END} --outdir /net/snowwhite/home/jedidiah/topmed_freeze5_dnms
