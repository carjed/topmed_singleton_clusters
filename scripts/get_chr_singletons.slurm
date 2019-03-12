#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jedidiah@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=2000
#SBATCH --time 24:00:00
#SBATCH --job-name=get_topmed_freeze5_singletons
#SBATCH --partition=bipolar
#SBATCH --array=1-300
#SBATCH --requeue
#SBATCH --output="/net/snowwhite/home/jedidiah/slurm/slurmJob-%J.out" --error="/net/snowwhite/home/jedidiah/slurm/slurmJob-%J.err"
# s#SBATCH --export=PATH=/net/snowwhite/home/jedidiah/anaconda3/envs/anno/bin:$PATH
# source /net/snowwhite/home/jedidiah/anaconda3/bin/activate /net/snowwhite/home/jedidiah/anaconda3/envs/anno

# set paths
POP=eur
PYTHONPATH=/net/snowwhite/home/jedidiah/anaconda3/envs/anno/bin/python
# SAMPLEFILE=/net/snowwhite/home/jedidiah/sandbox_data/${POP}_1k.txt
SAMPLEFILE=/net/snowwhite/home/jedidiah/sandbox_data/samples_55k_1k_${POP}.txt
SCRIPTDIR=/net/snowwhite/home/jedidiah/scripts
OUTDIR=/net/snowwhite/home/jedidiah/topmed_freeze5_singletons_mask

# Parse input parameters
WINDOWFILE=/net/snowwhite/home/jedidiah/sandbox_data/hg38.autosomes.10mb.sorted.bed
LINE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $WINDOWFILE)

CHR=$(echo $LINE | cut -d ' ' -f 1 | sed s/chr//)
BEGIN=$(echo $LINE | cut -d ' ' -f 2)
END=$(echo $LINE | cut -d ' ' -f 3)

# confirm outputs
echo $CHR
echo $BEGIN
echo $END

# run with params
srun ${PYTHONPATH} ${SCRIPTDIR}/get_chr_singletons.1.py --samplefile ${SAMPLEFILE} --chr ${CHR} --begin ${BEGIN} --end ${END} --outdir ${OUTDIR}/${POP} --mask
