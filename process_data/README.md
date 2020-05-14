The scripts in this directory perform the following preprocessing steps:

They scripts used here require that you create a local anaconda environment using the included `env.yml` file.

```{sh}
conda env create -n anno -f env.yml
source activate anno
```

1. Identify TOPMed subsamples to use, defined by global ancestry fraction (`define_samples.R`--requires various sample summary files)
2. Taking the full TOPMed VCF files as input, extract singletons within each ancestry subsample, annotate each singleton with the sequence context and mutation type, and output singleton data to a tab-delimited text file (`get_chr_singletons.1.py`—parallelized through use of `get_chr_singletions.slurm`)
3. Calculate the distance to the nearest singleton in each individual and add this as a new column in the text file (`sort_add_clid.py`—parallelized through use of `sort_add_clid.slurm`)
4. Finally, `merge_singletons.sh` merges the per-chromosome files to a single input file, to be read by `../scripts/analyze_clusters.R`.