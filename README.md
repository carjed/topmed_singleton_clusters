Modeling singleton distance distributions using exponential mixture models

## Installation

```{sh}
git clone https://github.com/carjed/doomsayer.git
cd doomsayer
```

## Preprocessing

There are several preprocessing steps required to filter and annotate singletons from the input VCF files. Preprocessing scripts and instructions are located in the `process_data` directory.

## Analysis

Once data is prepared, the scripts in the `scripts` directory will fit the exponential mixture models and run downstream analyses, generate plots, etc.