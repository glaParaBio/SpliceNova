<!-- vim-markdown-toc GFM -->

* [Set up environment](#set-up-environment)
* [Run](#run)

<!-- vim-markdown-toc -->

## Set up environment

```
conda create --yes -n diffsplice
mamba install --freeze-installed -n diffsplice --yes --file requirements.txt
```

## Run

Activate environment:

```
conda activate diffsplice
```

Run pipeline:

```
snakemake -p --dryrun --jobs 1 \
    --config fastqdir=/export/projects/III-data/wcmp_bioinformatics/db291g/data/20190411_kasia_diffsplice/fastq/ \
             sample_sheet=$PWD/sample_sheet.tsv \
             species=$PWD/species.tsv \
             contrasts=$PWD/contrasts.tsv \
             max_intron_len=5000 \
    --directory output
```
