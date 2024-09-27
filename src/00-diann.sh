#!/bin/bash
#SBATCH -J diann_hdrg
#SBATCH -N 1
#SBATCH --cpus-per-task=96

module purge

DATA=`\path\to\data\directory`
OUT=`\path\to\output\directory`

./bin/diann-1.8.1 --dir $DATA \
        --gen-spec-lib --fasta-search --use-quant \
        --predictor \
        --fasta hdrg/UP000005640_validated.fasta \
        --threads 96 --verbose 1 --qvalue 0.01 --matrices \
        --min-corr 2.0 --corr-diff 1.0 --time-corr-only \
        --min-fr-mz 100 --max-fr-mz 1700 \
        --out $OUT/report.tsv \
        --out-lib $OUT/report-lib.tsv \
        --cut K*,R* --missed-cleavages 2 --min-pep-len 5 --max-pep-len 52 \
        --min-pr-mz 350 --max-pr-mz 1200 --min-pr-charge 2 --max-pr-charge 4 \
        --var-mods 3 --var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n \
        --monitor-mod UniMod:1 --rt-profiling --reanalyse --met-excision \
        --pg-level 1 --no-ifs-removal