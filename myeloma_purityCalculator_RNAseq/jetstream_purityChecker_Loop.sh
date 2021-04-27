#!/usr/bin/env bash

# Assumes the script is initiated in a study output directory of jetstream pipeline

for project in `find . -maxdepth 1 -type d ! -name '.*' -printf '%f\n'`
do
  echo
  echo "-------------------------------------"
  echo "Processing project: $project"
  # navigate to RNA alignments folders
  cd $project/rna/alignment/star
  for sample in `find . -maxdepth 1 -type d ! -name '.*' -printf '%f\n'`
  do
    echo
    echo "Calculating results for: $sample"
    cd sample
    # Call purity checker script
    /home/jkeats/git_repositories/MMRF_CoMMpass/myeloma_purityCalculator_RNAseq/purityChecker_GRCh38_e98.sh \
      -t 10 \
      -b ${sample}.star.cram \
      -a CRAM
    # move back to sample directory
    cd ..
  done
  # move back to project directory
  cd ../../../..
done