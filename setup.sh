#!/usr/bin/env bash
set -euo pipefail

conda env create -vv -f environment.yml
conda env create -vv -f environment_py27.yml
conda activate rna_seq

# manually download gatk and register it
mkdir -p opt
if [ ! -f opt/*/GenomeAnalysisTK.jar ]; then
    curl 'https://software.broadinstitute.org/gatk/download/auth?package=GATK' | tar -xjf - -C opt
fi
gatk-register opt/*/GenomeAnalysisTK.jar

# link data
# pushd data
# for path in /data3/wkkg/pavitra/*
# do
#     test -L $(basename $path) || -ln -s $path
# done


