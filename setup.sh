#!/usr/bin/env bash
set -euo pipefail

conda env create -f environment.export.yml
conda env create -f environment_py27.yml
set +u
. activate rna_seq
set -u

# manually download gatk and register it
mkdir -p opt
if [ ! -f opt/*/GenomeAnalysisTK.jar ]; then
    curl 'https://software.broadinstitute.org/gatk/download/auth?package=GATK' | tar -xjf - -C opt
fi
gatk3-register opt/*/GenomeAnalysisTK.jar
