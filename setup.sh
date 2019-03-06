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

# download singularity image for R analysis
singularity pull library://wkretzsch/default/kannan2019:sha256.b23184cec3ce85a98f6d9cc56c8c1ceb2787216301a31b1686debc3b3755a456 kannan2019.sif
