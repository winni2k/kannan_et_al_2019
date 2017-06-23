#!/usr/bin/env bash
set -euo pipefail
set +u
. ./activate.sh
set -u

conda install -y --file requirements.txt -c conda-forge -c bioconda -c anaconda || \
    conda install -y --file requirements.txt -c conda-forge -c bioconda -c anaconda

pip install git+git://github.com/kvg/dmpy.git

. activate rna_seq_py27
pip install --requirement requirements-py27.txt
. activate rna_seq


# link data
pushd data
for path in /data3/wkkg/pavitra/*
do
    test -L $(basename $path) || -ln -s $path
done


