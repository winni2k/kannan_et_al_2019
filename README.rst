kannan_et_al_2019
=================

RNA-seq analyses for [kannan2019]_.

.. [kannan2019] Kannan, et al. (2019) Radiation resistant cancer cells enhance the survival and resistance of sensitive cells in prostate spheroids  Cell interactions paper RNA-seq analyses

Setup
-----

- Install `conda <https://conda.io/en/latest/miniconda.html>`_ and activate an environment.
- Run ``bash setup.sh``

Run analysis
------------

Please obtain a copy of th raw sequencing data and store it in `data`.

```
# activate the rna_seq conda evironment
conda activate rna_seq

# download dependencies
python download.py

# run sequence analysis
python analysis.py -r

# run differential expression analysis
# also include "--bind /srv,/data,.:/scratch" if running on saga
singularity exec --bind .:/scratch  kannan2019.simg bash -c 'cd /scratch && Rscript --vanilla  analysis.R'
```

Building the singularity container
----------------------------------

```
# tested with singularity version 2.4
sudo singularity build bioc3.8.simg Singularity.Bioc3.8
sudo singularity build kannan2019.simg Singularity
rm bioc.simg
```
