kannan_et_al_2019
=================

RNA-seq analyses for [kannan2019]_.

.. [kannan2019] Kannan, et al. (2019) Radiation resistant cancer cells enhance the survival and resistance of sensitive cells in prostate spheroids. bioRxiv: 564724. doi: https://doi.org/10.1101/564724

Setup
-----

- Install `conda <https://conda.io/en/latest/miniconda.html>`_ and activate an environment.
- Run ``bash setup.sh``

Run analysis
------------

Please obtain a copy of the raw sequencing data and store it in the ``data`` directory (to be updated).
::

   # activate the rna_seq conda evironment
   conda activate rna_seq

   # download dependencies
   python download.py

   # run sequence analysis
   python analysis.py -r

   # run differential expression analysis
   # also include "--bind /srv,/data,.:/scratch" if running on saga
   singularity exec --bind .:/scratch  kannan2019.sif bash -c 'cd /scratch && Rscript --vanilla  analysis.R'


Building the singularity container
----------------------------------

::

   # tested with singularity version 2.4
   sudo singularity build bioc3.8.simg Singularity.Bioc3.8
   sudo singularity build kannan2019.simg Singularity

   # for singularity version 3, convert to SIF
   singularity build kannan2019.sif knnan2019.simg

   # cleanup
   rm bioc.simg

