This repository contains the workflow and config files for the Dr. PRG paper.
See https://github.com/mbhall88/drprg for the Dr. PRG software repository.

To run the pipeline you will need [`conda`](https://docs.conda.io/en/latest/) and
[`singularity`](https://github.com/apptainer/singularity) installed. You can find a
conda environment specification file at `environment.yaml`.

## Samplesheets

The CSV files containing the accessions and their associated phenotypes
are [`config/illumina.samplesheet.csv`](./config/illumina.samplesheet.csv)
and [`config/nanopore.samplesheet.csv`](./config/nanopore.samplesheet.csv). The
versions of these files with the `.pass.csv` extension contain only those files that
passed the quality control step we outline in the paper.

The `config/samplesheets/` directory contains the individual samplesheets that were used
to construct these main samplesheets.
See [`workflow/notebook/notepad.ipynb`](./workflow/notebook/notepad.ipynb) for a
notebook that was used to clean and combine all of these data sources, along with the
links to the publications that generated them.