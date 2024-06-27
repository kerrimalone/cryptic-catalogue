# Analysis Overview

### Generation of training and validation datasets  
Please see [this html document](./creating_training_validation_sets/drprg-data-summary.html) for a detailed overview of how the training and validation datasets were constructed. A brief top level summary follows below.

#### Training
Training data: 25033 samples, file = [training_data_20231122.tsv](/./creating_training_validation_sets/training_data_20231122.tsv)  
Analysis file:  [creating_training_validation_datasets.Rmd](./creating_training_validation_sets/creating_training_validation_datasets.Rmd)

The training data was generated from CRyPTIC data ([CRyPTIC_reuse_table_20231208.csv](https://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/reuse/CRyPTIC_reuse_table_20231208.csv)) and mykrobe data ([mykrobe.20231121.tsv](./creating_training_validation_sets/mykrobe.20231121.tsv)). 501 samples removed due to duplication between the sets.  

131 samples removed from the resulting training data as they are part of the validation set.


#### Validation  
Validation data: 8914 samples, file = validation_set_20231110.pass.tsv


The validation set came from a curated a dataset of ~45k samples that includes samples used to train the WHO catalogue (n = 35.5k samples) and extra samples from publications/ENA (n = 9.5k samples). These were collected for the DrPRG publication by Hall et al. https://github.com/mbhall88/drprg.  

The samples have been phenotyped in various ways (MGIT, plates, proportional agar tests) and binary R/S phenotypes are reported.
As our goal is to compare the herein newly created CRyPTIC catalogue to the WHO catalogue, we excluded samples used to train the WHO catalogue from our validation set. Therefore, our validation dataset is represented by the 9.5k samples scraped from publications/ENA.
The unique identifiers for both the training and validation datasets are ENA accession numbers.  
 
 
### Catalogue construction
#### Variant calling
Test and validation data was downloaded from the ENA using ENA accessions and processed using Clockwork `regenotype`. FRS == 0.9 was chosen by default.

    nextflow run ../minos/nextflow/regenotype.nf \
        -resume \
        -with-singularity /nfs/research/zi/mhunt/Containers/minos_v0.12.5.img \
        -work-dir /hps/nobackup/iqbal/mhunt/tmp/regeno \
        -ansi-log false \
        -c ../minos/nextflow/regenotype.config \
        -with-trace nf.trace.txt \
        -profile large \
        --ref_fasta /nfs/research/zi/mhunt/Cryptic2/Refs/Ref.H37Rv/ref.fa \
        --manifest manifest.tsv \
        --outdir  Out

### Variant table preparation
Gnomonicus was used to generate the input mutations table for catalogue construction (Martin Hunt, output here:/hps/nobackup/iqbal/mhunt/Gnomicus/Out/).  
The VCFs were filtered beforehand to reduce RAM and run time to ~3GB per sample, 25min runtime. A catalogue was not supplied.  
The resulting x.variants.csv and x.mutations.csv were processed afterwards whereby
the FRS data from the x.variants.csv file was appended to the x.mutations.csv file for each sample. See `make_cat_builder_input.py` to do this.


### Catalogue generation
Catomatic was used to generate per drug catalogues. See `main.py` for this pipeline. Binomial testing was chosen as the statistical test. An experiment was run to choose the best background and CI values for each drug catalogue construction using a grid search with custom scoring that prioritises sensitivity, specificity and coverage* in that order.  

  *U predictions for samples that have phenotype R or S.


WHO catalogue

Training set evaluation

Test set evaluation
with gnomonicus



# Results  
See `results` folder for the main data files, results, graphs and tables.

1. Training data: 25033 samples, file = [training_data_20231122.tsv](creating_training_validation_sets/training_data_20231122.tsv)  
2. Validation data: 8914 samples, file = [validation_set_20231110.pass.tsv](creating_training_validation_sets/validation_set_20231110.pass.tsv)  

3. The ![training set phenotypes](creating_training_validation_sets/./creating_training_validation_sets/):   




### Training and validation results
I’ve been exploring the data and I attach a table showing the sources (n = 19) of these samples. You can see in the table that the majority of samples were phenotyped with MGIT (where phenotype method is stated).
 
I also attach a graph showing the breakdown of phenotypes for n = 22 drugs. Note, there is at least one drug phenotype for each sample in this set.


### Test set evaluation results
• Where variations in the summary metrics occur, they tend to be minimal.
• We out-perform WHOv1 for capreomycin and kanamycin predictions.
• We demonstrate better sensitivity but worse specificity for moxifloxacin and levofloxacin predictions. A quick glance places suspicions on different U and S classifications between the catalogues. 
• We demonstrate better specificity but worse sensitivity for ethambutol predictions and to a lesser extent for isoniazid and rifampicin predictions.
• We are slightly less sensitive for streptomycin predictions.
