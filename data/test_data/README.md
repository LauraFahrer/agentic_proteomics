# Test Data
The test datasets are obtained from quantms, a proteomics data processing ecosystem that provides standardized results from public mass spectrometry experiments. All datasets belong to the Human Cancer Cell Line Collection.

## Dataset structure
A dataset consists of several files. The following are particularly important to understand the experimental design and analyzing the data:

- *.sdrf.tsv: Sample metadata and experimental design, including sample characteristics, experimental conditions, and measurement file references
- *.feature.parquet: Mass spectrometry feature data
- *.run.parquet: Metadata describing individual measurement runs
- *.sample.parquet: Sample-level metadata
- *.pg.parquet: Protein-group information
- *_msstats_in.csv.gz: Quantitative proteomics data prepared for analysis with MSstats

## Initial Test Dataset

### Selection Criteria
To select the initial test dataset the following criterias have been applied:

- manageable number of samples
- sufficient biological replicates to estimate biological variability
- balanced group size where possible
- unambiguous labels
- metadata included
- sufficient number of identified proteins
- quantitative proteomic data included
- clear experimental design and interpretability

### Selection Process
- considered runs with sample size 15 to 50
-