# Genomic Copy Number Signatures Based Classifiers for Subtype Identification in Cancer

This repo hosts the scripts used in the study of [Signatures of Discriminative Copy Number Aberrations in 31 Cancer Subtypes](https://info.baudisgroup.org/publications/2020-12-18-publication-Bo-classifiers/).

## Data 
The open-access data from Progenetix and TCGA, and restricted data from PCAWG were used in the study.

 - The complete lists of the samples used in each data repository are provided in ```/data```.
 - The open-access data used in the study is available at [Progenetix](https://progenetix.org/gao-2021-signatures/search/).
 - In accordance with the data access policies of the ICGC, researchers need to apply to the [ICGC Data Access Compliance Office](http://icgc.org/daco) for PCWAG data access. Instruction on accessing restricted data from the ICGC/PCAWG is available at https://docs.icgc.org/pcawg/data/.


## File structure

```
alt-pipeline/			Scripts to generate .pkl files from the provided sample data at Progenetix.
classification/			Scripts of the classification experiments.
data/				The external & generated data used during the study.
integration/			Scripts to process the original data from Progenetix, TCGA, and PCAWG.
plots/				Scritps of all figures.
signatures/			Scripts for feature & signature generation using Autoencoder and LRP
```

## Workflow

### Data integration 
Copy number data from Progenetix, TCGA and PCAWG were preprocessed with the following steps, respectively:

- probe or segment data were lifted to hg38, if the original data was not in hg38.
- transformed to a uniform data structure and stored in mongodb.
- normalized using [mecan4CNA](https://github.com/baudisgroup/mecan4cna).

All data were combined in a single collection in mongodb (db:Rebased, collection:mecaned).

#### An example of the mongodb data structure 
```json
{
    "source" : "TCGA",
    "project" : "TCGA-BRCA",
    "sample_id" : "ae96c429-b221-4894-a45a-6aa4e8d32c71",
    "morphology" : "8500/3",
    "topography" : "Breast, NOS",
    "segments" : [
            {
            "chro" : "1",
            "start" : 3301765,
            "end" : 53333626,
            "probes" : 26594,
            "value" : -0.2022
        }
    ],
    "base" : 1.85,
    "level_distance" : 0.35,
    "normalized" : [
            {
            "chro" : "1",
            "start" : 3301765,
            "end" : 53333626,
            "probes" : 26594,
            "value" : -0.2504
        }
    ],
    "cytobands" : [
            {
            "start" : 0,
            "end" : 2300000,
            "name" : "p36.33",
            "note" : "gneg",
            "total_dup" : 0,
            "total_del" : 0,
            "dup_length" : 0,
            "del_length" : 0,
            "dup_count" : 0,
            "del_count" : 0,
            "chro" : "1",
            "ave_dup" : 0,
            "ave_del" : 0
        }
    ]
}
```
Here, `segments` stores the original data, `normalized` stores the normalized segments, and `cytobands` is used in the feature extraction procedure to store the summary of each cytoband. `base` and `level_distance` are parameters computed by mecan4CNA and are used for the normalization.

+

### Intermediate files
To facilitate the downstream pipelines, the following pickle files were created from the mongodb.

- `all_bands_meta.pkl`: the metadata of all samples.
- `all_bands.pkl`: the band features (weighted CNV average) of each sample.
- `all_bands_label.pkl`: the morphology, topography and organ labels of each sample.
- `all_bands_disease_label.pkl`: the morphology label of each sample.
- `all_bands_source_label.pkl`: the source of each sample.

### The alt-pipeline
When using the download data from Progenetix, please use the `alt-pipeline` instead of `integration` to preprocess data. Because of the difference in data, the `alt-pipeline` is not identical to the original pipeline.

Please run scripts in the following order:

1. load_data
2. calibration
3. combine_data
4. normalization
5. cytoband_data
6. gen_pickles

## Feature & signature generation
The procedure:

1. Build an autoencoder model using cytoband features
2. Extract high-weighting cytoband features
3. Build an autoencoder model using gene features (generated from high-weighting cytoband features)
4. Extract high-weight gene features
5. Generate signatures for cancer subtypes

## Classification
The procedure:

1. Filter data (with subtype signatures, enough samples)
2. Upsampling & downsampling during cross-validations
3. Multi-class classification of cancer subtypes
4. Extend classification results to organs of origin
