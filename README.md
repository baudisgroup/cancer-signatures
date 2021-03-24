## Genomic Copy Number Signatures Based Classifiers for Subtype Identification in Cancer

This repo hosts the scripts used in the study of [Signatures of Discriminative Copy Number Aberrations in 31 Cancer Subtypes](https://info.baudisgroup.org/publications/2020-12-18-publication-Bo-classifiers/).


The open-access data from Progenetix and TCGA, and restricted data from PCAWG were used in the study.

 - The complete lists of the samples used in each data repository are provided in ```/data```.
 - The open-access data used in the study is available at [Progenetix](https://progenetix.org/gao-2021-signatures/search/).
 - In accordance with the data access policies of the ICGC, researchers need to apply to the [ICGC Data Access Compliance Office](http://icgc.org/daco) for PCWAG data access. Instruction on accessing restricted data from the ICGC/PCAWG is available at https://docs.icgc.org/pcawg/data/.


### Data integration 
Copy number data from Progenetix, TCGA and PCAWG were preprocessed with the following steps, respectively:
- probe or segment data were lifted to hg38, if the original data was not in hg38.
- transformed to a uniform data structure and stored in mongodb.
- normalized using [mecan4CNA](https://github.com/baudisgroup/mecan4cna).

After preprocessing, all data were combined in a single collection in mongodb (db:Rebased, collection:mecaned).

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
To facilate the later pipelines, the following pickle files were created from the mongodb.

- cytoband summaies of each sample were computed with `signatures/cytoband_data.ipynb` and stored in the mongodb. Only samples will valid cytoband data were used in the following files.
- `all_bands_meta.pkl`: the meta data of all samples.
- `all_bands.pkl`: the band features (weighted CNV average) of each sample.
- `all_bands_label.pkl`: the morphology, topography and organ labels of each sample.
- `all_bands_disease_label.pkl`: the morphology label of each sample.
- `all_bands_source_label.pkl`: the source of each sample.


