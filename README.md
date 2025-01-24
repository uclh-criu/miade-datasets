# miade-datasets

This repository contains the scripts and data used to generate the concept lists for MiADE's MedCAT models.


The files in each folder may include:

- README.md
- R or python scripts
- Log files
- CSV files containing the data

Each folder in `model_files` contains the CDB data and lookup tables for a deployment of MiADE.

Folders:

```
.
├── README.md
├── model_files
│   └── uclh_trial_2023
├── sample_notes
├── scripts
├── synthetic_data
└── vocabularies
    └── wikipedia_vocab
```

- sample_notes - fake clinical notes for testing
- Vocabularies - sets of word embeddings to build medcat models with
  - wikipedia_vocab - an open vocabluary which has been trained on a wikipeda dump, and contains no other data.
- model_files - contains directories which contain CDB data and lookup tables for a MiADE deployment
  - uclh_trial_2023 - the lookup tables and CDB data for the MiADE version 1 and 2 trial deployments at UCLH.
- scripts - scripts used to generate files in cdb_and_model_files_sep_2023
