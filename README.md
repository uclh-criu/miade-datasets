# miade-datasets

This repository contains the scripts and data used to generate the concept lists for MiADE's MedCAT models.

Each folder contains the most up to date files in the root folder, and previous versions used for previous MedCAT models (along with scripts to generate them) in dated subfolders.

The files in each folder may include:

- README.md
- R or python scripts
- Log files
- CSV files containing the data

Folders:

- problems - datasets for MiADE problem list algorithm
- medication_and_allergies - datasets for MiADE medications and allergies algorithms
- sample_notes - fake clinical notes for testing
- Vocabularies - sets of word embeddings to build medcat models with
  - wikipedia_vocab - an open vocabluary which has been trained on a wikipeda dump, and contains no other data.
  
Folders for lookups and datasets currently in use:
- cdb_and_model_files_sep_2023 - lookups used for MiADE Version 1 (for problems/diagnoses only, implemented in UCLH Feb 2024) and Version 2 (problems/diagnoses only, updated blacklist for corrections, to be implemented in UCLH Apr 2024). Based on UK SNOMED CT May 2022 version.
- scripts - scripts used to generate files in cdb_and_model_files_sep_2023
