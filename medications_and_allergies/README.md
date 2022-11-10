# Medications and allergies

This folder contains the scripts and data used to generate the list of concepts used for the medications and allergies extractor in MiADE. You will also need to download a copy of SNOMED CT. These files were built using the May 2022 UK SNOMED CT release.

## Concept database for MedCAT

- med_allergy_reaction_MedCAT_cdb.csv = concept database for MedCAT
    - cui = SNOMED CT conceptID
    - name = term description
    - ontologies = all 'SNO' (SNOMED CT)
    - name_status = 'P' (preferred - one per concept) or 'A' (alternative)

## Lookups for post processing

- MiADE_med_allergy_reaction_lookup.csv = concepts in med_allergy_reaction_MedCAT_cdb including category, dose, units
- MiADE_VTM_VMP_dose_lookup.csv = for mapping VTM to VMP with the appropriate dose, for simple tablets only
- MiADE_ingredients_lookup.csv = ingredients (substances) included in each product
- MiADE_causative_agent_lookup.csv = for assigning a substance for the allergy
- MiADE_all_substances_lookup.csv = includes substances which are also products

## Script

- Generate_MiADE_med_allergy_cdb_v4.R = R script

## Additional files for reference

- allergies_training_sample.csv = sample CSV file in a format for training the meta-annotation algorithms

- Symptom_codelists.csv = Reference list of symptoms as Read codes from https://www.nature.com/articles/s41591-022-01909-w

- allergy_sample.csv = custom-created CSV file for converting to MedCATtrainer JSON

    Columns:

    - text = contains text with mapped phrase encased in opening "{X " and closing "}"
       - (e.g. he had a {p heart attack}).
       - {p ...} = problem
       - {m ...} = medication or allergy substance
       - {r ...} = reaction
    - p = problem SNOMED CT conceptId|description
    - m = medication SNOMED CT conceptId|description
    - r = reaction SNOMED CT conceptId|description
    - p_meta_relevance = Present | Historic | Irrelevant
    - p_meta_confirmed = Confirmed | Suspected | Negated
    - p_meta_laterality = No laterality | Left | Right | Bilateral
    - m_meta_category = Adverse reaction | Taking | Irrelevant
    - m_meta_allergytype = Allergy | Intolerance | Unspecified
    - m_meta_severity = Mild | Moderate | Severe | Unspecified 
    - r_meta_reactionpos = Not a reaction | Before | After
