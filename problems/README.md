# Problems

This folder contains the scripts and data used to generate the list of concepts used for the problem extractor in MiADE. You will also need to download a copy of SNOMED CT. These files were built using the Feb 2022 UK SNOMED CT release.

## Concept database for MedCAT

- full_condition_list_MedCAT_cdb.csv = concept database for MedCAT
    - cui = SNOMED CT conceptID
    - name = term description
    - ontologies = all 'SNO' (SNOMED CT)
    - name_statue = 'P' (preferred - one per concept) or 'A' (alternative)

## Lookups for post processing

Epic problem list entries do not include metadata attributes such as 'suspected' or 'negated'; they can only contain pre-coordinated SNOMED CT concepts. These lookup tables of base concepts and linked precoordinated negated, historic or suspected versions which can be provided as the MiADE output to Epic NoteReader.

- historic.csv = 'History of XXX' terms where available (otherwise they should be left as is)

- suspected.csv = 'Suspected XXX' terms (otherwise they should be ignored)

- negated.csv = 'Absence of XXX' terms  (otherwise they should be ignored)

- problem_blacklist.txt = List of SNOMED CT concept IDs (one per row) that should not be returned to NoteReader for presentation to the user. This consists of concepts that are too vague or generic to be useful.

## Script

Generate_MiADE_SNOMED_problems_cdb_v5.R = R script

## Additional files for reference

The 'historic_check' etc. documents are for manual checking, which include ambiguous SNOMED CT concepts (some of which might need to be fed back to SNOMED CT); these are not included in the 'historic', 'suspected' or 'negated' files for use by the algorithm.

- exclude_concepts.txt = suspected or historical concepts that are not to be used for annotation by MedCAT (MedCAT should use meta annotations with the base concept instead). 

- extra_condition_synonyms.csv = extra term descriptions that can be added to the vocab, derived from an old terminology system called OXMIS and a manually curated list of synonyms in the Freetext Matching Algorithm (https://github.com/anoopshah/freetext-matching-algorithm-lookups). Also includes alternative word orders for concepts with 2-3 words.

- historic_check.csv = 'History of XXX' terms where available

- suspected_check.csv = 'Suspected XXX' terms

- negated_check.csv = 'Absence of XXX' terms

- acronyms_to_check.csv = List of acronyms for manual check
