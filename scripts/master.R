# Master script for generating MiADE lookup tables
setwd('/home/anoop/pCloudDrive/anoop_work/2019to2021/MiADE/Terminologies/scripts_jul2023/')

# Specify file locations:

# REF = filepath for reference files
REF <- '~/pCloudDrive/anoop_work/2019to2021/MiADE/Terminologies/REF/'

# WORKING = filepath for intermediate files for manual checking
WORKING <- '~/pCloudDrive/anoop_work/2019to2021/MiADE/Terminologies/WORKING/'

# TEMPLATES = filepath for template files for synthetic metacat generation
TEMPLATES <- '~/pCloudDrive/anoop_work/2019to2021/MiADE/Terminologies/TEMPLATES/'

# OUTPUT = output files for use by MedCAT / MiADE
OUTPUT <- '~/pCloudDrive/anoop_work/2019to2021/MiADE/Terminologies/OUTPUT/'

# Load SNOMED dictionary
SNOMED <- readRDS('~/Terminologies/SNOMED_May2022.RDS')
READMAPS <- readRDS('~/Terminologies/READMAPS.RDS')

## Generate extra synonyms - problems

source('gen_synonyms_problems_v3.R')

## Generate CDB - problems

source('gen_cdb_problems_v3.R')

#~ - cdb_problems
#~ - suspected.csv: findingId, situationId
#~ - historic.csv: findingId, situationId
#~ - negated.csv: findingId, situationId
#~ - problem_blacklist.csv: [No headings] conceptId for concepts to exclude
#~ - REF/extra_condition_synonyms

## Generate CDB - meds and allergies

source('gen_cdb_medallerg_v3.R')

#~ - cdb_medallerg
#~ - valid_meds: [No headings] conceptId for valid medication
#~ - reactions_subset: reactionsId, subsetId (for subset of reactions that can be mapped to Epic)
#~ - allergens_subset: allergensId, subsetId, allergenType (for substances that can be mapped to Epic. allergenType is 'drug class', 'drug ingredient' etc.)
#~ - allergy_type: allergenType, adverseReactionType, adverseReactionId (adverseReactionType is 'intolerance' or 'allergy' or 'adverse reaction', and the ID is the SNOMED concept ID for the root concept)
#~ - vtm_to_vmp: vtmId, dose, unit, vmp
#~ - vtm_to_text: vtmId, fallbackText e.g. 'LANSOPRAZOLE ORAL' for VTMs not mapped. 

source('gen_cdb_medallerg_v3.R')

## Generate metacat

N <- 750
source('gen_metacat_v6.R')

############

#~ Folders:

#~ Lookup table specification

#~ https://www.notion.so/miade/Lookup-Table-Specifications-87037e3ec08240f4a73c912397af84cf

#~ REF
#~ - contains: SNOMED symptom list, SNOMED dictionary, SNOMED frequency distribution

#~ WORKING


#~ OUTPUT
#~ - contains CDB and lookups when finished
