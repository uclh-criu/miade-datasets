# This script updates the Sep 2023 MedCAT problems CDB to include
# case-sensitive synonyms for acronyms to ensure they are not
# misinterpreted

library(data.table)
library(Rdiagnosislist) # using 7 Oct 2024 version on GitHub
# Load May 2022 SNOMED
SNOMED <- readRDS('~/Terminologies/SNOMED_May2022.RDS')
data('MANUAL_SYNONYMS')
NEWCDB <- createCDB(WN = downloadWordnet(),
	MANUAL_SYNONYMS = MANUAL_SYNONYMS)
# NEWCDB <- createCDB(WN = WN, TRANSITIVE = CDB$TRANSITIVE, MANUAL_SYNONYMS = MANUAL_SYNONYMS)

# Find case-sensitive findings - if entirely upper case
# then exchange this in the CDB
NEW <- NEWCDB$FINDINGS[term == toupper(term),
	.(cui = conceptId, name = tolower(gsub('^ | $', '', term)),
	new_name = gsub('^ | $', '', term))]
NEW <- NEW[!duplicated(NEW)]

# Load existing CDB and de-duplicate
EXISTING <- fread('cdb_and_model_files_sep_2023/cdb_problems.csv')
EXISTING[, name := gsub(' +', ' ', name)]
EXISTING <- EXISTING[!duplicated(EXISTING)]

# Change capitalised names to all capitals in CDB (e.g.
# AF, COPD, TIA etc.)
EXISTING[, new_name := NEW[EXISTING, on = c('cui', 'name')]$new_name]
EXISTING[!is.na(new_name), name := new_name]
EXISTING[, new_name := NULL]

# Set name_status to N for 2 and 3 character concepts (this
# means they will always be disambiguated)
EXISTING[nchar(name) < 4, name_status := 'N']
fwrite(EXISTING, file = 'cdb_and_model_files_sep_2023/cdb_problems.csv')
