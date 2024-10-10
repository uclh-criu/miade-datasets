# This script updates the Sep 2023 MedCAT problems CDB to exclude
# blacklisted findings and also findings or events that are not
# in the Epic EDG file.

library(data.table)
library(Rdiagnosislist) # using 9 Oct 2024 version on GitHub
# Load May 2022 SNOMED
SNOMED <- readRDS('~/Terminologies/SNOMED_May2022.RDS')

# Load existing CDB and blacklist
EXISTING <- fread('cdb_and_model_files_sep_2023/cdb_problems.csv')
existing <- unique(as.SNOMEDconcept(EXISTING$cui))
blacklist <- as.SNOMEDconcept(scan('cdb_and_model_files_sep_2023/lookups/problem_blacklist.csv',
	'character'))

# Load list of SNOMED concepts in EDG
EDG <- fread('~/Terminologies/SNOMED Subset Update_v10_31.3.0_20210120000001 UK.csv')
edg <- as.SNOMEDconcept(EDG$conceptId)

findings <- descendants('Clinical finding')
disorders <- descendants('Disorder')
events <- descendants('Event')

# Exclude events and blacklisted findings
exclude <- intersect(existing, union(intersect(findings, blacklist), events))
new <- setdiff(existing, exclude)

# Find out if any findings in new are not in edg
findings_not_in_edg <- setdiff(intersect(new, findings), edg)

# Disorders not in edg (a subset of findings)
disorders_not_in_edg <- setdiff(intersect(new, disorders), edg)

# Removing findings not in edg from cdb
new <- setdiff(new, findings_not_in_edg)

# Write out
fwrite(EXISTING[cui %in% new], file = 'cdb_and_model_files_sep_2023/cdb_problems.csv')
