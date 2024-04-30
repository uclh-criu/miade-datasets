library(Rdiagnosislist)
library(data.table)
SNOMED <- readRDS('~/Terminologies/SNOMED_May2022.RDS')

# Load CDB
cdb_conceptIds <- unique(SNOMEDconcept(fread(
	'./cdb_and_model_files_sep_2023/cdb_problems.csv')$cui))

# Use SNOMED relationships table to find all ancestor / descendants pairs
SCT_isa <- SNOMEDconcept('Is a')

WORKING <- SNOMED$RELATIONSHIP[(sourceId %in% cdb_conceptIds |
	destinationId %in% cdb_conceptIds) & typeId == SCT_isa,
	.(childId = sourceId, parentId = destinationId)]

new_nrows <- nrow(WORKING)
old_nrows <- 0
while(new_nrows > old_nrows){
	WORKING <- rbind(WORKING, merge(
		WORKING[, .(childId, selfId = parentId)],
		WORKING[, .(selfId = childId, parentId)], by = 'selfId',
		allow.cartesian = TRUE)[,
		.(childId, parentId)])
	WORKING <- WORKING[!duplicated(WORKING)]
	old_nrows <- new_nrows
	new_nrows <- nrow(WORKING)
}

TRANSITIVE <- WORKING[childId %in% cdb_conceptIds &
	parentId %in% cdb_conceptIds,
	.(ancestorId = parentId, descendantId = childId)]

# Export transitive table
fwrite(TRANSITIVE,
	file = './cdb_and_model_files_sep_2023/lookups/transitive.csv',
	col.names = TRUE)
