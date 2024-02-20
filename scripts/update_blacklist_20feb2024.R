library(Rdiagnosislist)
library(data.table)

finding_observation <- SNOMEDconcept(unique(
	SNOMED$DESCRIPTION[active == TRUE, .(findingObservation =
	all(term %like% 'Finding|finding|Observation|observation')),
	by = conceptId][findingObservation == TRUE]$conceptId))

pt_condition <- setdiff(c(
	descendants('Body position finding', include_self = TRUE),
	descendants('General body state finding', include_self = TRUE),
	descendants('Patient condition finding', include_self = TRUE),
	descendants('Finding of temperature of skin', include_self = TRUE),
	SNOMEDconcept(c('Not pregnant', "Patient's condition stable"))),
	descendants('Disorder'))

specific_errors <- SNOMEDconcept(c(
	'SAM',
	'Genital finding',
	'Childhood granulomatous periorificial dermatitis',
	'Asymptomatic',
	'Bends',
	'Exercise state',
	'Acute ulcerative gingivitis',
	'No abnormality detected',
	'Internal malleolar torsion',         
	'Siti'))

dont_exclude <- descendants(c(
	'Amputee',
	'Somatic dysfunction',
	'Functional disease present',
	'Urobilinogenemia',
	'Pregnant',
	'Falls',
	'Recurrent falls', 
	'Elderly fall',
	'Falls infrequently',
	'Unexplained recurrent falls',
	'Unexplained falls',
	'Falls caused by medication'), include_self = TRUE)

# Load existing blacklist
existing_blacklist <- SNOMEDconcept(
	scan('./cdb_and_model_files_sep_2023/lookups/problem_blacklist.csv',
	what = 'character'))

# Load CDB
existing_cdb <- SNOMEDconcept(fread(
	'./cdb_and_model_files_sep_2023/cdb_problems.csv')$cui)

# Set new blacklist
new_blacklist <- sort(union(existing_blacklist,
	setdiff(intersect(existing_cdb,
		c(finding_observation, pt_condition, specific_errors)),
	dont_exclude)))

# Export new blacklist
fwrite(data.table(x = new_blacklist),
	file = './cdb_and_model_files_sep_2023/lookups/problem_blacklist.csv',
	col.names = FALSE)
