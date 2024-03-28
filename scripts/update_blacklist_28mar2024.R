library(Rdiagnosislist)
library(data.table)
SNOMED <- readRDS('~/Terminologies/SNOMED_May2022.RDS')

# Update blacklist to exclude all findings - for MiADE Version 2
# Also exclude specific disorders which currently cause errors
# because of training errors or ambiguous acronyms

findings <- setdiff(descendants('Clinical finding', include_self = TRUE),
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
	'Siti',
	'Ankle fracture', 
	'Myofascial pain syndrome',
	'Common acute lymphoblastic leukemia',
	'Secondary malignant neoplastic disease',
	'Primary open angle glaucoma'))

# Note that ankle fracture is an error because of training error
# to be corrected in future version. Others are due to
# disambiguation problems or acronyms not recognised - to
# address with additional training. coag = chronic open angle glaucoma.
# ankle = ankle fracture. siti = person's name. 'secondaries'
# aug (august / acute ulcerative gingivitis),
# imt (internal malleolar torsion).
# Colour findings should be excluded (descendants of 'Color finding')
# Ankle fracture is very common so this should be addressed in the
# training data as soon as possible.

# Only keep findings which are potentially important to include
# on the problem list to avoid concerning findings being
# forgotten e.g. potential cancers
dont_exclude <- descendants(c(
	'Amputee',
	'Mass of body region',
	'Functional disease present',
	'Seizure',
	'Syncope',
	'Falls',
	'Recurrent falls', 
	'Elderly fall',
	'Unexplained recurrent falls',
	'Unexplained falls',
	'Falls caused by medication'), include_self = TRUE)

# Load existing blacklist
existing_blacklist <- SNOMEDconcept(
	scan('./cdb_and_model_files_sep_2023/lookups/problem_blacklist.csv',
	what = 'character'))

# Save as 'version 1'
fwrite(data.table(x = sort(unique(existing_blacklist))),
	file = './cdb_and_model_files_sep_2023/lookups/problem_blacklist_V1.csv',
	col.names = FALSE)

# Load CDB
existing_cdb <- SNOMEDconcept(fread(
	'./cdb_and_model_files_sep_2023/cdb_problems.csv')$cui)

# Set new blacklist
new_blacklist <- union(existing_blacklist,
	setdiff(intersect(existing_cdb,
		c(findings, specific_errors)),
	dont_exclude))

# Export new blacklist
fwrite(data.table(x = sort(unique(new_blacklist))),
	file = './cdb_and_model_files_sep_2023/lookups/problem_blacklist.csv',
	col.names = FALSE)

