# Script to create medication and allergy concept databases and set of reactions
# for MiADE 

# Anoop Shah

# Version 2: 11 Sep 2023
#
# Improvements:
# - Uses additional synonyms from products
# - Includes additional synonyms with space between number and dose units
# - Includes additional synonyms for substances which are the substance names in FDB

# Version 3: 14 Nov 2023
#
# Improvements:
# - Includes additional synonyms for 'modified release'
#
# Bug fixes:
# - Includes multiple steps in hierarcy for complete set of VTM to VMP conversion
# - Avoids multiple rows in VTM to VMP conversion table

library(data.table)
library(Rdiagnosislist)

######################################
# Import relevant files

# Import SNOMED CT dictionary (created using loadSNOMED function)
# Using May 2022 SNOMED dictionary and Read maps
# Import SNOMED CT dictionary (created using loadSNOMED function)
if (!exists('SNOMED')){
	SNOMED <- readRDS('~/Terminologies/SNOMED_May2022.RDS')
}
# Import SNOMED CT dictionary (created using loadSNOMED function)
if (!exists('READMAPS')){
	READMAPS <- readRDS('~/Terminologies/READMAPS.RDS')
}

defaultfilepath <- '~/pCloudDrive/anoop_work/2019to2021/MiADE/Terminologies/'
if (!exists('WORKING')){ WORKING <- paste0(defaultfilepath, 'WORKING/') }
if (!exists('OUTPUT')){ OUTPUT <- paste0(defaultfilepath, 'OUTPUT/') }
if (!exists('REF')){ REF <- paste0(defaultfilepath, 'REF/') }
if (!exists('LOGS')){ LOGS <- paste0(defaultfilepath, 'LOGS/') }

# Load reference files

# Read symptoms derived from Long Covid project
READ_SYMPTOMS <- fread(paste0(REF, 'Symptom_codelists.csv'))

# Additional concepts from manual review
MED_EXTRA <- NULL
MED_EXTRA <- fread(paste0(REF, 'med_extra.csv'))

# Epic allergen files
ELG <- fread(paste0(REF, 'ELG_ALLERGEN.csv'))
EXM <- fread(paste0(REF, 'substance_subset.csv'))

# List of custom oral medication records
ORAL <- fread(paste0(REF, 'Oralmeds_with_auto_SCT_map.csv'))

# Medications, substances and reactions to exclude (manual list)
MED_SUB_REACTION_EXCLUDE <- fread(paste0(REF, 'med_sub_reaction_exclude.csv'))

# Commence log file
sink(paste0(LOGS, 'gen_cdb_medallerg_v3_', Sys.time(), '.log'), split = TRUE)

######################################

VTM <- children('Virtual therapeutic moiety') # 2944
VTM <- description(VTM, include_synonyms = TRUE)[
	type  == 'Synonym', .(conceptId, term, category = 'vtm')]

VMP <- children('Virtual medicinal product') # 18838
VMP <- description(VMP, include_synonyms = TRUE)[
	type  == 'Synonym', .(conceptId, term, category = 'vmp')]

# Substances
SUBS <- setdiff(descendants('Substance'), children('Substance')) # 27782
SUBS <- description(SUBS, include_synonyms = TRUE)[
	type  == 'Synonym', .(conceptId, term, category = 'substance')]

# Add additional synonyms from products in the product hierarchy
# that are similar to substances. Note that the SNOMED CT information
# model does not include the equivalence as a defined relationship,
# hence we need to use text matching.
PROD <- descendants('Medicinal product') # 2944
PROD <- description(PROD, include_synonyms = TRUE)[
	type  == 'Synonym', .(productId = conceptId, term)]
COMMON <- merge(SUBS[, .(substanceId = conceptId, term)], PROD,
	by = 'term')
# COMMON has columns: productId, term, substanceId
COMMON <- COMMON[, .(productId, substanceId)]
COMMON <- COMMON[!duplicated(COMMON)]
# nrow(COMMON)
# uniqueN(COMMON$productId)
# uniqueN(COMMON$substanceId)
PROD <- merge(PROD, COMMON, by = 'productId')
SUBS <- rbind(SUBS, PROD[, .(conceptId = substanceId, term, category = 'substance')])

# Add ALLERGEN_NAME as extra synonym
SUBS <- rbind(SUBS, merge(ELG[, .(ALLERGEN_ID, ALLERGEN_NAME)],
	EXM, by = 'ALLERGEN_ID')[,
	.(conceptId = bit64::as.integer64(SNOMEDCT),
	term = tolower(ALLERGEN_NAME),
	category = 'substance')])
SUBS <- SUBS[!is.na(conceptId)]
SUBS <- SUBS[!duplicated(SUBS)]

# Check substance synonyms 
# SUBS[term == 'eggs']
# SUBS[term == 'NSAID']
# SUBS[term == 'nsaids']

# Create the initial joint medication / substances dataset
M <- rbind(SUBS, VTM, VMP)

#### ADD EXTRA SUBSTANCES FROM MANUAL REVIEW

# This is temporary - needs to be rewritten to accept concepts + terms,
# check for versioning and set the category according to semantic type.
# list of useful trade names and med synonyms to include
if (!is.null(MED_EXTRA)){
	message('Adding substances from manual review')
	# Add extra substances or medications from manual review
	M <- rbind(M, MED_EXTRA[, .(conceptId, term = name, category = 'substance')])
	M <- rbind(M, description(MED_EXTRA$conceptId, include_synonyms = TRUE)[
		type  == 'Synonym', .(conceptId, term, category = 'substance')])
}

# Fully specified name for (clinical drug) may be never used in text
# e.g. "Product containing precisely bendroflumethiazide 2.5 milligram/1 each
#       conventional release oral tablet (clinical drug)"
# Hence keep only synonyms for clinical drugs

# Remove space between number and units
M[, term := sub('([0-9]) (ml|nanogram|nanograms|g|mg|microgram|micrograms)',
	'\\1\\2', term)]

# Add synonyms with space between number and unit
EXTRA <- M[term %like% '([0-9])(ml|nanogram|nanograms|g|mg|microgram|micrograms)']
EXTRA[, term := sub('([0-9])(ml|nanogram|nanograms|g|mg|microgram|micrograms)',
	'\\1 \\2', term)]
M <- rbind(M, EXTRA)

# Add synonyms for modified-release = MR
M <- rbind(M, M[term %like% 'modified-release' & category == 'vmp',
	.(conceptId, term = sub('modified-release|prolonged-release', 'mr',
	term), category = 'vmp')])

# Add alternate order for modified-release / MR
M <- rbind(M, M[term %like% 'modified-release|mr' & category == 'vmp',
	.(conceptId, term = sub(paste0(' ([[:digit:]]+ ?[microgranu]+) ',
	'(modified-release|mr) (capsule|capsules|tablet|tablets)'),
	' \\2 \\1', term), category = 'vmp')])

# Add synonyms for oral tablets
M <- rbind(M, M[term %like% ' oral tablet' &
	category == 'vmp', .(conceptId,
	term = sub(paste0(' oral tablet'), ' tablet', term),
	category = 'vmp')])

# Add synonyms for oral capsules
M <- rbind(M, M[term %like% ' oral capsule' &
	category == 'vmp', .(conceptId,
	term = sub(' oral capsule', ' capsule', term),
	category = 'vmp')])

# Remove substances that are also VTMs
SUBS_VTM_OVERLAP <- merge(SUBS[, .(subId = conceptId, term)],
	VTM[, .(vtmId = conceptId, term)], by = 'term')
M <- M[!(conceptId %in% SUBS_VTM_OVERLAP$subId)]

# Remove duplicates
M <- M[!duplicated(M)]

# Mark as 'preferred' the longest term for each concept
M <- M[order(conceptId, -nchar(term))]
M[, preferred := c(TRUE, rep(FALSE, .N - 1)), by = conceptId]

# Extract dose and units if available and if only a single dose
M[, dosecategory := 'NOT_APPLICABLE']
r_text_with_at_least_2_numbers <- '^[^0-9]*[0-9]+\\.{0,1}[0-9]*[^0-9\\.]+[0-9]+\\.{0,1}[0-9]*.*$'
M[category == 'vmp', dosecategory := sub(r_text_with_at_least_2_numbers,
	'MULTIPLE_NUMBERS', term)]
M[category == 'vmp' & dosecategory != 'MULTIPLE_NUMBERS',
	dosecategory := ifelse(term %like% ' (oral |)(capsule|tablet)(s|)$',
	'SIMPLE_TABLET', 'COMPLEX')]
M[dosecategory == 'SIMPLE_TABLET', dose := 
	as.numeric(sub('^[^0-9]+ ([0-9\\.]+)[ ]*(ml|nanogram|nanograms|g|mg|microgram|micrograms) [^0-9]*$',
	'\\1', term))]
M[dosecategory == 'SIMPLE_TABLET' & !is.na(dose), units :=
	sub('^[^0-9]+ ([0-9\\.]+)[ ]*(ml|nanogram|nanograms|g|mg|microgram|micrograms) [^0-9]*$',
	'\\2', term)]

MDOSE <- M[!is.na(dose), .(conceptId, dose, units)]
MDOSE <- MDOSE[!duplicated(MDOSE)]
MDOSE[, Ndoses := .N, by = conceptId]
# Limit to entries with a single unambiguous dose
# omit this one:
#   conceptId                                               term category
#1: 322519002 Dextropropoxyphene hydrochloride 65mg oral capsule      vmp
#2: 322519002       Propoxyphene hydrochloride 65mg oral capsule      vmp
#3: 322519002                   Dextropropoxyphene 60mg capsules      vmp
MDOSE <- MDOSE[Ndoses == 1]
MDOSE[, Ndoses := NULL]

#### VTM TO VMP conversion tool

# Create a VTM to VMP conversion tool based on dose
# VTM is ancestor of VMP and also is_a VTM
# Extract descendants of VTM and filter by VMPs
# Modified 1 Nov 2023 to include all descendants (not just immediate
# children)

VTM_R <- SNOMED$RELATIONSHIP[active == TRUE &
	typeId == SNOMEDconcept('Is a') & destinationId %in% VTM$conceptId,
		.(conceptId = sourceId, vtmId = destinationId)]
cat('\nVTM to VMP conversion - one level from SNOMED CT\n')
print(nrow(VTM_R))
max_gen <- FALSE
while (!(max_gen)){
	n_concepts <- uniqueN(VTM_R$conceptId)
	VTM_R <- merge(VTM_R, SNOMED$RELATIONSHIP[active == TRUE &
		typeId == SNOMEDconcept('Is a'),
		.(vmpId = sourceId, conceptId = destinationId)],
		by = 'conceptId', all.x = TRUE)
	VTM_R <- rbind(VTM_R[, .(conceptId = vmpId, vtmId)],
		VTM_R[, .(conceptId, vtmId)])
	VTM_R <- VTM_R[!duplicated(VTM_R)]
	if (uniqueN(VTM_R$conceptId) == n_concepts){
		max_gen <- TRUE
	}
}
cat('\nVTM to VMP conversion - multiple levels from SNOMED CT\n')
print(nrow(VTM_R))

VTM_VMP <- merge(M[category == 'vmp', .(conceptId)], VTM_R)[
	vtmId %in% VTM$conceptId]
VTM_VMP <- VTM_VMP[!duplicated(VTM_VMP)]
VTM_VMP[, dose := MDOSE[VTM_VMP, on = 'conceptId']$dose]
VTM_VMP[, units := MDOSE[VTM_VMP, on = 'conceptId']$units]
VTM_VMP <- merge(VTM_VMP, description(VTM_VMP$conceptId)[, .(conceptId, term)])
VTM_VMP <- VTM_VMP[!duplicated(VTM_VMP)]

# Keep only one VMP per VTM/dose/unit combination
# Keep conventional release and tablets rather than capsules
VTM_VMP[, mr := term %like% 'MR|SR|modified|sustained']
VTM_VMP[, capsule := term %like% 'capsule']
VTM_VMP <- VTM_VMP[!is.na(dose) & !is.na(units)][
	order(vtmId, dose, units, mr, capsule, -nchar(term))]
VTM_VMP[, keep := c(TRUE, rep(FALSE, .N - 1)), by = .(vtmId, dose, units)]

cat('\nFinal VTM to VMP table\n')
print(table(VTM_VMP$keep))
cat(nrow(VTM_VMP), 'rows\n')
rm(VTM_R)

##################################
# Extract ingredient concept to allow link from substance to product

# relatedConcepts(product, typeId = 'Has specific active ingredient (attribute)')
# is a substance
# relatedConcepts(product, typeId = 'Has active ingredient (attribute)')
INGREDIENTS <- SNOMED$RELATIONSHIP[active == TRUE &
	typeId %in% SNOMEDconcept(c('Has specific active ingredient (attribute)',
	'Has active ingredient (attribute)')), .(conceptId = sourceId,
	ingredientId = destinationId)]

##################################
# Findings for allergic reactions (not disorders)

# Set of symptoms from Long Covid study (symptoms mapped to Read)
# Load all SNOMED CT concepts and Read maps

snomed_symptoms <- as.SNOMEDconcept(unique(merge(
	READ_SYMPTOMS[, .(read2 = readcode)],
	READMAPS[, .(read2 = unlist(read2_code)), by = conceptId],
	by = 'read2')$conceptId))
cat('\nSNOMED symptoms:', length(snomed_symptoms))
#snomed_symptoms <- union(snomed_symptoms, descendants(snomed_symptoms))
#cat('SNOMED symptoms with descendants:', length(snomed_symptoms))

# Add specific concepts of interest
extra <- descendants(SNOMEDconcept(c(
'Anaphylaxis',
'Diarrhoea',
'Constipation',
'Nausea',
'Vomiting',
'Upset stomach',
'Stomach ache',
'Indigestion',
'Dyspepsia',
'Gastrointestinal bleeding',
'Rash',
'Itching',
'Hives',
'Shortness of breath',
'Swelling',
'Oedema')), include_self = TRUE)

#setdiff(extra, snomed_symptoms)
snomed_symptoms <- union(snomed_symptoms, extra)
cat('\nSymptoms with some specific additions:', length(snomed_symptoms))

# Limit to the findings hierarchy, and add adverse reactions and allergic reactions
findings <- descendants('Clinical finding')
adverse_reactions <- descendants(c('Adverse reaction',
	'Adverse drug interaction', 'Application site disorder', 'Dizziness due to drug',
	'Drug dependence', 'Drug-induced feminization', 'Drug-induced lesion',
	'Drug-induced virilization'), include_self = TRUE)
allergic_reactions <- SNOMEDconcept(SNOMED$RELATIONSHIP[active == TRUE &
	typeId %in% SNOMEDconcept(c('Due to', 'Pathological process')) &
	destinationId %in% SNOMEDconcept(c('Allergic reaction',
	'Allergic process', 'Atopic reaction'))]$sourceId)

# Create a list of symptoms and other allergic reactions
reaction_concepts <- intersect(snomed_symptoms, findings)
cat('\nSNOMED symptoms limited to findings:', length(reaction_concepts))
reaction_concepts <- union(reaction_concepts, adverse_reactions)
reaction_concepts <- union(reaction_concepts, allergic_reactions)
cat('\nAfter adding adverse & allergic reactions:', length(reaction_concepts))

# Remove disorders with a 'causative agent'
CAUSATIVE_AGENT <- SNOMED$RELATIONSHIP[active == TRUE &
	typeId %in% SNOMEDconcept('Causative agent (attribute)'),
	.(conceptId = sourceId, substanceId = destinationId)]
reaction_concepts <- reaction_concepts[!reaction_concepts %in%
	CAUSATIVE_AGENT$conceptId]
cat('\nAfter removing causative agents:', length(reaction_concepts))

REACTIONS <- description(reaction_concepts, include_synonyms = TRUE)[
	type  == 'Synonym', .(conceptId, term, category = 'reaction')]

# Remove scales and scores
exclude_scale <- REACTIONS[term %like%
	'Scale: grade|Scale grade|scale grade|Score|score']$conceptId
REACTIONS <- REACTIONS[!(conceptId %in% exclude_scale)]

# Add extra synonyms without 'On examination' and 'Complaining of'
OE_CO <- REACTIONS[term %like%
	'^On examination ([-] |)|^Complaining of ']
OE_CO[, term := sub(
	'^On examination ([-] |)|^Complaining of ', '', term)]

REACTIONS <- rbind(REACTIONS, OE_CO)

# Remove specific concepts that are not meaningful
exclude <- SNOMEDconcept(c(
'164462000',
'164463005',
'164650007',
'164688009',
'164690005',
'164691009',
'164692002',
'164306004',
'164285001',
'164299002',
'164169001',
'274303007',
'267104002',
'272027003'
))

# exclude
# [1] "164462000 | On examination - nails (finding)"                          
# [2] "164463005 | On examination - nails - no abnormality detected (finding)"
# [3] "164650007 | On examination - sign attached to organ (finding)"         
# [4] "164688009 | On examination - sign painful (finding)"                   
# [5] "164690005 | On examination - sign slightly painful (finding)"          
# [6] "164691009 | On examination - sign moderately painful (finding)"        
# [7] "164692002 | On examination - sign very painful (finding)"              
# [8] "164306004 | On examination - character of fever (finding)"             
# [9] "164285001 | On examination - fever - general (finding)"                
#[10] "164299002 | On examination - level of fever (finding)"                 
#[11] "164169001 | On examination - Hess test for purpura (finding)"          
#[12] "274303007 | On examination - lymph nodes (finding)"                    
#[13] "267104002 | Complaining of a pain (finding)"                           
#[14] "272027003 | Complaining of a headache (finding)"  

REACTIONS <- REACTIONS[!(conceptId %in% exclude)]
REACTIONS <- REACTIONS[!duplicated(REACTIONS)]

M <- rbind(M, REACTIONS, fill = TRUE)
M <- M[!duplicated(M)]

# Select primary concept (the one with the longest and most descriptive name)
M <- M[order(category, conceptId, -nchar(term))]
M[, name_status := c('P', rep('A', .N - 1)), by = conceptId]
# Ideally change this to the language Refset

# Add terms with punctuation (except .) converted to spaces and remove multi spaces
punc_to_spaces <- function(text){
	text <- gsub('([0-9])-([0-9])', '\\1HHYYPPHHEENN\\2', text)
	text <- gsub('[^[:alnum:]. ]', ' ', text)
	text <- gsub(' +', ' ', text)
	gsub('HHYYPPHHEENN', '-', text)
}

# Add terms with punctuation (except .) removed and remove multi spaces
punc_to_remove <- function(text){
	text <- gsub('([0-9])-([0-9])', '\\1HHYYPPHHEENN\\2', text)
	text <- gsub('[^[:alnum:]. ]', '', text)
	text <- gsub(' +', ' ', text)
	gsub('HHYYPPHHEENN', '-', text)
}

# Convert '+' to '/' for combination medications
convert_plus_to_slash <- function(text){
	gsub(' \\+ ', ' / ', text)
}

M <- rbind(M, copy(M)[, term := punc_to_spaces(term)],
	copy(M)[, term := punc_to_remove(term)],
	copy(M)[substance %in% c('vmp', 'vtm'), term := convert_plus_to_slash(term)])

# Trim excess spaces
M[, term := sub('^ +', '', sub(' +$', '', term))]

# Remove duplicates
M <- M[!duplicated(M)]

cat('\nFinal categories:\n')
print(M[, .N, by = .(category, dosecategory)])

#### VALID MEDICATION LOOKUP ####

# Is a valid medication i.e. not a substance or reaction
fwrite(data.table(a = sort(unique(M[category %in% c('vmp', 'vtm')]$conceptId))),
	file = paste0(OUTPUT, 'valid_meds.csv'), col.names = FALSE)

#### REACTION CONVERSION LOOKUP

# Reaction subset: the subset of reactions that Epic accommodates
# (all other reactions are mapped to 'Allergic reaction') 
R <- M[category == 'reaction' & name_status == 'P',
	.(reactionsId = conceptId,
	subsetId = SNOMEDconcept('Allergic reaction'), reactionTerm = term)]
R <- R[!duplicated(R)]

# Diarrhoea SCT 62315008 "Diarrhoea"
d <- SNOMEDconcept('62315008')
R[reactionsId %in% c(d, descendants(d)), subsetId := d]

# Nausea and/or vomiting  SCT 16932000 "Nausea and vomiting"
d <- SNOMEDconcept('16932000')
e <- c(SNOMEDconcept('Vomiting'), SNOMEDconcept('Nausea'))
R[reactionsId %in% c(d, descendants(d), e, descendants(e)), subsetId := d]

# Rash, itching or hives SCT 126485001 "Urticaria"
d <- SNOMEDconcept('126485001')
e <- c(SNOMEDconcept('Eruption of skin'), SNOMEDconcept('Itching'),
	SNOMEDconcept('Wheal'))
R[reactionsId %in% c(d, descendants(d), e, descendants(e)), subsetId := d]

# Gastrointestinal bleeding SCT 74474003 "Gastrointestinal haemorrhage"
d <- SNOMEDconcept('74474003')
e <- c(SNOMEDconcept('Haematemesis'), SNOMEDconcept('Melaena'))
R[reactionsId %in% c(d, descendants(d), e, descendants(e)), subsetId := d]

# Anaphylaxis SCT 39579001 "Anaphylaxis"
d <- SNOMEDconcept('39579001')
R[reactionsId %in% c(d, descendants(d)), subsetId := d]

# Shortness of breath SCT 267036007 "Dyspnoea"
d <- SNOMEDconcept('267036007')
R[reactionsId %in% c(d, descendants(d)), subsetId := d]

# Swelling SCT 65124004 "Swelling"
d <- SNOMEDconcept('65124004')
e <- c(SNOMEDconcept('Swelling'), SNOMEDconcept('Angioedema'))
R[reactionsId %in% c(d, descendants(d), e, descendants(e)), subsetId := d]

# Other (see comments) SCT 419076005,"Allergic reaction"
# All other reactions

R[, subsetTerm := description(subsetId)$term, by = subsetId]
fwrite(R, file = paste0(WORKING, 'reactions_subset_check.csv'))

fwrite(R[, .(reactionsId, subsetId)], file = paste0(OUTPUT, 'reactions_subset.csv'))

#### ALLERGY TYPE

# Allergy type
#              NAME    N
#1:      Drug Class  245
#2: Drug Ingredient 2425
#3:            Drug  762
#4:            Food   41
#5:   Environmental   12
#6:          Animal    7
#7:        Chemical    1

A <- fread('
AdverseReactionTerm|adverseReactionId|allergenType|adverseReactionType
Food Allergy|414285001|Food|Allergy
Drug Allergy|416098002|Drug|Allergy
Drug Allergy|416098002|Drug Class|Allergy
Drug Allergy|416098002|Drug Ingredient|Allergy
Allergy to substance|419199007|Chemical|Allergy
Environmental allergy|426232007|Environmental|Allergy
Environmental allergy|426232007|Animal|Allergy
Food Intolerance|235719002|Food|Intolerance
Drug Intolerance|59037007|Drug|Intolerance
Drug Intolerance|59037007|Drug Class|Intolerance
Drug Intolerance|59037007|Drug Ingredient|Intolerance
Propensity to adverse reactions to substance|418038007|Chemical|Intolerance
Propensity to adverse reactions|420134006|Environmental|Intolerance
Propensity to adverse reactions|420134006|Animal|Intolerance
Propensity to adverse reactions to food|418471000|Food|Unspecified
Propensity to adverse reactions to drug|419511003|Drug|Unspecified
Propensity to adverse reactions to drug|419511003|Drug Class|Unspecified
Propensity to adverse reactions to drug|419511003|Drug Ingredient|Unspecified
Propensity to adverse reactions to substance|418038007|Chemical|Unspecified
Propensity to adverse reactions|420134006|Environmental|Unspecified
Propensity to adverse reactions|420134006|Animal|Unspecified
')

# adverseReactionType = Allergy, Intolerance, Unspecified

fwrite(A[, .(allergenType, adverseReactionType, adverseReactionId)],
	paste0(OUTPUT, 'allergy_type.csv'))

#### VTM TO VMP LOOKUP

# VTM to VMP: select VMP with appropriate dose

fwrite(VTM_VMP[keep == TRUE, .(vtmId, dose,
	unit = units, vmpId = conceptId)][order(vtmId, dose)],
	paste0(OUTPUT, 'vtm_to_vmp.csv'))

# VTM to text: VTM convert to 'XXX ORAL' text if possible, as a 
# fallback

# Get list of concept : name matches
ORAL <- ORAL[!(Name %like% ' \\(GENERIC MANUF\\)'), .(Name, conceptId)]
ORAL <- ORAL[!duplicated(ORAL)]
# Check concept ID is unique
stopifnot(uniqueN(ORAL$conceptId) == nrow(ORAL))
M[, ORALNAME := ORAL[M, on = 'conceptId']$Name]

# Exclude common oral meds that are usually taken by another route
# e.g. salbutamol
M[ORALNAME %in% c('BUDESONIDE ORAL', 'SALBUTAMOL ORAL',
	'SODIUM BENZOATE ORAL', 'POTASSIUM IODIDE ORAL',
	'SEMAGLUTIDE ORAL', 'BECLOMETASONE ORAL'), ORALNAME := NA]

fwrite(M[name_status == 'P' & !is.na(ORALNAME),
	.(conceptId, term, ORALNAME)],
	paste0(WORKING, 'vtm_to_text_check.csv'))

fwrite(M[name_status == 'P' & !is.na(ORALNAME),
	.(vtmId = conceptId, text = ORALNAME)],
	paste0(OUTPUT, 'vtm_to_text.csv'))

#### CREATION OF LOOKUPS IN CORRECT FORMAT ####

# Output files
fwrite(M, paste0(WORKING, 'med_working_table_full.csv'))
fwrite(CAUSATIVE_AGENT, paste0(WORKING, 'causative_agent_lookup.csv'))
fwrite(INGREDIENTS, paste0(WORKING, 'ingredients_lookup.csv'))

TEMP <- description(SUBS$conceptId)[, .(conceptId, term)]
TEMP <- TEMP[!duplicated(TEMP)]
TEMP <- TEMP[order(conceptId)]
fwrite(TEMP, paste0(WORKING, 'all_substances_lookup.csv'))
rm(TEMP)

#### EXCLUDE INAPPROPRIATE CONCEPTS AND WRITE OUT FINAL CDB

TEMP <- merge(MED_SUB_REACTION_EXCLUDE[is.na(conceptId), .(name)],
	SNOMED$DESCRIPTION[, .(conceptId, name = tolower(term))], by = 'name')
TEMP <- rbind(TEMP, MED_SUB_REACTION_EXCLUDE[!is.na(conceptId), .(conceptId, name)])
TEMP[, term := description(conceptId)$term[1], by = conceptId]
TEMP[, semanticType := semanticType(conceptId)]
# use substance, disorder or finding for exclusion list
to_exclude <- as.SNOMEDconcept(TEMP[semanticType %in% c('substance', 'disorder',
	'finding')]$conceptId)
M[, exclude := conceptId %in% to_exclude]
rm(TEMP)

# Remove duplicates
M <- M[!duplicated(M)]

# Write out final CDB
OUT <- M[exclude == FALSE, .(cui = conceptId, name = tolower(term),
	ontologies = 'SNO', name_status)][order(cui, name_status, name)]
OUT <- OUT[!duplicated(OUT)]
fwrite(OUT, paste0(OUTPUT, 'cdb_medallerg.csv'))

#### MAP ALLERGENS TO SUBSET ####

# Map allergens to subset: maps vtm, vmp etc. to substance
# and then substance to parent substance to 

# Correction to 'alcohol' (food) - map to 'alcoholic beverage' instead 
# EXM[ALLERGEN_ID == 88014, SNOMEDCT := 53527002]
# Only make the correction here if it is updated in Epic itself
S2 <- EXM[, .(conceptId = as.SNOMEDconcept(strsplit(SNOMEDCT, '\n')[[1]])),
	by = ALLERGEN_ID]
S2 <- S2[!duplicated(S2)]

# First map everything to a substance
S1 <- M[category %in% c('vtm', 'vmp', 'substance') & name_status == 'P']

# Identify the substance if a single ingredient
INGREDIENTS[, one := .N == 1, by = conceptId]
S1[, substanceId := INGREDIENTS[one == TRUE][S1, on = 'conceptId']$ingredientId]
S1[category == 'substance', substanceId := conceptId]
# S1[, substanceName := description(substanceId)$term[1], by = substanceId]

# Generate the descendants of each substance 
SLOOKUP <- merge(S2, S2[, .(subtree = list(c(conceptId, descendants(conceptId)))),
	by = conceptId], by = 'conceptId', all.x = TRUE)
SLOOKUP[, nsub := sapply(subtree, length)]

# Apply to master allergen list - in order from large to small
# this means that more specific substances should overwrite broader substances
SLOOKUP <- SLOOKUP[order(-nsub)]

S1[, ALLERGEN_ID := as.integer(NA)]
S1[, subsetId := bit64::as.integer64(NA)]
# This loop takes a long time (~1/2 hour) to run
for (i in 1:nrow(SLOOKUP)){
	S1[substanceId %in% SLOOKUP[i]$subtree[[1]],
		ALLERGEN_ID := SLOOKUP[i]$ALLERGEN_ID]
	S1[substanceId %in% SLOOKUP[i]$subtree[[1]],
		subsetId := SLOOKUP[i]$conceptId]
}

# Mapping to generic 'Substance' is not allowed as it is not informative
S1[subsetId == SNOMEDconcept('Substance (substance)'), subsetId := NA]

# Now merge the allergen type
S1[, allergenType := ELG[S1, on = 'ALLERGEN_ID']$NAME]

# Add 'Drug' if allergen type is blank
S1[is.na(allergenType) & subsetId %in% descendants('Drug or medicament'),
	allergenType := 'Drug']
S1[is.na(allergenType) & subsetId %in% descendants('Chemical'),
	allergenType := 'Chemical']

# Export
fwrite(S1[!is.na(subsetId), .(allergensId = conceptId, subsetId,
	allergenType)], paste0(OUTPUT, 'allergens_subset.csv'))

fwrite(S1, paste0(WORKING, 'allergens_subset_check.csv'))
# allergensId, subsetId, allergenType	csv	

cat('\nDONE\n')
sink()

