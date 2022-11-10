# Script to create additional synonyms for SNOMED CT concepts for
# MedCat concept database for MiADE

# Anoop Shah, October 2022
# https://docs.google.com/document/d/1suIa4NoFZ92rwFzPi6juNkqWFQ-0tGtng8nxR874lJs/edit?pli=1#

library(data.table)
library(Rdiagnosislist)

# To use:
# Substances hierarchy
# Product VMPs
# Product VTMs
# Remove substances that are also VTMs (to avoid duplication)

# Extract dose and strength from VMPs

# Log file
sink(paste0('Generate_MiADE_meds_allergy_cdb_v4', Sys.time(), '.log'), split = TRUE)

######################################
# Import relevant files

# Import SNOMED CT dictionary (created using loadSNOMED function)
# Using May 2022 SNOMED dictionary and Read maps
SNOMED <- readRDS('~/Terminologies/SNOMED_May2022.RDS')
READMAPS <- readRDS('~/Terminologies/READMAPS.RDS')
READ_SYMPTOMS <- fread('Symptom_codelists.csv') # derived from Long Covid project

######################################

SUBS <- setdiff(descendants('Substance'), children('Substance')) # 27782
SUBS <- description(SUBS, include_synonyms = TRUE)[
	type  == 'Synonym', .(conceptId, term, category = 'substance')]

# Ensure that environmental allergens are included as in Epic ELG file
# PERMEABLE ADHESIVE TAPES
# TAPE AND TAPE GLUE
# DRESSINGS AND DRESSING ADHESIVE
# ADHESIVE/SURGICAL TAPE AND SURGICAL STRIPS

# Closest match for these is:
SNOMEDconcept('Adhesive agent')
SNOMEDconcept('Glue')
# Otherwise there are lots of specific options. We should not map to
# 'physical object' individual dressings 

VTM <- children('Virtual therapeutic moiety') # 2944
VTM <- description(VTM, include_synonyms = TRUE)[
	type  == 'Synonym', .(conceptId, term, category = 'vtm')]

VMP <- children('Virtual medicinal product') # 18838
VMP <- description(VMP, include_synonyms = TRUE)[
	type  == 'Synonym', .(conceptId, term, category = 'vmp')]

M <- rbind(SUBS, VTM, VMP)

# Fully specified name for (clinical drug) may be never used in text
# e.g. "Product containing precisely bendroflumethiazide 2.5 milligram/1 each
#       conventional release oral tablet (clinical drug)"
# Hence keep only synonyms for clinical drugs

# Remove space between number and units
M[, term := sub('([0-9]) (ml|nanogram|nanograms|g|mg|microgram|micrograms)',
	'\\1\\2', term)]

# Add synonyms for modified-release = MR
M <- rbind(M, M[term %like% 'modified-release' & category == 'vmp',
	.(conceptId, term = sub('modified-release', 'MR', term), category = 'vmp')])

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
M[category == 'vmp', dosecategory := sub('^.*[0-9]+[^\\.].*[0-9]+.*$',
	'MULTIPLE_NUMBERS', term)]
M[category == 'vmp' & dosecategory != 'MULTIPLE_NUMBERS',
	dosecategory := ifelse(term %like% ' (oral |)(capsule|tablet)(s|)$',
	'SIMPLE_TABLET', 'COMPLEX')]
M[dosecategory == 'SIMPLE_TABLET', dose := 
	as.numeric(sub('^.*[^0-9]([0-9\\.]+)[ ]*(ml|nanogram|nanograms|g|mg|microgram|micrograms).*$',
	'\\1', term))]
M[dosecategory == 'SIMPLE_TABLET' & !is.na(dose), units :=
	sub('^.*[^0-9]([0-9\\.]+)[ ]*(ml|nanogram|nanograms|g|mg|microgram|micrograms).*$',
	'\\2', term)]

MDOSE <- M[!is.na(dose), .(conceptId, dose, units)]
MDOSE <- MDOSE[!duplicated(MDOSE)]
MDOSE[, Ndoses := .N, by = conceptId]
# Limit to entries with a single unambiguous dose
MDOSE <- MDOSE[Ndoses == 1]
MDOSE[, Ndoses := NULL]

# Create a VTM to VMP conversion tool based on dose
# VTM is parent of VMP and also is_a VTM
# Extract table of parents
# Intersect parents with VTMs
# Allow only one VTM per VMP
VTM_VMP <- merge(M[category == 'vmp', .(conceptId)],
	SNOMED$RELATIONSHIP[active == TRUE &
	typeId == SNOMEDconcept('Is a'), .(conceptId = sourceId,
	vtmId = destinationId)])[vtmId %in% VTM$conceptId]
VTM_VMP <- VTM_VMP[!duplicated(VTM_VMP)]
VTM_VMP[, dose := MDOSE[VTM_VMP, on = 'conceptId']$dose]
VTM_VMP[, units := MDOSE[VTM_VMP, on = 'conceptId']$units]
VTM_VMP <- merge(VTM_VMP, description(VTM_VMP$conceptId)[, .(conceptId, term)])


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
M <- M[order(category, dosecategory, conceptId, -nchar(term))]
M[, name_status := c('P', rep('A', .N - 1)), by = conceptId]

cat('\nFinal categories:\n')
print(M[, .N, by = .(category, dosecategory)])

# Output files
fwrite(description(SUBS$conceptId)[, .(conceptId, term)],
	'MiADE_all_substances_lookup.csv')
fwrite(CAUSATIVE_AGENT, 'MiADE_causative_agent_lookup.csv')
fwrite(INGREDIENTS, 'MiADE_ingredients_lookup.csv')
fwrite(VTM_VMP, 'MiADE_VTM_VMP_dose_lookup.csv')
fwrite(M, 'MiADE_med_allergy_reaction_lookup.csv')
fwrite(M[, .(cui = conceptId, name = tolower(term), ontologies = 'SNO',
	name_status)], 'med_allergy_reaction_MedCAT_cdb.csv')

cat('\nDONE\n')
sink()
