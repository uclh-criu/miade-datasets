# Script to create additional synonyms for MedCat concept database
# for disorder / finding terms.

# Format of output file for MedCAT:
# cui, name, ontologies, name_status

#### LOAD FILES ####

# Create a script for conversion of Read to SNOMED
library(data.table)
library(Rdiagnosislist)
library(parallel)

defaultfilepath <- '~/pCloudDrive/anoop_work/2019to2021/MiADE/Terminologies/'
if (!exists('WORKING')){ WORKING <- paste0(defaultfilepath, 'WORKING/') }
if (!exists('OUTPUT')){ OUTPUT <- paste0(defaultfilepath, 'OUTPUT/') }
if (!exists('REF')){ REF <- paste0(defaultfilepath, 'REF/') }
if (!exists('LOGS')){ LOGS <- paste0(defaultfilepath, 'LOGS/') }

# Log file
sink(paste0(LOGS, 'gen_cdb_problems_v2_', Sys.time(), '.log'),
	split = TRUE)

file_extras <- paste0(WORKING, 'extra_synonyms_28sep.csv')

# Import SNOMED CT dictionary (created using loadSNOMED function)
if (!exists('SNOMED')){
	SNOMED <- readRDS('~/Terminologies/SNOMED_May2022.RDS')
}

# Manual synonyms (e.g. for mandatory comordbidities)
MAN <- rbind(fread(paste0(REF, 'mandatory_comorbidities.csv')),
	fread(paste0(REF, 'manual_synonyms.csv')))
	
# SNOMED CT code usage
FREQ <- fread(paste0(REF, 'SNOMED_code_usage_2021_22.txt'))

# Load up the language refset
# Identify UK preferred synonyms
# https://confluence.ihtsdotools.org/display/DOCRFSPG/5.8.+Language+Reference+Set
snomed_filepath <- '~/Terminologies/May2022/'
x <- list.files(snomed_filepath, recursive = TRUE, full.names = TRUE)
lang_refset <- x[x %like% 'Refset_LanguageUK.+Snapshot']

#### PROCESS DICTIONARIES ####

LANG <- data.table(id = character(0), effectiveTime = integer(0),
	active = integer(0), moduleId = bit64::integer64(0),
	refsetId = bit64::integer64(0),
	referencedComponentId = bit64::integer64(0),
	acceptabilityId = bit64::integer64(0))
for (i in lang_refset) LANG <- rbind(LANG, fread(i))
# Keep only the most recent value for active concepts
LANG[, keep := active == 1 & effectiveTime == max(effectiveTime),
	by = referencedComponentId]
LANG <- LANG[keep == TRUE]
LANG <- LANG[!duplicated(LANG)]

# Check that referencedComponentId appears only once
# stopifnot(nrow(LANG) == length(unique(LANG$referencedComponentId)))

# Mark the UK preferred terms
sct_pref <- bit64::as.integer64(SNOMEDconcept('Preferred'))
SNOMED$DESCRIPTION[, UKpref := 
	LANG[, .(UKpref = any(acceptabilityId == sct_pref)),
	by = .(id = referencedComponentId)][
	SNOMED$DESCRIPTION, on = 'id']$UKpref]
synonym_type <- SNOMEDconcept('Synonym')
SNOMED$DESCRIPTION[, Synonym := typeId == synonym_type]

# Select a single UK preferred concept
setkey(SNOMED$DESCRIPTION, conceptId, UKpref, Synonym)
SNOMED$DESCRIPTION[, choose := c(rep(FALSE, .N - 1), TRUE), by = conceptId]

#### SUSPECTED, HISTORIC, NEGATED, BLACKLIST ####

# Suspected, Historic, Negated version of concepts (where available)
relation_source <- function(destination){
	sort(unique(as.SNOMEDconcept(SNOMED$RELATIONSHIP[
		destinationId == SNOMEDconcept(destination)]$sourceId)))
}

# Procedures can be used for historic concepts for problem list
# but not for current procedures (yet)
disorders <- descendants('Disorder')
acute_diseases <- descendants('Acute disease')
procedures <- descendants('Procedure')
suspected <- relation_source('Suspected')
known_absent <- relation_source('Known absent')
known_present <- relation_source('Known present')
subject_of_record <- relation_source('Subject of record')
family_person <- relation_source('Person in family of subject')
in_the_past <- relation_source('In the past')
done <- relation_source('Done')
current_or_specified_time <- relation_source('Current or specified time')

allergy <- descendants('Propensity to adverse reactions to substance')
ignore_findings <- SNOMEDconcept(c('Disease', 'Clinical finding', 'Problem',
	'Impairment', 'Chief complaint', 'Sign', 'Complaint', 'Sequela',
	'Early complication', 'Co-morbid conditions', 'Pre-existing condition',
	'Acute disease', 'Subacute disease', 'Chronic disease',
	'General problem AND/OR complaint', 'Evaluation finding',
	'Administrative statuses', 'Finding by site', 'Finding by method',
	'Clinical history and observation findings', 'Behaviour',
	'Adverse incident outcome categories', 'Prognosis/outlook finding',
	'General clinical state finding', 'Disorder by body site',
	'Failure', 'Acute failure', 'Subacute failure', 'Chronic failure',
	'Decompensation', 'Discrepancy', 'Idiosyncrasy', 'Inefficiency',
	'General body state finding', 'General clinical state finding',
	'Pressure', 'Disease related state', 'Absence of pressure',
	'Decreased pressure', 'Increased pressure', 'Swelling',
	'Disease condition finding', 'Allergic disposition', 'Pain',
	'Values (community)', 'Fit and well', 'No sensitivity to pain'))

ASSOC_FINDPROC <- SNOMED$RELATIONSHIP[
	typeId %in% SNOMEDconcept(c('Associated finding', 'Associated procedure')),
	.(situationId = sourceId, sitTerm = description(sourceId)$term,
	findingId = destinationId, finTerm = description(destinationId)$term)]

HISTORIC <- ASSOC_FINDPROC[situationId %in% intersect(c(known_present, done),
	intersect(subject_of_record, in_the_past)) &
	findingId %in% c(procedures)][order(findingId)]

# Currently only mapping to procedures - 27 Oct 2023
#	findingId %in% c(acute_diseases, procedures)][order(findingId)]

# Keep if only one situation per finding/procedure, to avoid errors
# e.g. 'Implantation procedure' --> 'History of heart valve recipient'
# For diseases we are only interested in acute diseases
# not 'history of hypothyrodism etc.'
HISTORIC[, Nsit := .N, by = findingId]
HISTORIC <- HISTORIC[Nsit == 1 & findingId != SNOMEDconcept('Aftercare')]

NEGATED <- ASSOC_FINDPROC[situationId %in% intersect(known_absent,
	intersect(subject_of_record, current_or_specified_time))][order(findingId)]
# Keep if only one situation per finding, to avoid errors
NEGATED[, Nsit := .N, by = findingId]
NEGATED <- NEGATED[Nsit == 1]

SUSPECTED <- ASSOC_FINDPROC[situationId %in% intersect(suspected,
	intersect(subject_of_record, current_or_specified_time))][order(findingId)]
# Keep if only one situation per finding, to avoid errors
SUSPECTED[, Nsit := .N, by = findingId]
SUSPECTED <- SUSPECTED[Nsit == 1]

# List of suspected, negated, historic, allergy, family history and
# procedures that do not have an associated historic concepts,
# to exclude from MedCAT annotation training
exclude_conceptId <- sort(unique(c(suspected, in_the_past, known_absent,
	family_person, allergy, ignore_findings,
	setdiff(procedures, ASSOC_FINDPROC$findingId))))

#### FUNCTIONS FOR ADDING SPECIFIC SNOMED CT SYNONYMS ####

ACRONYMS <- fread('
acronym|expansion
af|atrial fibrillation
aaa|abdominal aortic aneurysm
cad|coronary artery disease
ccf|congestive cardiac failure
chd|coronary heart disease
dm|diabetes mellitus
flu|influenza
pah|pulmonary artery hypertension
mi|myocardial infarction
ms|multiple sclerosis
pid|pelvic inflammatory disease
rapd|relative afferent pupillary defect
nfa|no fixed abode
tb|tuberculosis
uc|ulcerative colitis
aom|acute otitis media
itp|immune thrombocytopenia
nof|fracture of neck of femur
svt|supraventricular tachycardia
ibs|irritable bowel syndrome
aki|acute renal failure syndrome
sah|subarachnoid intracranial hemorrhage
hf|heart failure
psp|progressive supranuclear palsy
vf|ventricular fibrillation
htn|hypertension
cva|cerebrovascular accident
cvd|cardiovascular disease
mca|middle cerebral artery
aca|anterior cerebral artery
pca|posterior cerebral artery
lad|left anterior descending
rca|right coronary artery
nof|neck of femur
oa|osteoarthrosis
oa|osteoarthritis
hoh|hard of hearing
hoh|hard hearing
gord|gastroesophageal reflux disease
ca|cancer
pmr|polymyalgia rheumatica
gist|gastrointestinal stromal tumor
gist|gastrointestinal stromal tumour
gi|gastrointestinal
git|gastrointestinal tract
gih|gastrointestinal haemorrhage
hcc|hepatocellular carcinoma
stemi|st segment elevation myocardial infarction
nstemi|non-st segment elevation mi
t2dm|diabetes mellitus type 2
hfref|heart failure with reduced ejection fraction
hfpef|heart failure with preserved ejection fraction
hfmef|heart failure with mid range ejection fraction')

#### LOAD EXTRA SYNONYMS ####

EXTRA_SYNONYMS <- fread(file_extras)[!is.na(synonym) & synonym != '']
# How many were sucessfully converted
print(uniqueN(EXTRA_SYNONYMS$conceptId))
print(nrow(EXTRA_SYNONYMS))

#### FUNCTIONS TO ADD ACRONYMS ####

# Add acronyms for SNOMED CT disorders, findings and events
get_acronym <- function(terms){
	# check for pattern ABCD - Another Bland Cardiovascular Disease
	# allow for word order to be different and for one word to be omitted
	# or vice versa
	acronym <- rep(NA_character_, length(terms))
	words <- strsplit(tolower(terms), split = ' ')
	maybe_acronym <- sapply(words, length) >= 3 & terms %like% '^([[:alnum:]]+) - '
	stated_acronym <- strsplit(sapply(words, function(x) x[1]),
		split = character(0))
	stated_expansion <- sub('^([[:alnum:]]+) - (.*)$', '\\2', terms)
	created_acronym <- lapply(words, function(x){
		substr(x[3:length(x)], 1, 1)
	})
	
	if (any(maybe_acronym)){
		for (i in seq_along(acronym)[maybe_acronym == TRUE]){
			if (((length(setdiff(stated_acronym[[i]], created_acronym[[i]])) <= 1 |
				(length(setdiff(stated_acronym[[i]], created_acronym[[i]])) <= 2 &
				length(stated_acronym[[i]]) > 1 & length(created_acronym[[i]]) > 1)) &
				(length(setdiff(created_acronym[[i]], stated_acronym[[i]])) <= 1)) |
				(nrow(ACRONYMS[acronym == words[[i]][1] &
				expansion == stated_expansion[i]]) > 0)){
				acronym[i] <- words[[i]][1]
			}
		}
	}
	
	acronym
}

# Testing
terms <- c('ca - cancer', 'nof - fracture of neck of femur',
	'abc - aa bb cc dd', 'abc - aa bb', 'abc - def ghi', 'abc',
	'ivdu - intravenous drug user', 'oe - fever',
	'nstemi - non st elevation mi', 'non-venomous')
print(get_acronym(terms))
# should be 
#[1] "ca"     "nof"    "abc"    "abc"    NA       NA       "ivdu"   NA      
#[9] "nstemi"

remove_stopwords2 <- function(terms,
	stopwords = c('the', 'of', 'from', 'with')){
	# remove [X] preface
	text <- paste0(' ', tolower(sub('^\\[X\\]', '', terms)), ' ')
	# remove stopwords
	for (word in stopwords){
		text <- gsub(paste0(' ', word, ' '), ' ', text) 
	}
	# remove leading and trailing spaces
	text <- sub('^ +', '', sub(' +$', '', text))
	# remove acronyms (format abc - aa bb cc)
	acronyms <- get_acronym(text)
	text[!is.na(acronyms)] <- sub('^([[:alnum:]]+) - (.*)$', '\\2',
		text[!is.na(acronyms)])
	# remove non-alphanumeric characters and multiple spaces
	text <- gsub('[^[:alnum:] ]', ' ', text)
	text <- gsub(' +', ' ', text)
	text
}
print(remove_stopwords2(terms))
# should be
#[1] "cancer"                "facture neck femur"    "aa bb cc dd"          
#[4] "aa bb"                 "abc def ghi"           "abc"                  
#[7] "intravenous drug user" "oe fever"              "non st elevation mi"  

remove_hyphens <- function(terms){
	gsub('-','', terms)
}

add_rows <- function(SCT, NEW){
	message('SCT has ', nrow(SCT), ' rows and NEW has ', nrow(NEW), ' rows.')
	SCT <- rbind(SCT[, .(conceptId, term)], NEW[, .(conceptId, term)])
	SCT <- SCT[!duplicated(SCT)]
	message(nrow(SCT), ' rows after deduplication.')
	SCT
}

acronym_rows <- function(SCT, outfile = paste0(WORKING, 'acronyms_to_check.csv')){
	SCT[, acronym := get_acronym(term)]

	# If the acronym is identical to an actual term, add the term as well
	# e.g. COLD = cold or 'chronic obstructive lung disease'
	# This is important for disambiguation
	all_acronyms <- sort(unique(SCT$acronym))
	all_acronyms <- all_acronyms[!is.na(all_acronyms)]
	
	# If term is already an acronym
	SCT[tolower(term) %in% all_acronyms, acronym := tolower(term)]

	# Create a table of acronyms (acronym capitalised)
	A <- SCT[!is.na(acronym), .(conceptId,
		term = acronym, origterm = tolower(term))]
	A[, nconcepts := uniqueN(conceptId), by = term]
	A[, nterms := uniqueN(tolower(origterm)), by = term]
	if (!is.null(outfile)){
		fwrite(A[(nterms > 1 & nconcepts > 1)][order(term)],
			file = outfile)
	}
	SCT[, acronym := NULL]
	
	# Add acronyms where there is only one unique concept linked to the
	# acronym 
	A[nconcepts == 1][, .(conceptId, term)]
}

#### PROCESS SNOMED CONCEPTS - ADD SYNONYMS ####

# Initial set of concepts of interest for problems using UK preferred synonyms
concepts_to_use <- setdiff(c(descendants('Clinical finding'),
	descendants('Event'), procedures,
	descendants('Situation with explicit context')), exclude_conceptId)

SCT <- SNOMED$DESCRIPTION[conceptId %in% concepts_to_use &
	typeId == SNOMEDconcept('Synonym'),
	list(conceptId, origterm = term, name_status =
	ifelse(choose, 'P', 'A'))]
SCT[, term := sub('^\\[X\\]', '', tolower(origterm))]
# Note that procedures are included at this stage but not included
# in final output unless they are historic and have been converted
# into a 'history' concept
 
SCTPREF <- SCT[name_status == 'P']

message('SCT has ', nrow(SCT), ' rows (', uniqueN(SCT$conceptId),
	' concepts).')

# Verify only one preferred per concept Id
stopifnot(nrow(SCTPREF) == SCTPREF[, uniqueN(conceptId)])

# Add versions of concepts with hyphens removed
SCT <- add_rows(SCT, SCT[, .(conceptId, term = remove_hyphens(term))])

# Add versions of concepts with stopwords removed
SCT <- add_rows(SCT, SCT[, .(conceptId, term = remove_stopwords2(term))])

# Add unambiguous acronyms specified in SNOMED CT 
SCT <- add_rows(SCT, acronym_rows(SCT))

# Add manual synonyms
SCT <- add_rows(SCT, MAN[, .(term = sub('^ +| +$', '',
	unlist(strsplit(manual_synonyms, ',')))), by = conceptId])

# Add extra synonyms for most frequent concepts
SCT <- add_rows(SCT, EXTRA_SYNONYMS[, .(conceptId, term = synonym)])

# Mark which terms are preferred and restore original form
# (sentence case) for preferred terms
SCT2 <- merge(SCT[conceptId %in% SCTPREF$conceptId,
	.(conceptId, term)], SCTPREF, by = c('conceptId', 'term'), all = TRUE)
SCT2[is.na(name_status), name_status := 'A']
SCT2[name_status == 'P', term := origterm] 
SCT2[, origterm := NULL]
SCT2 <- SCT2[!duplicated(SCT2)]

# Check that each concept has exactly one preferred term
stopifnot(all(SCT2[, .(numP = sum(name_status == 'P')), by = conceptId]$numP) == 1)

# Remove all one character terms and remove 2-3 character terms which
# are common English words and not in the list of specified acronyms.
BRIF <- fread('https://raw.githubusercontent.com/anoopshah/freetext-matching-algorithm-lookups/master/2of4brif.txt',
	header = FALSE, col.names = 'word')

SHORT <- SCT2[name_status == 'A' & nchar(term) < 4, .(conceptId, term)]
SHORT <- SHORT[(term %in% BRIF[[1]] & !(term %in% ACRONYMS$acronym)) |
	nchar(term) == 1]
SHORT[, full := description(conceptId)$term[1], by = conceptId]

# Estimate words frequencies in SNOMED CT
allwords <- unlist(strsplit(tolower(SNOMED$DESCRIPTION$term), ' '))
WORDFREQ <- as.data.table(table(allwords))
SHORT[, wordfreq := WORDFREQ[, .(term = allwords, freq = N)][
	SHORT, on = 'term']$freq]
SHORT[is.na(wordfreq), wordfreq := 0]
SHORT[, freq:= FREQ[, .(conceptId = SNOMED_Concept_ID,
	freq = as.numeric(Usage))][SHORT, on = 'conceptId']$freq]
SHORT[is.na(freq), freq := 0]
SHORT <- SHORT[order(-wordfreq, term)]

fwrite(SHORT, file = paste0(WORKING, 'short_terms_to_exclude.csv'))

# Remove short terms not in the acronym list with 3 words and more than 10
# frequency among SNOMED term descriptions, or 1-2 words regardless of
# frequency.
SCT2 <- SCT2[!(term %in% SHORT[wordfreq > 10 | nchar(term) <= 2]$term &
	name_status == 'A')]

# Remove 'situation' concepts not in UCLH problem list subset
SUBSET <- fread(paste0(REF, 'kawsar_snomed_subset.csv'))
SUBSET[, semanticType := semanticType(cui)]
SCT2[, semanticType := semanticType(conceptId)]

SCT2[, exclude := FALSE]
SCT2[semanticType == 'situation' & !(conceptId %in% SUBSET$cui),
	exclude := TRUE]
SCT2 <- SCT2[exclude == FALSE]
SCT2[, semanticType := NULL]
SCT2[, exclude := NULL]

# Remove generic body system findings
BODY_DESC <- description(descendants('Body structure'),
	include_synonyms = T)

body_desc_o <- c(paste(BODY_DESC$term, 'finding'),
	paste('Finding of', tolower(BODY_DESC$term)),
	paste(BODY_DESC$term, 'system finding'),
	paste('Finding of', tolower(BODY_DESC$term), 'system'),
	paste(BODY_DESC$term, 'observation'),
	paste('Observation of', tolower(BODY_DESC$term)),
	paste(BODY_DESC$term, 'system observation'),
	paste('Observation of', tolower(BODY_DESC$term), 'system'))

BODY_DESC_O <- SNOMED$DESCRIPTION[term %in% body_desc_o]
SCT2 <- SCT2[!(conceptId %in% BODY_DESC_O$conceptId)]
rm(BODY_DESC_O)
rm(body_desc_o)

# Remove social history (except housing problems and care needs),
# administrative statuses (except registered disabled) and normal findings
social_history <- descendants('Social and personal history finding',
	include_self = TRUE)
housing_and_care <- descendants(c('Homeless', 'No fixed abode',
	'Lives alone',
	'Lives in supported home', 'Unsatisfactory living conditions',
	'Finding related to care and support circumstances and networks'),
	include_self = TRUE)
admin <- descendants('Administrative statuses',
	include_self = TRUE)
regstatus <- SNOMEDconcept(c('On adult protection register',
	'On learning disability register',
	'On social services disability register',
	'Registered blind',
	'Registered partially sighted',
	'Registered deaf',
	'Registered hearing impaired',
	'Registered sight impaired'))
normal <- SNOMEDconcept('^Normal| normal', exact = FALSE)
normal <- normal[semanticType(normal) == 'finding']

SCT2 <- SCT2[!(conceptId %in% c(setdiff(admin, regstatus),
	setdiff(social_history, housing_and_care), normal))]
rm(admin)
rm(housing_and_care)
rm(social_history)
rm(normal)

#### EXPORT DATASETS ####

# Export full vocab for MiADE 
setkey(SCT2, conceptId)
fwrite(SCT2, file = paste0(WORKING, 'full_condition_list.csv'))

# Format of output file for MedCAT:
# cui, name, ontologies, name_status
fwrite(SCT2[, .(cui = conceptId, name = term,
	ontologies = 'SNO', name_status)], file = paste0(OUTPUT, 'cdb_problems.csv'))

# Remove any synonyms that are already in SNOMED CT
# (whether in original or lower case form)
OUT <- merge(SCT2[, .(conceptId, term)],
	rbind(SNOMED$DESCRIPTION[, .(conceptId, term = term,
	already = TRUE)], SNOMED$DESCRIPTION[, .(conceptId, term = tolower(term),
	already = TRUE)]), by = c('conceptId', 'term'), all.x = TRUE)
OUT <- OUT[is.na(already)]
OUT[, already := NULL]

# Export list of additional concepts
fwrite(OUT, file = paste0(WORKING, 'extra_condition_synonyms.csv'))

# Negated

fwrite(NEGATED[, .(findingId, situationId)][order(findingId)],
	file = paste0(OUTPUT, 'negated.csv')) 
fwrite(NEGATED, file = paste0(WORKING, 'negated_check.csv'))

# Historic (Oct 2023 - currently procedures only)

fwrite(HISTORIC[, .(findingId, situationId)][order(findingId)],
	file = paste0(OUTPUT, 'historic.csv')) 
fwrite(HISTORIC, file = paste0(WORKING, 'historic_check.csv'))

# Suspected

fwrite(SUSPECTED[, .(findingId, situationId)][order(findingId)],
	file = paste0(OUTPUT, 'suspected.csv')) 
fwrite(SUSPECTED, file = paste0(WORKING, 'suspected_check.csv'))

# Exclusions - concepts not to use in MedCAT training
write(as.character(exclude_conceptId),
	file = paste0(WORKING, 'exclude_concepts.txt'), ncolumns = 1)

blacklist <- intersect(c(ignore_findings, procedures), SCT2$conceptId)

# Blacklist of ignorable concepts not to present as final output 
write(as.character(sort(unique(blacklist)),
	file = paste0(OUTPUT, 'problem_blacklist.csv'), ncolumns = 1)

cat('\nDONE\n')
sink()

