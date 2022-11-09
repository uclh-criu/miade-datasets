# Script to create additional synonyms for MedCat concept database
# for disorder / finding terms.
# Using synonyms table and OXMIS table from the Freetext Matching
# Algorithm lookups

#### NEED TO INCLUDE PROCEDURES THAT HAVE 'HISTORY OF' FOR PROBLEM LIST
# TODO
# Also some other vague failure / impairment concepts to exclude 
# Remove hyphen

#### NEED TO EXCLUDE ALLERGY / INTOLERANCE CONCEPTS? ####
# Descendants of 'Propensity to adverse reactions to substance'

# Now without Read and OXMIS maps, just using manual synonyms.

# Options / steps:
# 1. Acronym expansion from SNOMED CT descriptions
# 2. Manual synonyms
# 3. Manual alternate terms (e.g. FMA alternate terms)
# 4. Alternate word order -- [I believe MedCAT does not handle this]

# How does MedCAT hande stopwords and punctuation?

# Format of output file for MedCAT:
# cui, name, ontologies, name_status

# Anoop Shah, July 2022

# Documentation: https://docs.google.com/document/d/1suIa4NoFZ92rwFzPi6juNkqWFQ-0tGtng8nxR874lJs/edit?pli=1#

# Create a script for conversion of Read to SNOMED
library(data.table)
library(Rdiagnosislist)
library(parallel)
library(CALIBERfma)

setwd('~/Terminologies/')

# Log file
sink(paste0('Generate_MiADE_SNOMED_problems_cdb_v5', Sys.time(), '.log'),
	split = TRUE)

# Import SNOMED CT dictionary (created using loadSNOMED function)
SNOMED <- readRDS('~/Terminologies/Feb2022/SNOMED_Feb2022_active_only.RDS')

# Initial set of concepts of interest for problems
SCT <- SNOMED$DESCRIPTION[semanticType(conceptId) %in%
	c('finding', 'disorder', 'event', 'situation') &
	SNOMED$DESCRIPTION$conceptId %in% SNOMED$CONCEPT$id,
	list(conceptId, term, name_status = ifelse(
	typeId == bit64::as.integer64('900000000000003001'), 'P', 'A'))]

# Suspected, Historic, Negated version of concepts (where available)
relation_source <- function(destination){
	sort(unique(as.SNOMEDconcept(SNOMED$RELATIONSHIP[
		destinationId == SNOMEDconcept(destination)]$sourceId)))
}

suspected <- relation_source('Suspected')
known_absent <- relation_source('Known absent')
known_present <- relation_source('Known present')
subject_of_record <- relation_source('Subject of record')
family_person <- relation_source('Person in family of subject')
in_the_past <- relation_source('In the past')
current_or_specified_time <- relation_source('Current or specified time')

allergy <- descendants('Propensity to adverse reactions to substance')
ignore <- SNOMEDconcept(c('Disease', 'Clinical finding', 'Problem',
	'Impairment', 'Chief complaint', 'Sign', 'Complaint', 'Sequela',
	'Early complication', 'Co-morbid conditions', 'Pre-existing condition',
	'Acute disease', 'Subacute disease', 'Chronic disease',
	'General problem AND/OR complaint', 'Evaluation finding',
	'Administrative statuses', 'Finding by site', 'Finding by method',
	'Clinical history and observation findings',
	'Adverse incident outcome categories', 'Prognosis/outlook finding',
	'General clinical state finding', 'Disorder by body site',
	'Failure', 'Acute failure', 'Subacute failure', 'Chronic failure',
	'Decompensation', 'Discrepancy', 'Idiosyncrasy', 'Inefficiency',
	'General body state finding', 'General clinical state finding'))

ASSOC_FINDING <- SNOMED$RELATIONSHIP[
	typeId == SNOMEDconcept('Associated finding'),
	.(situationId = sourceId, sitTerm = description(sourceId)$term,
	findingId = destinationId, finTerm = description(destinationId)$term)][
	situationId %in% SCT$conceptId & findingId %in% SCT$conceptId]

HISTORIC <- ASSOC_FINDING[situationId %in% intersect(known_present,
	intersect(subject_of_record, in_the_past))][order(findingId)]
# Keep if only one situation per finding, to avoid errors
HISTORIC[, Nsit := .N, by = findingId]
fwrite(HISTORIC[Nsit == 1, .(findingId, situationId)],
	file = 'MiADE/historic.csv') 
fwrite(HISTORIC, file = 'MiADE/historic_check.csv') 

NEGATED <- ASSOC_FINDING[situationId %in% intersect(known_absent,
	intersect(subject_of_record, current_or_specified_time))][order(findingId)]
# Keep if only one situation per finding, to avoid errors
NEGATED[, Nsit := .N, by = findingId]
fwrite(NEGATED[Nsit == 1, .(findingId, situationId)],
	file = 'MiADE/negated.csv') 
fwrite(NEGATED, file = 'MiADE/negated_check.csv')

SUSPECTED <- ASSOC_FINDING[situationId %in% intersect(suspected,
	intersect(subject_of_record, current_or_specified_time))][order(findingId)]
# Keep if only one situation per finding, to avoid errors
SUSPECTED[, Nsit := .N, by = findingId]
fwrite(SUSPECTED[Nsit == 1, .(findingId, situationId)],
	file = 'MiADE/suspected.csv') 
fwrite(SUSPECTED, file = 'MiADE/suspected_check.csv')

# List of suspected, negated, historic, allergy and family history concepts to
# exclude from MedCAT annotation training
exclude_conceptId <- sort(unique(c(suspected, in_the_past, known_absent,
	family_person, allergy, ignore)))
write(as.character(exclude_conceptId),
	file = 'MiADE/exclude_concepts.txt', ncolumns = 1)

#### FUNCTIONS FOR ADDING SPECIFIC SNOMED CT SYNONYMS ####

# Add acronyms for SNOMED CT disorders, findings and events
get_acronym <- function(terms){
	# returns a lower case acronym for the term if it is in the form
	# ABCD - Another Bland Cardiovascular Disease
	words <- strsplit(tolower(terms), ' ')
	is_acronym <- sapply(words, function(x) length(x) >= 3 & x[2] == '-')
	stated_acronym <- sapply(words[is_acronym], function(x) x[1])
	created_acronym <- sapply(words[is_acronym], function(x){
		paste(substr(x[3:length(x)], 1, 1), collapse = '')
	})
	acronym <- rep(NA_character_, length(terms))
	acronym[is_acronym] <- ifelse(stated_acronym == created_acronym,
		stated_acronym, NA_character_)
	acronym
}

add_rows <- function(SCT, NEW){
	message('SCT has ', nrow(SCT), ' rows and NEW has ', nrow(NEW), ' rows.')
	SCT <- rbind(SCT[, .(conceptId, term)], NEW[, .(conceptId, term)])
	SCT <- SCT[!duplicated(SCT)]
	message(nrow(SCT), ' rows after deduplication.')
	SCT
}

add_specific_acronyms<- function(SCT){
	# Specific acronyms which have a usual meaning even though it may
	# be ambiguous in SNOMED CT  
	specific_acronyms <- c(
	'af - atrial fibrillation',
	'cad - coronary artery disease',
	'ccf - congestive cardiac failure',
	'chd - coronary heart disease',
	'dm - diabetes mellitus',
	'mi - myocardial infarction',
	'ms - multiple sclerosis',
	'pid - pelvic inflammatory disease',
	'rapd - relative afferent pupillary defect',
	'rapd - relative afferent pupil defect',
	'nfa - no fixed abode',
	'ld - learning difficulties',
	'hf - heart failure',
	'psp - progressive supranuclear palsy',
	'vf - ventricular fibrillation')
	add_rows(SCT, add_acronyms(SCT[tolower(term) %in% specific_acronyms]))
}

add_acronyms <- function(SCT, outfile = 'MiADE/acronyms_to_check.csv'){
	SCT[, acronym := get_acronym(term)]

	# If the acronym is identical to an actual term, add the term as well
	# e.g. COLD = cold or 'chronic obstructive lung disease'
	# This is important for disambiguation
	all_acronyms <- sort(unique(SCT$acronym))
	all_acronyms <- all_acronyms[!is.na(all_acronyms)]
	
	# If term is already an acronym
	SCT[tolower(term) %in% all_acronyms, acronym := tolower(term)]

	ACRONYMS <- SCT[!is.na(acronym), .(conceptId,
		term = acronym, origterm = tolower(term))]
	ACRONYMS[, nconcepts := uniqueN(conceptId), by = term]
	ACRONYMS[, nterms := uniqueN(tolower(origterm)), by = term]
	if (!is.null(outfile)){
		fwrite(ACRONYMS[(nterms > 1 & nconcepts > 1)][order(term)],
			file = outfile)
	}
	SCT[, acronym := NULL]
	
	# Add acronyms where there is only one unique concept linked to the
	# acronym 
	add_rows(SCT, ACRONYMS[nconcepts == 1])
}

add_alternate_word_order <- function(SCT){
	# Allow alternative word order for terms with up to 3
	# non-numeric words
	
	reorder <- function(words, y){
	# function to reorder words and create new phrases
	# words is a list of character vectors
	# y is an integer vector representing the order in which the
	# words are to be combined 
	out <- sapply(words, function(x) x[y[1]])
		if (length(y) > 1){
			for (i in y[2:length(y)]){
				out <- paste(out, sapply(words, function(x) x[i]))
			}
		}
		out
	}

	perm <- function(v) {
		# function to generate permutations
		n <- length(v)
		if (n == 1) v
		else {
			X <- NULL
			for (i in 1:n) X <- rbind(X, cbind(v[i], perm(v[-i])))
			X
		}
	}
	
	SCT[, words := strsplit(tolower(gsub('[-\\,\\.]*', '', term)), ' ')]
	# Strip 'the', 'and', 'of', 'from', 'by', 'in', 'on'
	SCT[, words := lapply(words, function(x)
		x[!x %in% c('and', 'of', 'the', 'from', 'by', 'in', 'on')])]
	# Remove words that are entirely non-alphanumeric
	SCT[, words := lapply(words, function(x)
		x[x %like% '[[:alpha:][:digit:]]'])]	
	SCT[, numwords := sapply(words, length)]
	
	# Alternate word order for 2-4 non-numeric terms
	ALT2 <- SCT[numwords == 2 & !(term %like% '[0-9]')]
	ALT3 <- SCT[numwords == 3 & !(term %like% '[0-9]')]
	ALT4 <- SCT[numwords == 4 & !(term %like% '[0-9]')]
	
	add_rows(SCT, rbind(
		rbindlist(lapply(1:factorial(2), function(x){
			ALT2[, .(conceptId, term = reorder(words, perm(1:2)[x, ]))]
			}))
		, rbindlist(lapply(1:factorial(3), function(x){
			ALT3[, .(conceptId, term = reorder(words, perm(1:3)[x, ]))]
			}))
#		, rbindlist(lapply(1:factorial(4), function(x){
#			ALT4[, .(conceptId, term = reorder(words, perm(1:4)[x, ]))]
#			}))
	))
}

add_synonyms <- function(SCT, text, sct){
	NEW <- copy(SCT[, .(conceptId, term = tolower(term))])
	
	# Search for and replace synonyms 
	for (i in seq_along(text)){
		scti <- paste0(' ', sct[i], ' ')
		texti <- paste0(' ', text[i], ' ')
		has <- NEW[, grepl(scti, paste0(' ', NEW$term, ' '))]
		if (any(has)){
			message('Replacing ', scti, ' with ', texti,
				' in ', sum(has), ' terms.')
			SYN <- NEW[has]
			SYN[, term := sub('^ ', '', sub(' $', '',
				sub(scti, texti, paste0(' ', term, ' '))))]
			NEW <- rbind(NEW, SYN)
			NEW <- NEW[!duplicated(NEW)]
			gc()
		}
	}
	NEW <- NEW[!(term %in% tolower(SCT$term))]
	
	add_rows(SCT, NEW)
}

#### PROCESS SNOMED CONCEPTS - ADD SYNONYMS ####

message('SCT has ', nrow(SCT), ' rows (', uniqueN(SCT$conceptId),
	' concepts).')
SCT <- SCT[!(conceptId %in% exclude_conceptId)]
message('SCT has ', nrow(SCT), ' rows after removing exclusion concepts (',
	uniqueN(SCT$conceptId), ' concepts).')

# SCT now contains all finding, event and situation concepts
# (without exclusions e.g. family history)

# Lower case and remove semantic tag
SCT[, term := sub(' \\((event|disorder|situation|finding)\\)$', '',
	tolower(term))]
PREFERRED <- SCT[name_status == 'P']
SCT[, name_status := NULL]
SCT <- SCT[!duplicated(SCT)]
message('SCT has ', nrow(SCT), ' rows after removing semantic tags')

# Add specific synonyms
SCT <- add_synonyms(SCT, 'cancer', 'malignant neoplasm')

# Allow alternate word order for 2-3 words (excluding the, of and etc.)
SCT <- add_alternate_word_order(SCT)

# Add specific common acronyms
SCT <- add_specific_acronyms(SCT)

# Add any other unambiguous acronyms specified in SNOMED CT 
SCT <- add_acronyms(SCT)

# Function to quickly check for a term 
find <- function(text){
	SCT[conceptId %in% SCT[term == text]$conceptId]
}

# Export full vocab for MiADE 
setkey(SCT, conceptId)
fwrite(SCT, file = 'MiADE/full_condition_list.csv')

# Format of output file for MedCAT:
# cui, name, ontologies, name_status
SCT[, name_status := PREFERRED[SCT, on = c('conceptId', 'term')]$name_status]
SCT[is.na(name_status), name_status := 'A']
fwrite(SCT[, .(cui = paste0('S', conceptId), name = term,
	ontologies = 'SNO', name_status)], file = 'MiADE/full_condition_list_MedCAT_cdb.csv')

# Remove any synonyms that are already in SNOMED CT
# (whether in original or lower case form)
OUT <- merge(SCT[, .(conceptId, term)], rbind(SNOMED$DESCRIPTION[, .(conceptId, term = term,
	already = TRUE)], SNOMED$DESCRIPTION[, .(conceptId, term = tolower(term),
	already = TRUE)]), by = c('conceptId', 'term'), all.x = TRUE)
OUT <- OUT[is.na(already)]
OUT[, already := NULL]

# Export list of additional concepts
fwrite(OUT, file = 'MiADE/extra_condition_synonyms.csv')

cat('\nDONE\n')
sink()

###

#SNOMED error:

#Associated finding of 
#16697921000119109  Suspected amblyopia of bilateral eyes
#is the pair of:
#Amblyopia of left eye and Amblyopia of right eye
#and not
#347241000119107 | Amblyopia of bilateral eyes (disorder)
