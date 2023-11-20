# Script to create additional synonyms for MedCat concept database
# for disorder / finding terms.

# Anoop Shah

# Steps for CDB expansion
# - Processes each concept individually:
#   - Standard expansion (default) using SNOMED synonyms only
#   - Option for 'deep' expansion also using manual synonyms from
#     Freetext Matching Algorithm
# - Choice based on SNOMED CT usage in primary care 2021-2022:
#   >= 500 uses: generate extra synonyms, standard
#   >= 10000 uses or in the list of mandataory comorbidities: deep 
# 
# Algorithm:
# - Use other SNOMED concepts (finding site, associated morphology,
#   organism, attribute etc.) to find synonyms
# - Remove ignorable words (and, the, of, from)
# - Within each finding or disorder concept and split into SNOMED concept parts 
#   e.g. 'Fracture of femur' --> 'Fracture' 'Femur' --> 'Femur' 'fracture'
# - Split into sections by 'due to', 'without',  etc.
#   e.g. 'Gestational proteinuria without hypertension'
#        'Gestational proteinuria' 'WITHOUT' 'hypertension'
# - Allow reordering of the SNOMED concept parts within section. If there are
#   words between two SNOMED CT concepts keep the whole set together
# - For lateralisable body parts, optionally add 'Left' or 'Right' before
#   the term to generate a 'fake' lateralised term

######################################
# Import relevant files

# Import SNOMED CT dictionary (created using loadSNOMED function)
# Using May 2022 SNOMED dictionary and Read maps
# Import SNOMED CT dictionary (created using loadSNOMED function)
if (!exists('SNOMED')){
	SNOMED <- readRDS('~/Terminologies/SNOMED_May2022.RDS')
}

defaultfilepath <- '~/pCloudDrive/anoop_work/2019to2021/MiADE/Terminologies/'
if (!exists('WORKING')){ WORKING <- paste0(defaultfilepath, 'WORKING/') }
if (!exists('OUTPUT')){ OUTPUT <- paste0(defaultfilepath, 'OUTPUT/') }
if (!exists('REF')){ REF <- paste0(defaultfilepath, 'REF/') }
if (!exists('LOGS')){ LOGS <- paste0(defaultfilepath, 'LOGS/') }
# Log file
sink(paste0(LOGS, 'gen_synonyms_problems_v3_', Sys.time(), '.log'),
	split = TRUE)

# File name for extra synonym file
file_extras <- paste0(WORKING, 'extra_synonyms_28sep.csv')

# List of mandatory comorbidities for deep expansion
MANDATORY <- fread(paste0(REF, 'mandatory_comorbidities.csv'))

# Import SNOMED CT dictionary (created using loadSNOMED function)
if (!exists('SNOMED')){
	SNOMED <- readRDS('~/Terminologies/SNOMED_May2022.RDS')
}

# SNOMED disorders and findings
disorders <- descendants('Disorder')
findings <- descendants('Clinical finding')

# Load frequency table of SNOMED CT concepts from UK primary care data
# FREQ <- fread('https://files.digital.nhs.uk/63/992555/SNOMED_code_usage_2021_22.txt')
FREQ <- fread(paste0(REF, 'SNOMED_code_usage_2021_22.txt'))
SNOMED$CONCEPT[, freq := 
	FREQ[, .(id = SNOMED_Concept_ID, freq = as.integer(FREQ$Usage))][
	SNOMED$CONCEPT, on = 'id']$freq]
SNOMED$CONCEPT[is.na(freq), freq := 0] # 0 means less than 10

#### SNONYMS FROM MANUAL TABLE ####
# Load FMA synonyms table
# Use FMA synonyms: only if priority is >0, substitute 'text ' for
# 'read' in part of a SNOMED description, only if number of words
# in 'text' is equal to or less than the number of words in 'read',
# and avoid spaced out one letter acronyms (do not want to generate
# unnecessary / unwarranted expansions or incorrectly interpret words as acronyms)
# From https://github.com/anoopshah/freetext-matching-algorithm-lookups/blob/master/synonyms.txt
SYN <- fread('https://raw.githubusercontent.com/anoopshah/freetext-matching-algorithm-lookups/master/synonyms.txt')
SYN[, nword_text := 1 + nchar(gsub('[^ ]', '', text))]
SYN[, spaced_out_acronym := (nchar(text) + 1) / 2 == nword_text]
SYN[, nword_read := 1 + nchar(gsub('[^ ]', '', read))]
# Keep only entries with priority > 2 (i.e. similar) 
SYN <- SYN[priority > 2 & nword_text <= nword_read & spaced_out_acronym == FALSE]
setkey(SYN, read)

#### LATERALITY ####

# Create a list of lateralised body structures
body_structures <- descendants('Body structure')
left_structures <- relatedConcepts('Left', typeId = 'Laterality',
	reverse = TRUE)
right_structures <- relatedConcepts('Right', typeId = 'Laterality',
	reverse = TRUE)
lateralisable_structures <- relatedConcepts('Side', typeId = 'Laterality',
	reverse = TRUE)

# Create a list of lateralised findings
left_findings <- relatedConcepts(left_structures, typeId = 'Finding site',
	reverse = TRUE)
right_findings <- relatedConcepts(right_structures, typeId = 'Finding site',
	reverse = TRUE)
bilateral_findings <- intersect(left_findings, right_findings)
left_findings <-  setdiff(left_findings, bilateral_findings)
right_findings <-  setdiff(right_findings, bilateral_findings)
lateralisable_findings <- relatedConcepts(lateralisable_structures,
	typeId = 'Finding site', reverse = TRUE)

# Create a laterality flag
SNOMED$CONCEPT[, laterality := "No laterality"]
SNOMED$CONCEPT[id %in% lateralisable_findings, laterality := "Lateralisable"]
SNOMED$CONCEPT[id %in% bilateral_findings, laterality := 'Bilateral']
SNOMED$CONCEPT[id %in% left_findings, laterality := 'Left']
SNOMED$CONCEPT[id %in% right_findings, laterality := 'Right']
SNOMED$CONCEPT[, .N, by = laterality]

#### FUNCTIONS FOR ADDING SPECIFIC SNOMED CT SYNONYMS ####

# Options / steps for CDB expansion
# - Acronym expansion from SNOMED CT descriptions
# - Use other SNOMED concepts (finding site, associated morphology,
#   organism, attribute etc.) to find synonyms
# - Remove ignorable words (and, the, of, from, with)
# - Within each finding or disorder concept and split into SNOMED concept parts 
#   e.g. 'Fracture of femur' --> 'Fracture' 'Femur' --> 'Femur' 'fracture'
#   - single words initially
#   - then double words
# - Section by 'attribute' concepts + 'caused by' (not a synonym of 'Due to')
#   e.g. 'Gestational proteinuria without hypertension'
#        'Gestational proteinuria' 'WITHOUT' 'hypertension' (i.e. do not reorder)
# - Allow reordering of the SNOMED concept parts within section. If there are
#   words between two SNOMED CT concepts keep the whole set together
# - Reverse causality (e.g. XXX due to YYY --> YYY causing XXX)
# - Remove duplicates
# 
# - For lateralisable body parts, optionally add 'Left' or 'Right' before
#   the term to generate a 'fake' lateralised term
# - In the post processing stage, lateralised terms are converted to native
#   terms + laterality attribute

remove_stopwords <- function(terms,
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
	text <- remove_acronyms(text)
	# remove non-alphanumeric characters and multiple spaces
	text <- gsub('[^[:alnum:] ]', ' ', text)
	text <- gsub(' +', ' ', text)
	text
}

remove_acronyms <- function(terms){
	# if x is of the form abc another bland cardiac then only the
	# expansion is retained, otherwise it is unchanged. This is because
	# many of the acronyms in SNOMED CT are non-standard and confusing.
	# Useful acronyms are included at a later stage if they are unique
	# or if they have been added manually.
	# Must be in lower case already
	
	without_acronym <- terms
	words <- strsplit(tolower(terms), split = ' ')
	maybe_acronym <- sapply(words, length) >= 3 & terms %like% '^([[:alnum:]]+) - '
	stated_acronym <- strsplit(sapply(words, function(x) x[1]),
		split = character(0))
	stated_expansion <- sub('^([[:alnum:]]+) - (.*)$', '\\2', terms)
	created_acronym <- sapply(words, function(x){
		substr(x[3:length(x)], 1, 1)
	})
	
	if (any(maybe_acronym)){
		for (i in seq_along(without_acronym)[maybe_acronym == TRUE]){
			if ((length(setdiff(stated_acronym[[i]], created_acronym[[i]])) <= 2) &
				(length(setdiff(created_acronym[[i]], stated_acronym[[i]])) <= 1)){
				without_acronym[i] <- stated_expansion[i]
			}
		}
	}
	
	without_acronym
}

# Testing acronym removal
terms <- c('nof - fracture of neck of femur', 'ca - cancer',
	'nstemi - non-st segment elevation mi',
	'abc - aa bb cc', 'abc - aa bb cc dd', 'abc - aa bb',
	'abc - def ghi', 'abc , . / ; a', 'ivdu - intravenous drug user')
print(remove_acronyms(terms))
print(remove_stopwords(terms))

add_manual_synonyms <- function(x,
	search = paste0(' ', SYN$read, ' '),
	replacement = paste0(' ', SYN$text, ' ')){
	# x is a character vector of concept names to expand
	# search, replacement are character vectors of same length
	# e.g. search <- c('fever', 'femur'); replacement <- c('pyrexia', 'femoral')
	# returns x concatenated with versions in which search is
	# replaced with the corresponding replacement
	out <- mapply(function(search, replacement){
		sub(' $', '', sub('^ ', '', sub(search,
		replacement, paste0(' ', x, ' '))))
	}, search, replacement)
	unique(c(x, unlist(out)))
}

add_manual_synonyms_deep <- function(x,
	search = paste0(' ', SYN$read, ' '),
	replacement = paste0(' ', SYN$text, ' ')){
	#### NOT USED - results in excessive expansion of synonyms
	# x is a character vector of concept names to expand
	# search, replacement are character vectors of same length
	# e.g. search <- c('fever', 'femur'); replacement <- c('pyrexia', 'femoral')
	# returns x concatenated with versions in which search is
	# replaced with the corresponding replacement
	x <- paste0(' ', x, ' ')
	for (i in seq_along(search)){
		if (any(grepl(search[i], x))){
			x <- unique(c(x, sub(search[i], replacement[i], x)))
		}
	}
	sub('^ ', '', sub(' $', '', x))
}

get_synonyms <- function(conceptIds, SNOMED, limit_length = 4,
	words_only = FALSE, deep = FALSE){
	# returns a list with synonyms for each synonym of each concept
	# designed for disorders, findings and body structures
	out <- NULL
	for (i in seq_along(conceptIds)){
		x <- sort(unique(remove_stopwords(description(conceptIds[i], SNOMED = SNOMED,
			include_synonyms = TRUE)[type == 'Synonym']$term)))
		# remove synonyms longer than limit_length words
		x <- x[nchar(gsub('[a-z]', '', x)) < limit_length]
		# remove the word 'finding', 'observation' or 'structure'
		# at the end of phrase
		x <- sub('( finding| findings| structure| observation| observations)$', '', x)
		# remove the word 'general' / 'entire' at the beginning of phrase
		# e.g. 'general bone observation'
		x <- sub('^(finding |observation |general |entire |structure )', '', x)
		# add manual synonyms using SYN table (if extra snonyms required)
		if (deep) { x <- add_manual_synonyms(x) }
		# keep unique values only
		x <- unique(x)
		# create a list of synonyms for each option
		# avoid synonyms where an extra word is added, e.g.
		# 'fracture' --> 'fracture bone' because this can lead to
		# unwanted expansion; also exclude self from each set in list
		if (words_only){
			# return a list of empty vectors
			y <- lapply(x, function(j) character(0))
		} else {
			y <- lapply(x, function(j){
				x[sapply(x, function(k) {
					!(paste0(' ', k, ' ') %like% paste0(' ', j, ' ')) &
					(k %like% ' ' | j %like% ' ')
					# root word/phrase is not entirely contained within
					# synonym and root and synonym are not both single words
					# in which case the synonym table is not needed because
					# this can be matched by MedCAT's disambiguation
				})]
			})
		} 
		names(y) <- x
		if (is.null(out)){
			out <- y
		} else {
			out <- c(out, y)
		}
	}
	out
}

any_word_duplicate <- function(x){
	# whether any word is duplicated in x
	y <- strsplit(x, ' ')
	sapply(y, function(z) any(duplicated(z[z != ''])))
}

expand_concept <- function(conceptId, SNOMED, deep = FALSE){
	# Provides an enhanced list of synonyms including lateralised
	# versions if appropriate for a single concept
	
	# Restrict the maximum number of words in synonym to 8, or
	# number of words in FSN + 1 (to allow for left/right)
	conceptId <- as.SNOMEDconcept(conceptId)
	maxlength <- min(8, nchar(gsub('[^ ]', '', description(conceptId,
		SNOMED = SNOMED)$term[1])))
	
	#### Get list of synonyms ####
	
	# If lateralisable, add 'left' and 'right' options
	get_synonyms_left_right <- function(conceptId, deep = FALSE){
		finding_site <- get_synonyms(
			relatedConcepts(conceptId, 'Finding site', SNOMED = SNOMED),
			SNOMED = SNOMED, deep = deep)
		# finding_site is a synonym list
		if (SNOMED$CONCEPT[id == conceptId]$laterality == 'Lateralisable' &
			!is.null(finding_site)){
			y <- names(finding_site)
			finding_site <- lapply(y, function(z){
				c(finding_site[[z]],
				paste('left', c(z, finding_site[[z]])),
				paste('right', c(z, finding_site[[z]])))
			})
			names(finding_site) <- y
		}
		finding_site
	}
	
	# Finding site
	finding_site <- get_synonyms_left_right(conceptId, deep = deep)
	
	# Include underlying cause ('Due to' attribute)
	cause <- relatedConcepts(conceptId, 'Due to', SNOMED = SNOMED)
	
	# Finding site (including finding site of 'Due to' linked concepts)
	if (length(cause) > 0){
		for (i in seq_along(cause)){
			finding_site <- c(finding_site,
				get_synonyms_left_right(cause[i], deep = deep))
		}
	}
	
	# Putative synonyms from causative agent
	# Include all ancestors for self and cause
	# Related concepts for a disorder
	related <- c(intersect(findings, ancestors(c(cause, conceptId),
		SNOMED = SNOMED)),
		relatedConcepts(union(cause, conceptId),
		'Causative agent', SNOMED = SNOMED))
	
	# Include morphology but only for words in SNOMED, not synonyms
	morph <- relatedConcepts(union(cause, conceptId),
		'Associated morphology', SNOMED = SNOMED)
	
	# Synonyms of self and parents
	synonyms <- c(finding_site,
		get_synonyms(related, SNOMED = SNOMED, deep = deep),
		get_synonyms(morph, SNOMED = SNOMED, words_only = TRUE))
	
	#### Process concept ####
		
	S <- data.table(orig = names(get_synonyms(conceptId,
		SNOMED = SNOMED, limit_length = maxlength, deep = deep)),
		prelink = NA_character_, link = NA_character_,
		postlink = NA_character_)
	
	links <- c('due to', 'caused by', 'because of',
		'causing', 'resulting in', 'leading to',
		'without')
	
	# Create a data.table with components of synonyms:
	for (i in links){
		S[orig %like% i, prelink := sub(
			paste0('^(.*) ', i, ' (.*)$'), ' \\1 ', orig)]
		S[orig %like% i, link := i]
		S[orig %like% i, postlink := sub(
			paste0('^(.*) ', i, ' (.*)$'), ' \\2 ', orig)]
		S[is.na(prelink), prelink := paste0(' ', orig, ' ')]
	}
	S[, orig := NULL]
	
	# Reversal of phrases (e.g. swollen arm --> arm swollen) 
	regex <- paste0('^ (', paste(names(synonyms), collapse = '|'),
		') (', paste(names(synonyms), collapse = '|'), ') $')
	MATCH <- S[prelink %like% regex]
	if (nrow(MATCH) > 0){
		# reverse the phrases
			S <- rbind(S, MATCH[, .(
				prelink = sub(regex, ' \\2 \\1 ', prelink), link, postlink)])
		}
	MATCH <- S[postlink %like% regex]
	if (nrow(MATCH) > 0){
		# reverse the phrases
			S <- rbind(S, MATCH[, .(
				prelink, link, postlink = sub(regex, ' \\2 \\1 ', postlink))])
		}
	S <- S[!duplicated(S)]
	
	# Now expand the synonyms in table S
	if (length(synonyms) > 0){
		for (i in 1:length(synonyms)){
			# prelink
			MATCH <- S[prelink %like% paste0(' ', names(synonyms)[i], ' ')]
			if (nrow(MATCH) > 0 & length(synonyms[[i]]) > 0){
				for (j in 1:length(synonyms[[i]])){
					NEW <- MATCH[, .(
						prelink = sub(paste0(' ', names(synonyms)[i], ' '),
						paste0(' ', synonyms[[i]][j], ' '), prelink),
						link, postlink)]
					S <- rbind(S, NEW)
				}
			}
			S <- S[!duplicated(S)]
			
			# postlink
			MATCH <- S[postlink %like% paste0(' ', names(synonyms)[i], ' ')]
			if (nrow(MATCH) > 0 & length(synonyms[[i]]) > 0){
				for (j in 1:length(synonyms[[i]])){
					NEW <- MATCH[, .(prelink, link,
						postlink = sub(paste0(' ', names(synonyms)[i], ' '),
						paste0(' ', synonyms[[i]][j], ' '), postlink))]
					S <- rbind(S, NEW)
				}
			}
			S <- S[!duplicated(S)]
		}
	}
	
	# Reinstate entire phrase
	S[, final := gsub('  ', ' ', gsub('^ | $| NA|NA', '',
		paste0(prelink, link, postlink)))]
		
	# Remove any with duplicate words
	S <- S[!(any_word_duplicate(final))]
	
	# Remove any with left/right duplication
	S <- S[!(final %like% 'left' & final %like% 'right')]
	
	# Check if the original term includes 'finding' 'observation',
	#   'disorder' 'disease' etc., and exclude any terms that don't
	#   (e.g. to prevent 'cardiac finding' being converted to 'cardiac')
	if (all(tolower(description(conceptId, SNOMED = SNOMED,
		include_synonyms = TRUE)$term) %like%
		'finding|observation|sympyom|sign')){
		S <- S[final %like% 'finding|observation|disorder|disease|illness|sign|symptom|problem']
	}
	
	# Remove any that are longer than maxlength
	S <- S[nchar(gsub('[^ ]', '', final)) + 1 <= maxlength]

	unique(S$final)
}

#### TESTING ####

test <- function(concept){
	conceptId <- SNOMEDconcept(concept)
	cat('\nTEST: Original synonyms:\n')
	out1 <- description(conceptId, include_synonyms = TRUE)$term
	print(out1)
	
	cat('\nAdditional synonyms (deep = FALSE):\n')
	out2 <- expand_concept(conceptId, SNOMED, deep = FALSE)
	print(setdiff(out2, tolower(out1)))
	
	cat('\nAdditional synonyms (deep = TRUE):\n')
	out3 <- expand_concept(conceptId, SNOMED, deep = TRUE)
	print(setdiff(out3, tolower(out2)))
	invisible(out3)
}

test('Cardiac finding')

test('Hepatitis A')

test('Fracture of radius')

test('Swelling of upper limb')

test('5913000')
# 5913000 |Fracture of neck of femur (disorder)|

test('704330001')
# 704330001 | Osteoporotic fracture of femur (disorder)

test('16000511000119103')
# 16000511000119103 |Cerebrovascular accident due to occlusion of left middle cerebral artery (disorder)|

test('10625711000119105')
# Bronchopneumonia caused by Streptococcus pneumoniae (disorder)|

test('765330003')
# Autosomal dominant polycystic kidney disease

test('449618007')
# 449618007 |Swelling of upper limb (finding)|

test('Mild depression')

test('838450006') # Severe mitral valve stenosis (disorder)|

test('61142002') # Microphthalmos

h <- test('363393007') # Tonsil cancer
H <- data.table(h, nwords = 1 + nchar(gsub('[a-z]', '', h)))
print(H[order(-nwords)])

test('21719001') # Hay fever - has multiple causes

#### BATCH PROCESS EXTRA SYNONYMS ####

# To expand concepts used over 500 times in the UK or those in the
# mandatory list
TOEXPAND <- data.table(id = union(findings, MANDATORY$conceptId))
TOEXPAND[, freq := SNOMED$CONCEPT[TOEXPAND, on = 'id']$freq]
TOEXPAND[, mandatory := id %in% MANDATORY$conceptId]
TOEXPAND <- TOEXPAND[freq >= 500 | mandatory][order(-freq)]
TOEXPAND[, todo := TRUE]

# Find out which are remaining to do (in case need to restart)
if (!(file.exists(file_extras))){
	fwrite(data.table(conceptId = bit64::integer64(0),
		synonym = character(0)), file = file_extras)
} else {
	EXISTING <- fread(file_extras)
	TOEXPAND[id %in% EXISTING$conceptId, todo := FALSE]
}

# Function to process a set of concepts
process_concepts <- function(conceptIds, deep = FALSE){
	message(paste('Processing', length(conceptIds),
		'concepts with option deep =', deep))
	for (i in seq_along(conceptIds)){
		synonyms <- NULL
		try(synonyms <- expand_concept(conceptIds[i],
			SNOMED = SNOMED, deep = deep))
		if (is.null(synonyms)){
			fwrite(data.table(conceptId = conceptIds[i],
			synonym = NA_character_),
			file = file_extras, append = TRUE)
		} else {
			fwrite(data.table(conceptId = conceptIds[i],
			synonym = synonyms),
			file = file_extras, append = TRUE)
		}
		if (i %% 20 == 0){message(paste('Done', i, 'at', Sys.time()))}
	}
}

# Deep processing for mandatory comorbidities
cat('\nMandatory comorbidities\n')
process_concepts(TOEXPAND[mandatory == TRUE & todo == TRUE]$id, deep = TRUE)
TOEXPAND[id %in% fread(file_extras)$conceptId, todo := FALSE]

# Deep processing for findings with >= 10000 uses
cat('\nFindings used 10000 times or more\n')
process_concepts(TOEXPAND[freq >= 10000 & todo == TRUE]$id, deep = TRUE)
TOEXPAND[id %in% fread(file_extras)$conceptId, todo := FALSE]

# Standard processing for remaining findings with >= 500 uses
cat('\nFindings used 500 times or more\n')
process_concepts(TOEXPAND[freq >= 500 & freq < 10000 & todo == TRUE]$id)
TOEXPAND[id %in% fread(file_extras)$conceptId, todo := FALSE]

cat('\nDONE\n')
sink()

