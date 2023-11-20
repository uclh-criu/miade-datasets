# Script to create fake MetaCat training data for meds and allergies
# Anoop Shah, August 2023

library(data.table)
library(Rdiagnosislist)

# pre-specify number of reps (N)
# N <- 5
if (!exists('N')) {N <- 5}

defaultfilepath <- '~/pCloudDrive/anoop_work/2019to2021/MiADE/Terminologies/'
if (!exists('TEMPLATES')){ TEMPLATES <- paste0(defaultfilepath, 'TEMPLATES/') }
if (!exists('WORKING')){ WORKING <- paste0(defaultfilepath, 'WORKING/') }
if (!exists('OUTPUT')){ OUTPUT <- paste0(defaultfilepath, 'OUTPUT/') }
if (!exists('REF')){ REF <- paste0(defaultfilepath, 'REF/') }
if (!exists('LOGS')){ LOGS <- paste0(defaultfilepath, 'LOGS/') }

# Log file
sink(paste0(LOGS, 'gen_metacat_v6_', Sys.time(), '.log'),
	split = TRUE)

# Import SNOMED CT dictionary (created using loadSNOMED function)
if (!exists('SNOMED')){
	SNOMED <- readRDS('~/Terminologies/SNOMED_May2022.RDS')
}


# Frequency distribution of SNOMED CT concepts
FREQ <- fread(paste0(REF, 'SNOMED_code_usage_2021_22.txt'))
FREQ[, conceptId := SNOMED_Concept_ID]
suppressWarnings(FREQ[, freq := as.numeric(Usage)]) # NAs if Usage is <10
FREQ[is.na(freq), freq := 0]

# Import working CDB for medications and allergies 
CDB <- fread(paste0(WORKING, 'med_working_table_full.csv'))
CDB[, freq := FREQ[CDB, on = 'conceptId']$freq]
CDB[, numwords := sapply(strsplit(term, ' '), length)]

#### USEFUL FUNCTIONS ####

textclean <- function(text, keepgroup = NULL, keepas = NULL){
	# Text cleaning regex
	# keepgroup is a character vector of items in curly brackets to keep
	# e.g. p, r, m
	# keepas allows converstion e.g. r to p
	
	if (is.null(keepas)){
		keepas <- keepgroup
	}
	try(setattr(keepas, 'names', keepgroup))

	# Remove all punctuation except full stops, question marks, commas,
	# hyphen and line breaks
	if (!is.null(keepgroup)){
		for (i in keepgroup){
			text <- sub(paste0('\\{', i, ' ([^\\}]+)\\}'),
				paste0('ZZZKEEP', i, ' \\1KEEPZZZ'), text)
		}
	}
	# Remove {x yyyy}
	text <- gsub('\\{[[:alpha:]] ([^\\}]+)\\}', '\\1', text)
	text <- gsub('[^[:alnum:] \\.,?\n\\-]', ' ', text)
	# Replace brackets
	if (!is.null(keepgroup)){
		for (i in keepgroup){
			text <- gsub(paste0('ZZZKEEP', i), paste0('{', keepas[i]), text)
		}
	}
	text <- gsub('KEEPZZZ', '\\}', text)
	
	# Replace multiple spaces with a single normal space
	# python expression r'(?:(?!\n)\s)+'
	text = gsub(' +', ' ', text)
	
	# Remove spaces if the entire line (between two line breaks) is just spaces
	# python expression r'(?<=\n)\s+(?=\n)'
	text = gsub('\n +\n', '\n\n', text)
	
	return(text)
} 

# Test textclean
textclean('ab
cd') # should preserve line breaks
textclean('a {m another} ,./^&*?} h {r another} {n another}')
# a {x b} c -> a b c
textclean('a {m another} } h {r another} {n another}', keepgroup = c('m', 'n'))
textclean('a {m another} } h {r another} {n another}', keepgroup = c('m', 'n'), keepas = c('k', 'l'))

consam <- function(concepts, n = overallN){
	# Generate a vector of concepts of length n with random selection from
	# a SNOMEDconcept object
	sample(concepts, n, replace = TRUE)
}

desc <- function(concepts){
	# simple descriptions without semantic type suffix
	# randomly choose a synonym, preferring shorter synonyms
	# ensure only one per concept
	if (!bit64::is.integer64(concepts)){concepts <- as.SNOMEDconcept(concepts)}
	OUT <- SNOMED$DESCRIPTION[conceptId %in% concepts, .(conceptId, term)]
	OUT[, term := tolower(sub(' \\([[:alnum:]\\/\\+ ]+\\)$', '', term))]
	OUT[, term := gsub('\\[[xvd]\\]|^product containing |^product containing only ',
		'', term)]
	OUT[, term := sub('^precisely |\\[|\\]', '', term)]
	OUT[, priority := runif(nrow(OUT)) * countwords(term)]
	setkey(OUT, conceptId, priority)
	OUT[, keep := c(TRUE, rep(FALSE, .N - 1)), by = conceptId]
	OUT[keep == TRUE]$term
}

limit_numwords <- function(conceptIds, numwords){
	# Limits a set of concepts to those with fewer or equal to numwords
	conceptIds <- unique(conceptIds)
	conceptIds[countwords(desc(conceptIds)) <= numwords]
}

countwords <- function(x, sep = ' '){
	sapply(strsplit(x, sep), length)
}

#### OBTAIN RELEVANT CONCEPTS ####

allergy_concepts <- FREQ[Description %like%
	'Allergy|allergy|Adverse reaction|adverse reaction']$conceptId

SUBSTANCE_FREQ <- SNOMED$RELATIONSHIP[
	typeId == SNOMEDconcept('Causative agent'),
	.(adrId = sourceId, conceptId = destinationId)][adrId %in% allergy_concepts]
SUBSTANCE_FREQ <- merge(SUBSTANCE_FREQ,
	FREQ[, .(adrId = conceptId, freq)], on = 'adrId')
DESC <- description(SUBSTANCE_FREQ$conceptId)[, .(desc = term[1]), by = conceptId]
SUBSTANCE_FREQ[, Description := DESC[SUBSTANCE_FREQ, on = 'conceptId']$desc]
SUBSTANCE_FREQ <- SUBSTANCE_FREQ[, .(freq = sum(freq)), by = .(conceptId, Description)]
CDB[is.na(freq), freq := SUBSTANCE_FREQ[CDB[is.na(freq)], on = 'conceptId']$freq]
CDB[is.na(freq), freq := 0]

# Conversion between VTMs and substances - use exact text matching
TEMP <- merge(SUBSTANCE_FREQ, description(SUBSTANCE_FREQ$conceptId,
	include_synonyms = TRUE), by = 'conceptId')[, .(term, freq)]
CDB[, freqterm := TEMP[CDB, on = 'term']$freq]
CDB[freq == 0 & !is.na(freqterm), freq := freqterm]
TEMP <- merge(FREQ[], description(FREQ$conceptId,
	include_synonyms = TRUE), by = 'conceptId')[, .(term, freq)]
TEMP <- TEMP[, .(maxfreq = max(freq)), by = term]
CDB[, freqterm := TEMP[CDB, on = 'term']$maxfreq]
CDB[, freqterm := NULL]

# Keep the shortest term per concept ID
setkey(CDB, conceptId, numwords)
CDB[, keep := c(TRUE, rep(FALSE, .N - 1)), by  = conceptId]
CDB <- CDB[keep == TRUE]

# Short reactions, allergies, drugs, disorders etc. to use
TEMP <- CDB[category == 'reaction' & numwords < 3 & freq > 0]
reactions <- as.SNOMEDconcept(TEMP$conceptId)
setattr(reactions, 'weights', log(TEMP$freq))

TEMP <- CDB[category  %in% c('vtm', 'substance') & numwords < 3 & freq > 0]
substances <- as.SNOMEDconcept(TEMP$conceptId)
setattr(substances, 'weights', log(TEMP$freq))

TEMP <- CDB[category  %in% c('vtm', 'vmp') & numwords < 6 & freq > 0]
drugs <- as.SNOMEDconcept(TEMP$conceptId)
setattr(drugs, 'weights', log(TEMP$freq))

TEMP <- FREQ[freq > 0 & conceptId %in% descendants('Disorder')]
TEMP <- TEMP[conceptId %in% limit_numwords(TEMP$conceptId, 5)]
disorders <- as.SNOMEDconcept(TEMP$conceptId)
setattr(disorders, 'weights', log(TEMP$freq))

#### CREATE SAMPLE TEXTS FOR COMMON PHRASES

#food_text <- desc(intersect(substances, descendants('Food')))
relation_text <- desc(limit_numwords(descendants('Person in the family'), 1))
situation_text <- desc(FREQ[Description %like% '\\(situation\\)$' &
	!(Description %like% 'History|Suspected|Family') &
	freq > 100]$conceptId)
advice_text <- desc(FREQ[Description %like% '\\(situation\\)$' &
	(Description %like% 'Advise|Advice|advise|advice|Recommend|recommend') &
	freq > 10]$conceptId)
exam_text <- sub('^on examination[ \\-]*', '',
	desc(FREQ[Description %like% '^On examination' & 
	freq > 100]$conceptId))
situation_text <- setdiff(situation_text, exam_text)

#COMMON <- data.table(codename = c('FOOD', 'RELATION', 'SITUATION', 'EXAM', 'ADVICE'),
#	text = paste0('[', sapply(list(food_text, relation_text, situation_text, exam_text),
#	paste, collapse = '|'), ']'))

COMMON <- data.table(codename = c('RELATION', 'SITUATION', 'EXAM', 'ADVICE'),
	text = paste0('[', sapply(list(relation_text, situation_text, exam_text,
	advice_text), paste, collapse = '|'), ']'))

# Export to template in same format as other sample phrases
fwrite(COMMON, paste0(TEMPLATES, 'common_template.csv'))

#### CREATE ARTIFICIAL TRAINING DATA FOR CONTEXT ####

# Template format

# Paragraph template

# text = contains text with repeating phrases encased in opening "{0-4 "
#   and closing "}" with numbers (e.g. 0-4) specifying that this is to
#   repeat a random number of times selected from a uniform distribution
#   with specified limits. Can also contain codenames in capitals e.g.
#   RELATION, DRUG, which may have their own definition or may be defined
#   in a data template
# n = number of repetitions to create. Can be mutiplied by 'reps' when
#   calling the function

# Data template

# Can load multiple data templates

# codename (all capitals e.g. RELATION)
# text = contains text with mapped phrase encased in opening "{X " and closing "}"
#        (e.g. he had a {p heart attack}).
#   {p ...} = problem
#   {m ...} = medication or allergy substance
#   {r ...} = reaction
# m = medication SNOMED CT conceptId|description
# r = reaction SNOMED CT conceptId|description
# p = problem SNOMED CT conceptId|description
# p_meta_relevance = Present | Historic | Irrelevant
# p_meta_confirmed = Confirmed | Suspected | Negated
# p_meta_laterality = No laterality | Left | Right | Bilateral
# m_meta_category
# m_meta_allergytype
# m_meta_severity
# r_meta_reactionpos

# To translate text phrases in the form of pseudo-regexes into options
# [] = group a set of options
# | = separator for options
# DRUG, REACTION, SUBSTANCE, DISORDER = random choice from dictionary
# EXAM, SITUATION, RELATION, ADVICE
#
# e.g.
# [word|letter][ option |]MEDICATION fgdfg
# word {m paracetamol} fgdfg
# letter {m clopidogrel} fgdfg
# word option {m aspirin} fgdfg

# Line breaks should use the newline character (\n) in CSV documents
# (visible as new lines in a spreadsheet, shown as \n in R string literals,
# and rendered as a new line by cat)

##### TEMPLATE FORMAT

# Guidance on creation of templates

# Paragraph template must refer to data-containing elements and should
# contain the punctuation or line breaks etc. between elements.
# Any elements nested within a data template will not be selected.

# No line breaks in data template (all punctuation in paragraph template)
# Codename format:
# ...LINE = intended to be a line in a list (e.g. PLANLINE)
# ...HEADER = section header for problem list, meds list etc.

#### OUTPUT PATTERN FILE FORMAT ####

# MED / ALLERGY:

# text = contains text with mapped phrase encased in opening "{X " and closing "}"
#        (e.g. he had a {p heart attack}).
#   {m ...} = medication or allergy substance
#   {r ...} = reaction
# m = medication SNOMED CT conceptId|description
# r = reaction SNOMED CT conceptId|description
# m_meta_category = Adverse reaction | Taking | Irrelevant
# m_meta_allergytype = Allergy | Intolerance | Unspecified
# m_meta_severity = Mild | Moderate | Severe | Unspecified
# r_meta_reactionpos = Not a reaction | Before | After

# PROBLEMS

# text = contains text with mapped phrase encased in opening "{X " and closing "}"
#        (e.g. he had a {p heart attack}).
#   {p ...} = problem
# p = problem SNOMED CT conceptId|description
# p_meta_relevance = Present | Historic | Irrelevant
# p_meta_confirmed = Confirmed | Suspected | Negated
# p_meta_laterality = No laterality | Left | Right | Bilateral

#### RANDOM SELECTION OF CONCEPTS

generate_from_template <- function(paragraph_template_filepath,
	data_template_filepaths, output_medallerg, output_problems, reps = 1){
	# Processes a paragraph template
	
	if (!exists('paragraph_template_filepath')) {paragraph_template_filepath <-
		paste0(TEMPLATES, 'test_template_para.csv')}
	if (!exists('data_template_filepaths')) {data_template_filepaths <- 
		paste0(TEMPLATES, c('common_template.csv', 'test_template.csv'))}
	
	PARA_TEMPLATE <- fread(paragraph_template_filepath, sep = ',')
	DATA_TEMPLATE <- rbindlist(lapply(data_template_filepaths, fread),
		fill = TRUE)
		
	if (!exists('reps')) {reps <- 1}
	totalrows <- sum(PARA_TEMPLATE$n) * reps
	message(paste('Total number of rows to be created =', totalrows))
	
	DATA_TEMPLATE[, is_data := FALSE]
	if ('p' %in% names(DATA_TEMPLATE)){
		DATA_TEMPLATE[gsub(' ', '', p) == '', p := NA]
		DATA_TEMPLATE[, is_data := is_data | !is.na(p)]
	}
	if ('m' %in% names(DATA_TEMPLATE)){
		DATA_TEMPLATE[gsub(' ', '', m) == '', m := NA]
		DATA_TEMPLATE[, is_data := is_data | !is.na(m)]
	}
	if ('r' %in% names(DATA_TEMPLATE)){
		DATA_TEMPLATE[gsub(' ', '', r) == '', r := NA]
		DATA_TEMPLATE[, is_data := is_data | !is.na(r)]
	}
	datafields <- unique(DATA_TEMPLATE[is_data == TRUE]$codename)
	
	# Check that all codenames are included in templates
	observed_codenames <- setdiff(unique(unlist(strsplit(gsub('[^A-Z]', ' ',
		c(DATA_TEMPLATE$text, PARA_TEMPLATE$text)), ' '))), '')
	codenames <- c(DATA_TEMPLATE$codename, 'DISORDER', 'REACTION',
		'SUBSTANCE', 'DRUG', 'OTHERDISORDER', 'OTHERREACTION',
		'OTHERSUBSTANCE', 'OTHERDRUG', 'INT')
	message('The following codenames do not have an entry in templates:')
	print(observed_codenames[which(!(observed_codenames %in% codenames))])
	
	# Prepare output training data file
	TRAIN <- data.table(text = character(totalrows),
		p = character(totalrows),
		m = character(totalrows),
		r = character(totalrows),
		p_meta_relevance = character(totalrows),
		p_meta_confirmed = character(totalrows),
		p_meta_laterality = character(totalrows),
		m_meta_category = character(totalrows),
		m_meta_allergytype = character(totalrows),
		m_meta_severity = character(totalrows),
		r_meta_reactionpos = character(totalrows)
	)
		
	# Output headers
	fwrite(TRAIN[0, .(text, m, r, m_meta_category,
		m_meta_allergytype, m_meta_severity, r_meta_reactionpos)],
		file = output_medallerg)
	fwrite(TRAIN[0, .(text, p, p_meta_relevance, p_meta_confirmed,
		p_meta_laterality)], file = output_problems)
	
	process_token <- function(text, extract_data = FALSE){
		# Process a token. text is either a codename or text, in which
		# case it is split. codename refers to a name in DATA_TEMPLATE
		whichdata <- which(DATA_TEMPLATE$codename == text)
		
		if (length(whichdata) == 1){
			dtrow <- whichdata
			text <- DATA_TEMPLATE[dtrow]$text
		} else if (length(whichdata) > 1) {
			dtrow <- sample(x = whichdata, size = 1)
			text <- DATA_TEMPLATE[dtrow]$text
		} else {
			# default: text remains as is
			dtrow <- NULL
		}
		split <- unlist(strsplit(text, '\\]|\\['))
		out <- list()
		
		if (length(text) == 0){
			out$text <- ''
			return(out)
		}
		if (is.na(text) | text == ''){
			out$text <- ''
			return(out)
		}
		
		if (extract_data){
			for (colname_ in setdiff(names(TRAIN), 'text')){
				if (colname_ %in% names(DATA_TEMPLATE)){
					out[[colname_]] <- as.character(DATA_TEMPLATE[dtrow][[colname_]])
				}
			}
			if (any(split == 'SUBSTANCE')){
				theconcept <- sample(substances, 1,
					prob = attr(substances, 'weight'))
				split[split == 'SUBSTANCE'] <- paste0('{m ', desc(theconcept), '}')
				out$m <- showdesc(theconcept)
			}
			if (any(split == 'DRUG')){
				theconcept <- sample(drugs, 1, prob = attr(drugs, 'weight'))
				split[split == 'DRUG'] <- paste0('{m ', desc(theconcept), '}')
				out$m <- showdesc(theconcept)
			}
			if (any(split == 'REACTION')){
				theconcept <- sample(reactions, 1, prob = attr(reactions, 'weight'))
				split[split == 'REACTION'] <- paste0('{r ', desc(theconcept), '}')
				out$r <- showdesc(theconcept)
			}
			if (any(split == 'DISORDER')){
				theconcept <- sample(disorders, 1, prob = attr(disorders, 'weight'))
				split[split == 'DISORDER'] <- paste0('{p ', desc(theconcept), '}')
				out$p <- showdesc(theconcept)
			}
		} else {
			# Remove curly brackets around {m ...}, {r ...} etc.
			# as not extracting these data
			split <- sub('^([^\\{]*)\\{[mrp] ([^\\}]+)\\}(.*)$', '\\1\\2\\3', split)
		}
		# Data items which are not the focus in this template
		for (i in which(split %in% c('OTHERSUBSTANCE', 'SUBSTANCE'))){
			split[i] <- desc(sample(substances, size = 1,
				prob = attr(substances, 'weight')))
		}
		for (i in which(split %in% c('OTHERDRUG', 'DRUG'))){
			split[i] <- desc(sample(drugs, 1, prob = attr(drugs, 'weight')))
		}
		for (i in which(split %in% c('OTHERREACTION', 'REACTION'))){
			split[i] <- desc(sample(reactions, 1, prob = attr(reactions, 'weight')))
		}
		for (i in which(split %in% c('OTHERDISORDER', 'DISORDER'))){
			split[i] <- desc(sample(disorders, 1, prob = attr(disorders, 'weight')))
		}
		# Numbers e.g. INT_20_50
		for (i in which(split %like% '^INT_')){
			limits <- as.integer(strsplit(split[i], '_')[[1]][2:3])
			split[i] <- as.character(floor(runif(1) *
				(max(limits) + 1 - min(limits)) + min(limits)))
		}
		# Options
		for (i in which(split %like% '\\|')){
			split[i] <- sample(unlist(strsplit(split[i], '\\|')), size = 1)
		}
		# Remaining codenames are capitalised
		for (i in which(split %like% '[A-Z]$')){
			# Recursion to expand any nested data elements
			if (split[i] %in% DATA_TEMPLATE$codename){
				split[i] <- process_token(split[i], extract_data = FALSE)$text
			} else {
				warning('Codename ', split[i], ' not found in template')
			}
		}
		# Reassemble text into a single string
		out$text <- paste(split, collapse = ' ')
		out
	}

	showdesc <- function(x){
		# Show SNOMED CT concept with description
		paste0(x, " | ", description(x, SNOMED = SNOMED)$term)
	}
	
	makereps <- function(text, n, sep = ' '){
		# Repeat text n times, and return as a string
		mapply(function(text, n){
			if (is.na(n) | n == 0){
				return('')
			} else {
				paste(rep(text, n), collapse = ' ')
			}
		}, text, n)
	}
	
	expand_numbers <- function(x){
		# Convert {a-b XXX} to a set of repeats XXX XXX XXX where the number
		# of repeats is a random number between a and b inclusive
		# Use this function to create problem lists or medication lists
		# with a random number of elements.
		x <- gsub('\\{([0-9]+\\-[0-9]+ [^\\}]+)\\}', '{REP}{REP\\1REP}{REP}', x)
		x <- strsplit(x, '\\{REP\\}')[[1]]
		suppressWarnings({
			lower <- as.numeric(sub('^\\{REP([0-9]+)\\-[0-9]+ .*REP\\}$', '\\1', x));
			upper <- as.numeric(sub('^\\{REP[0-9]+\\-([0-9]+) .*REP\\}$', '\\1', x));
			random <- floor(runif(length(lower)) * (upper - lower + 1)) + lower
		})
		random[is.na(random)] <- 1
		x <- sub('^\\{REP[0-9]+\\-[0-9]+ (.*)REP\\}$', '\\1', x)
		paste(makereps(x, random), collapse = ' ')
	}
	
	trim <- function(x){
		# Function to remove leading and trailing spaces, remove
		# spaces around line breaks, and convert multiple spaces to one
		gsub(' ,', ',', gsub('^ | $', '', gsub(' \\n|\\n | \\n ', '\n',
			gsub(' +', ' ', x))))
	}
	
	k <- 1
	lastwritten <- 0
	# Loop through paragraph template rows
	for (i in 1:nrow(PARA_TEMPLATE)){
		# Original text is the paragraph template
		thetext <- gsub('\\n', ' \n ', PARA_TEMPLATE[i]$text)
		nreps <- reps * PARA_TEMPLATE[i]$n
		message('Starting para ', i, ' (', nreps, ' reps).')
		if (nreps >= 1){
			for (j in 1:nreps){
				
				# First expand numbers {1-4 xxx}
				# e.g. x <- 'Diagnoses: {2-3 DISEASE [presnt|present]} \n'
				x <- expand_numbers(thetext)
				
				# Split into a vector of tokens
				x <- trim(unlist(strsplit(x, '\\]|\\[')))
				
				# Find out if any of the tokens are a codename for a data
				# template row containing data. If there are multiple
				# possibilities, choose one at random.
				if (any(x %in% datafields)){
					chosen_data_item <- sample(which(x %in% datafields), 1)
				}
				
				# The chosen token (codename) will be the central phrase for
				# this training text. Extract data and transfer data elements
				thedata <- process_token(x[chosen_data_item], extract_data = TRUE)
				x[chosen_data_item] <- thedata$text
				for (colname_ in setdiff(names(TRAIN), 'text')){
					if (colname_ %in% names(thedata) &
						length(thedata[[colname_]]) == 1){
						TRAIN[k, (colname_) := thedata[[colname_]]]
					}
				}
				
				# Process other (non-data) tokens
				for (m in (setdiff(seq_along(x), chosen_data_item))){
					x[m] <- process_token(x[m])$text
				}
				
				# Paste the text back together
				TRAIN[k, text := trim(paste(x, collapse = ' '))]
				k <- k + 1
			}
		}
		
		TRAIN[p == '', p := NA]
		TRAIN[m == '', m := NA]
		TRAIN[r == '', r := NA]

		# Incremental write out (so that partial files are available if the
		# connection is lost
		# Medications and allergies
		fwrite(TRAIN[(lastwritten + 1):k][!is.na(m) | !is.na(r),
			.(text = textclean(text, keepgroup = c('m', 'r')), m, r,
			m_meta_category, m_meta_allergytype,
			m_meta_severity, r_meta_reactionpos)],
			file = output_medallerg, append = TRUE, col.names = FALSE)
		# Export reactions as irrelevant problems, and remove medications
		fwrite(TRAIN[(lastwritten + 1):k][!is.na(r) &
			r_meta_reactionpos %in% c('Before', 'After'),
			.(text = textclean(text, keepgroup = 'r', keepas = 'p'), p = r,
			p_meta_relevance = 'Irrelevant',
			p_meta_confirmed = 'Confirmed', p_meta_laterality)],
			file = output_problems, append = TRUE, col.names = FALSE)
		# Problems
		fwrite(TRAIN[(lastwritten + 1):k][!is.na(p),
			.(text = textclean(text, keepgroup = 'p'), p, p_meta_relevance,
			p_meta_confirmed, p_meta_laterality)],
			file = output_problems, append = TRUE, col.names = FALSE)
		lastwritten <- k
	}
	
	# Report distribution of outputs
	cat('\nProblems:\n')
	print(TRAIN[!is.na(p), .N, by = p_meta_relevance]); cat('\n')
	print(TRAIN[!is.na(p), .N, by = p_meta_confirmed]); cat('\n')
	print(TRAIN[!is.na(p), .N, by = p_meta_laterality]); cat('\n')
	cat('\nMedallergy:\n')
	print(TRAIN[!is.na(m) | !is.na(r), .N, by = m_meta_category]); cat('\n')
	print(TRAIN[!is.na(m) | !is.na(r), .N, by = m_meta_allergytype]); cat('\n')
	print(TRAIN[!is.na(m) | !is.na(r), .N, by = m_meta_severity]); cat('\n')
	print(TRAIN[!is.na(m) | !is.na(r), .N, by = r_meta_reactionpos]); cat('\n')
	
	invisible(TRAIN)
}

#### TEST CONVERSION

paragraph_template_filepath <- paste0(TEMPLATES, 'test_template_para.csv')
data_template_filepaths <- paste0(TEMPLATES, c('common_template.csv',
	'test_template_data.csv'))
output_medallerg <- paste0(TEMPLATES, 'test_pattern_medallerg.csv')
output_problems <- paste0(TEMPLATES, 'test_pattern_problems.csv')
reps <- 1

OUT <- generate_from_template(paragraph_template_filepath,
	data_template_filepaths, output_medallerg, output_problems, reps)

message('Done test')

#### ACTUAL CONVERSION

# use N=500-1000 reps for actual conversion
message('\nDoing actual conversion with number of reps = ', N)
timecode <- format(Sys.time(), '%d.%b.%H.%M.%S')
message('Saving using timecode: ', timecode)
gc()

# Creating fake notes similar to real notes
OUT <- generate_from_template(
	paragraph_template_filepath = paste0(TEMPLATES, 'template_para_v1.csv'),
	data_template_filepaths = paste0(TEMPLATES, c('common_template.csv',
		'template_medallerg_v1.csv', 'template_problems_v1.csv',
		'template_headers_v1.csv')),
	output_medallerg <- paste0(OUTPUT, 'patterns_medallerg', timecode, '.csv'),
	output_problems <- paste0(OUTPUT, 'patterns_problems', timecode, '.csv'),
	reps = N)

#### ADD EXAMPLE DATA AND MIX UP

A <- fread(paste0(OUTPUT, 'patterns_medallerg', timecode, '.csv'))
A <- rbind(A, fread(paste0(OUTPUT, 'example1_medallerg.csv')))
A <- rbind(A, fread(paste0(OUTPUT, 'example2_medallerg.csv')))
A <- A[sample(.N, size = .N, replace = FALSE)]
fwrite(A, paste0(OUTPUT, 'patterns_medallerg.csv'))

B <- fread(paste0(OUTPUT, 'patterns_problems', timecode, '.csv'))
B <- B[sample(.N, size = .N, replace = FALSE)]
fwrite(B, paste0(OUTPUT, 'patterns_problems.csv'))

#### OPTIONAL: SANITIZE METACAT DATA (NOT USED)

# Remove all punctuation except ? and .
# Standardise the spacing (one space before and after)
# No longer used - sanitization happens at the point of export

sanitize_metamodeldata <- function(filename){
	A <- fread(filename)
	fwrite(A, paste0(filename, 'BACKUP.csv'))
	A[, text := gsub('[^[:alnum:]\\.\\?\n ]', '', text)]
	A[, text := gsub(' +', ' ', text)]
	fwrite(A, filename)
}

#sanitize_metamodeldata(paste0(OUTPUT, 'patterns_medallerg27.Sep.17.31.46.csv'))
#sanitize_metamodeldata(paste0(OUTPUT, 'patterns_problems27.Sep.17.31.46.csv'))

#### COPY OUTPUTS TO GITHUB

file.copy(from = paste0(OUTPUT, 'patterns_medallerg.csv'),
	to = '~/GitHub/miade-datasets/cdb_and_model_files_sep_2023/')
file.copy(from = paste0(OUTPUT, 'patterns_problems.csv'),
	to = '~/GitHub/miade-datasets/cdb_and_model_files_sep_2023/')

cat('\nDONE\n')
sink()
