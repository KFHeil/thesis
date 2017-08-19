#script listing functions related to data sorting
#the different functions read and filter synaptic data
#based on the given command line input
#functions are used by generate_input_data.py

import os
from string import strip
import sys

"""
function to read the input data
gene_IDs are the key
values are a list of studies where the gene appears
"""
def generate_protein_input_list(infiles):
	# print "infiles", infiles
	#we generate a dictionary
	genes_studies = {}
	#we generate a list for all the studies in the dataset:
	studies = []
	"""
	we read each of the input files - column by column
	if a 1 is present (compared to a 0) and the gene is
	not yet present in the dictionary we add it
	we add the study name (same column, first row)
	to the value-list 
	"""
	for infile in infiles:
		print infile, "infile"
		for line in open(infile, "ro"):
			#we check to identify the first row
			row = map(strip, line.split("\t"))
			if row[0].startswith("MOUSE.ENTREZ.ID"):
				#we safe the studies of interest in a list
				studies = row[4:]
				#we exit the iteration
				continue

			# print "studies", studies
			"""
			we iter over the number of studies to add genes to the dictionary (if present - indicated through "1") - study names are added to values
			"""
			for x in range(0,len(studies)):
				#we check if the gene is present in the study of interest
				if row[4+x] == "1":
					#if so we add the gene as a key and a list of studies as the value
					genes_studies.setdefault(row[2], [])
					genes_studies[row[2]].append(studies[x])
	return genes_studies

def filter_studies(studies, dictionary):
	#we generate a second dictionary with all the positive hits (of interest)
	genes_selected_studies = {}
	for study in studies:
		# print study, "study"
		# print dictionary.values()
		if not study in [x for v in dictionary.values() for x in v]:
			sys.exit("study: not in dictionary\ncheck if the correct input file was read or study correctly indicated")
		for key, value in dictionary.items():
			if study in value:
				genes_selected_studies[key] = value
	return genes_selected_studies

def filter_coverage(coverage, dictionary):
	#we filter the data dictionary based on the number of studies the genes of interest apear in
	#by default coverage is set to 1
	#the output is saved in a list
	genes_interest = []
	for key, value in dictionary.items():
		if len(value) >= coverage:
			genes_interest.append(key)
		#try to add warning if coverage > number of studies 
			#possibly add to function above..
	return genes_interest

"""
function to read the PPI input data
dictionaries have PPI-pairs as keys and subdictionaries with information interesting for filtering the data
After the filtering a PPI-list will be generated
PPIs are not sored (A-B and B-A are two entries - possibly different information annotated to each)
"""

def read_PPI_input(PPI_infile):
	#we read the file and generate a PPI_dictionary
	PPI_dict = {}
	for line in open(PPI_infile, "ro"):
		row = map(strip, line.split("\t"))
		if row[0].startswith("#ID(s) interactor A"):
			continue
		PPI = "_".join([row[0]] + [row[1]])
		PPI_dict.setdefault(PPI, {})
		#we set up the subdictionaries
		PPI_dict[PPI].setdefault("int_det_met", row[6])
		PPI_dict[PPI].setdefault("taxid", "_".join([row[9]] + [row[10]]))
		PPI_dict[PPI].setdefault("int_type", row[11])
		PPI_dict[PPI].setdefault("pubmedID", row[8])
		PPI_dict[PPI].setdefault("source_dbs", row[12])

		# this filtering step is now done later on
		# PPI_dict[PPI].setdefault("int_det_met", separator_to_list(row[6]))
		# PPI_dict[PPI].setdefault("taxid", "_".join([row[9]] + [row[10]]))
		# PPI_dict[PPI].setdefault("int_type", separator_to_list(row[11]))
		# PPI_dict[PPI].setdefault("source_dbs", separator_to_list(row[12]))


	return PPI_dict

#the following function splits "|" separated single values to lists
# def separator_to_list(entry):
# 	elements = entry.split("|")
# 	return elements

#we generate functions to filter the PPI data for:
	#int_det_met - interaction detection method
	#taxid - this is generally already done, but can be repeated here
	#int_type - interaction type
	#source database
	#database coverage

def PPI_filter_int_det_method(int_det_list, dictionary):
	#we generate a new dictionary for the filtered output
	dict_filt = {}
	#we iter over the different interaction detection methods to extract all the entries that contain at least one of them
	# this function might need the previous step of single entries per row - but filtering is currently carried out in a previous step
	for int_det_method in int_det_list:
		# print "int_det_method", int_det_method
		#check to add a warning if the ID exists or is misspelled...
		for key, value in dictionary.items():
			# print key, value, "key, value"
			for key_1, value_1 in value.items():
				# print "val", val
				#we check for the correct value and check the presence of our interaction detection method ID
				if key_1 == "int_det_met":
					# print "key_1", key_1
					# print "value_1", value_1
					if int_det_method in value_1:
						# print "key_1, value_1", key_1, value_1
						dict_filt[key] = value
	#we return the filtered dictionary
	return dict_filt

def PPI_filter_taxID(taxID, dictionary):
	# print taxID
	#we generate a new dictionary for the filtered output
	dict_filt = {}
	#we check if the taxID is as provided/the one we want to filter for
	#currently this only works if both proteins/genes have the same taxid - can be adjusted to two later on!
	for key, value in dictionary.items():
		for key_1, value_1 in value.items():
			#we check for the correct value and check the presence of our interaction detection method ID
			if key_1 == "taxid":
				#we check the taxid
				if value_1 == "_".join([taxID] + [taxID]):
					dict_filt[key] = value
				else: 
					print "no valid taxID supplied"
				#modify this message/adjust once using data from other species!
	#we return the filtered dictionary
	return dict_filt

def PPI_filter_int_type(int_type_list, dictionary):
	#we generate a new dictionary for the filtered output
	dict_filt = {}
	#we iter over the different interaction types to extract all the entries that contain at least one of the ones of interest
	for int_type in int_type_list:
		# print "int_type", int_type
		#check to add a warning if the ID exists or is misspelled...
		for key, value in dictionary.items():
			# print key, value, "key, value"
			for key_1, value_1 in value.items():
				# print "key_1, value_1", key_1, value_1
				#we check for the correct value and check the presence of our interaction detection type ID
				if key_1 == "int_type":
					# print "key_1", key_1
					# print "value_1", value_1
					if int_type in value_1:
						# print "key_1, value_1", key_1, value_1
						dict_filt[key] = value
	#we return the filtered dictionary
	return dict_filt

def PPI_filter_sourcedb_coverage(coverage, dictionary):
	#we generate a new dictionary for the filtered output
	dict_filt = {}
	#we iter over the different interaction types to extract all the entries that appear in equal or more than "coverage" databases
	for key, value in dictionary.items():
		for key_1, value_1 in value.items():
			# print "key_1, value_1", key_1, value_1
			#we check for the correct value/here source_dbs and check in how many databases the interaction is listed
			if key_1 == "source_dbs":
				# print "key_1", key_1
				# print "value_1", value_1
				if len(value_1) >= coverage:
					# print "key_1, value_1, len(value_1)", key_1, value_1, len(value_1)
					dict_filt[key] = value
	#we return the filtered dictionary
	return dict_filt

def PPI_filter_sourcedb(db_int_list, dictionary):
	#we generate a new dictionary for the filtered output
	dict_filt = {}
	#we iter over the different database IDs to extract all the entries that contain at least one of the ones of interest
	for DB in db_int_list:
		# print "DB", DB
		#check to add a warning
		for key, value in dictionary.items():
			# print key, value, "key, value"
			for key_1, value_1 in value.items():
				# print "key_1, value_1", key_1, value_1
				#we check for the correct value and check the presence of our interaction detection type ID
				if key_1 == "source_dbs":
					# print "key_1", key_1
					# print "value_1", value_1
					if DB in value_1:
						# print "key_1, value_1", key_1, value_1
						dict_filt[key] = value
	#we return the filtered dictionary
	return dict_filt
		#, final_PPIs_interest

###########
#the following function filters all internatl interaction for the geneset of interest, based on the defined filter parameters
#it generates a PPI list with "A_B" interactions
#and safes this output in the PPI_int_outfile in two columns, tab separated
def get_PPIs_interest(genes_interest, PPIs, output, logfile):
	#we generate a list for PPIs of interest
	PPI_list = []
	#we generate a list to collect all genes present in an internal interaction
	gene_list = set()
	#we iter over the PPI list and keep interactions if both geneIDs appera in out gene list of interest
	for PPI in PPIs:
		IDs = PPI.split("_")
		if IDs[0] in genes_interest and IDs[1] in genes_interest:
			PPI_list.append(PPI)
			gene_list.add(IDs[0])
			gene_list.add(IDs[1])
	#we generate the path to the output file 
	out_file = open(os.path.join(output, "PPIs_interest.csv"), "wo")
	print >> logfile, "The final PPI list is stored in the out_file: ", out_file, "\n"
	#we print the PPIs to an output file
	for PPI in set(PPI_list):
		IDs = PPI.split("_")
		print >> out_file, ("%s\t%s") % (IDs[0], IDs[1])
	out_file.close()
	return PPI_list, gene_list
