#!/usr/bin/python

#script to generate input data for synaptic network analysis
#command line input is required
#functions are saved in data_reading_functions.py

from data_reading_functions import *
#we import required modules 
import sys, getopt
import os

def main(argv):
	print "running script 1"
	#we define the homedirectory - currently hard-coded
	home = ""

	#the path to the pre, post and synaptosome data are hardcoded
	pre = "111_data/Presynaptic(Published)_4.csv"
	post = "111_data/PSD0_full_Published_final_4.csv"
	synapt = "111_data/synaptosome_Published_4.csv"
	#the following lists saves the datainput locations of interest (see links above)
	in_genes = []

	try:
		opts, args = getopt.getopt(sys.argv[1:],"hi:d:g:p:m:t:n:c:b:o:e:",["incompartment=", "datasets=", "genecoverage=", "ppiinfile=", "intmethod=", "taxid=", "inttype=", "intcoverage=", "dbtypes=", "PPIoutput=", "home="])
		print opts, args, "opts, args"
	#except appears in case the flag is given, but no argument supplied
	except getopt.GetoptError:
		print 'test_command_line_input.py -i <--incompartment> takes the following arguments: pre, post, synapt - input in form of comma separated list \n-d <--datasets> - specify datasets to be included - comma separated, but no spaces \n-g <--genecoverage> - number of studies in which the genes need to appear to be included in the gene list\n -p <--ppiinfile> - specification of an altered PPI infile - otherwise this is hard coded\n -m <--intmethod> - comma separated list (no spaces) of interaction method filtering - if not given no filters are applied\n -t <--taxid> - if not specified differently the taxid is set to human for both interacting genes\n -n <--inttype> - list with interaction types of interest - comma separated - no spaces\n -c <--intcoverage> - coverage/number of databases listing a PPI - default set to 1\n -b <--dbtypes> - list of databases we filter PPIs for - comma separated, no spaeces\n -o <--PPIoutput> - folder name for the output file with PPIs of interest\n -e <--home> - home folder - where folder to output folders and data folder as well as CDM suits are stored (no dash at the end of path)'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'test_command_line_input.py\n-i <--incompartment> takes the following arguments: pre, post, synapt - input in form of comma separated list \n-d <--datasets> - specify datasets to be included - comma separated, but no spaces \n-g <--genecoverage> - number of studies in which the genes need to appear to be included in the gene list\n-p <--ppiinfile> - specification of an altered PPI infile - otherwise this is hard coded\n-m <--intmethod> - comma separated list (no spaces) of interaction method filtering - if not given no filters are applied\n-t <--taxid> - if not specified differently the taxid is set to human for both interacting genes\n-n <--inttype> - list with interaction types of interest - comma separated - no spaces\n-c <--intcoverage> - coverage/number of databases listing a PPI - default set to 1\n-b <--dbtypes> - list of databases we filter PPIs for - comma separated, no spaeces\n-o <--PPIoutput> - folder name for the output file with PPIs of interest\n -e <--home> - home folder - where folder to output folders and data folder as well as CDM suits are stored (no dash at the end of path)'
			sys.exit()
		# elif opt in ("-i", "--incompartment"):
		# elif opt == "-i" or opt == "--incompartment":
		elif opt in ("-i", "--incompartment"):
			if "pre" in arg:
				in_genes.append(pre)
			if "post" in arg:
				in_genes.append(post)
			if "synapt" in arg: 
				in_genes.append(synapt)
			# 	print "Error: no correct input file type supplied"
			# 	sys.exit()
		elif opt in ("-d", "--datasets"):
			#we check if the parameter is empty - indicatd with "-" - if so we set the datasets to empty, else we add the arguments to the datasets variable
			if arg == "-":
				datasets = []
			else:
				datasets = arg.split(",")
		elif opt in ("-g" , "--genecoverage"):
			coverage_gene = int(arg)
		#PPI filtering
		elif opt in ("-p", "--ppiinfile"):
			PPI_infile = arg
		elif opt in ("-m", "--intmethod"):
			#we check as for datasets
			if arg == "-":
				int_method_interest = []
			else:
				int_method_interest = arg.split(",")
		elif opt in ("-t", "--taxid"):
			taxID = arg
		elif opt in ("-n", "--inttype"):
			if arg == "-":
				interaction_type_interest = []
			else:
				interaction_type_interest = arg.split(",")
		elif opt in ("-c", "--intcoverage"):
			coverage_interaction = int(arg)
		elif opt in ("-b", "--dbtypes"):
			#we check as for databases
			if arg == "-":
				DB_types = []
			else: 
				DB_types = arg.split(",")
		elif opt in ("-o", "--PPIoutput"):
			PPI_int_output = arg
		elif opt in ("-e", "--home"):
			home = arg

	#we generate the path to the outfolder
	outfolder = os.path.join(home, "002_network_out", PPI_int_output)
	print "outfolder - file 1", outfolder
	
	# we start printing information to the log file
	log_file = open(os.path.join(outfolder + "/log_file.txt"), "w")
	print "log_file", log_file

	print >> log_file, ("%s") % ("we are running script 1 - 1_generate_filtered_PPI_list_interest.py")
	print >> log_file, ("%s %s") % (log_file, "log_file")

	print >> log_file, "\ninput parameters"
	print >> log_file, "###################"
	print >> log_file, ("%s %s") % (in_genes, "in_genes")
	print >> log_file, datasets, "datasets (default no filter)"
	print >> log_file, coverage_gene, "coverage_gene (default: 1)"
	print >> log_file, "\nthe following arguments modify/filter the PPI data"
	print >> log_file, PPI_infile, "PPI_infile (default set to most recent clean PPI-set"
	print >> log_file, int_method_interest, "interaction methods specified for filtering (default no filter)"
	print >> log_file, taxID, "taxID (default: human: 9606)"
	print >> log_file, interaction_type_interest, "interaction_type_interest (default filter for all direct interactions) - if no filter - set command line input to []"
	print >> log_file, coverage_interaction, "coverage_interaction (default: 1)"
	print >> log_file, DB_types, "DB_types (default no filter)"
	print >> log_file, "\n the following parameter describes the name of the output folder for PPIs of interest (internal) and the log file, as well as output from following scripts"
	print >> log_file, PPI_int_output, "PPI_int_output"
	print >> log_file, home, "path to home directory"
	print >> log_file, "###################\n"
	########################################################
	#we generate the dictionary with genes of interest

	#the input to this function is the list of input files
	if not len(in_genes) == 0:
		print "running generate_protein_input_list"
		#we modify the in_genes path to add the home directory
		in_genes_full = []
		for entry in in_genes:
			new_path = os.path.join(home, entry)
			in_genes_full.append(new_path)
		print "in_genes_full", in_genes_full
		out_dict = generate_protein_input_list(in_genes_full)
		# print "upper level file - out_dict", out_dict
		print >> log_file, "initial dictionay generated"
		print >> log_file, "it includes genes from the following datasets"
		print >> log_file, in_genes, "\n"
		print >> log_file, "The number of genes in the initial dictionary is", len(out_dict.keys()), "\n"
	else: 
		print "WARNING!!! \nno gene list of interest provided"
		print "see filename -h for help"
		print >> log_file, "WARNING!!! \nno gene list of interest provided"
		print >> log_file, "see filename -h for help"

	#we pass the dataset list to the filter function
	#if none is specified we assume, that no filter is applied and all studies are kept
	#otherwise genes from the studies of interest are extracted
	# print "datasets", datasets, "\n"
	if not len(datasets) == 0:
		out_dict_selected_studies = filter_studies(datasets, out_dict)
		print >> log_file, "dictionary filtered for the studies of interest"
		print >> log_file, "the included studies are:", datasets
		print >> log_file, "this minimizes the list of genes of interest to: ", len(out_dict_selected_studies.keys()), "genes\n"
	else: 
		out_dict_selected_studies = out_dict
		print >> log_file, "not filtered for studies - all studies included"
		print >> log_file, "The number of genes remains at:", len(out_dict_selected_studies.keys()), "genes\n"

	#we generate the final gene list
	#one further filtering step, if wanted, can be applied and the coverage of a studies listing the gene of interest can be checked
	final_gene_list = filter_coverage(coverage_gene, out_dict_selected_studies)

	print len(final_gene_list), "len(final_gene_list)"

	print >> log_file, "the gene list was filtered for a certain coverage of studies per gene entry"
	print >> log_file, "that coverage was set to ", coverage_gene
	print >> log_file, "the final_gene_list - after applying all the filtering options consists of", len(final_gene_list), "genes\n"

	"""
	the next part of code reads and filters the PPI set of interest. 
	This can then be used to generate the edge list of interest
	"""

	PPI_first_dict = read_PPI_input(PPI_infile)
	print >> log_file, "###################"
	print >> log_file, "the first PPI interaction dictionary is build"
	print >> log_file, "The number of PPIs is:", len(PPI_first_dict.values()), "\n"

	print len(PPI_first_dict.keys()), "len(PPI_first_dict.keys())"
	#we filter the PPI data for interaction methods of interest
	#if no interaction methods are specified no filtering occurs
	if not len(int_method_interest) == 0:
		PPI_first_filtered = PPI_filter_int_det_method(int_method_interest, PPI_first_dict)
		print >> log_file, "The PPI list was filtered for interaction methods of interest"
		print >> log_file, "The number of remaining PPIs is:", len(PPI_first_filtered), "\n"
	else: 
		PPI_first_filtered = PPI_first_dict
		print >> log_file, "The PPI list was not filtered for interaction methods (no filters supplied)"
		print >> log_file, "The number of PPIs remains the same:", len(PPI_first_filtered.values()), "\n"

	#we apply a filter for taxIDs
	PPI_second_filtered = PPI_filter_taxID(taxID, PPI_first_filtered)
	print >> log_file, "after filtering for the taxid (default - human:", taxID, ") the number of PPIs is: ", len(PPI_second_filtered), "\n"

	#we apply a filter for the interaction type of interest
	#if filters are given
	if not len(interaction_type_interest) == 0:
		PPI_third_filtered = PPI_filter_int_type(interaction_type_interest, PPI_second_filtered)
		print >> log_file, "The PPIs were filtered for the following interaction types:", interaction_type_interest
		#add option to filter e.g. for direct interactions! 
		print >> log_file, "The number of PPIs after filtering for interaction types of interest is:", len(PPI_third_filtered), "\n"
		print "The PPIs were filtered for the following interaction types:", interaction_type_interest
		print "The number of PPIs after filtering for interaction types of interest is:", len(PPI_third_filtered), "\n"
	else:
		PPI_third_filtered = PPI_second_filtered
		print >> log_file, "No interaction type filter was applied - the number of PPIs remains at:", len(PPI_third_filtered), "\n"

	print len(PPI_third_filtered), "len(PPI_third_filtered)"

	#we filter for the number of databases that backup/list an interaction
	#the default if not specified differently is 1
	PPI_fourth_filtered = PPI_filter_sourcedb_coverage(coverage_interaction, PPI_third_filtered)
	print >> log_file, "The data is filtered for the coverage:", coverage_interaction,  "- how many databases list the PPI of interest"
	print >> log_file, "The number of PPIs after applying a filter on the coverage is:", len(PPI_fourth_filtered), "\n"
	print "The number of PPIs after applying a filter on the coverage is:", len(PPI_fourth_filtered), "\n"

	#filtering for presence of PPIs in specific databases
	if not len(DB_types) == 0:
		PPI_fifth_filtered = PPI_filter_sourcedb(DB_types, PPI_fourth_filtered)
		print >> log_file, "A filter for specific databases was applied - the included databases are:", DB_types
		print >> log_file, "The number of PPIs after filtering for presence in specific databases is: ", len(PPI_fifth_filtered), "\n"
		print "The number of PPIs after filtering for presence in specific databases is: ", len(PPI_fifth_filtered), "\n"
	else:
		PPI_fifth_filtered = PPI_fourth_filtered
		print >> log_file, "No filtering for specific databases was applied - the number of PPIs remains at: ", len(PPI_fifth_filtered), "\n"
		print "No filtering for specific databases was applied - the number of PPIs remains at: ", len(PPI_fifth_filtered), "\n"

	#we cross the two lists to generate a PPI list of interest
	#this is a list with "A_B" interactions that is also printed to an output file (edge_list.csv)
	PPIs_interest, genes_in_PPI_interest = get_PPIs_interest(final_gene_list, PPI_fifth_filtered.keys(), outfolder, log_file)
	print >> log_file, "The PPI list of interest based on the gene list of interest was generated\nit contains:", len(PPIs_interest), "interactions"
	print >> log_file, len(genes_in_PPI_interest), "genes are present within these interactions."
	# print "PPIs_interest[:10]", PPIs_interest[:10]

	print "script 1 finished - move on"
if __name__ == "__main__":
	main(sys.argv[1:])
