#!/usr/bin/python

#importing the subprocess function and using it to call the imported script!

#script is adjusted to 

#we also add getopt to be able to read command line input
import subprocess
import sys, os
import getopt

def main(argv):
	print "running general script"
	#we define the homedirectory
	#this is currently hard coded - some other parts as well
	home = ""
	#needs to be hardcoded when on cluster
	# home = "/exports/eddie/scratch/s1368042"
	#the following list saves the data input of interest
	#input lists (pre/post/synaptosome) of interest
	in_genes = []

	#the following list saves the datasets of interest
	datasets = []
	datasets.append("-")

	#we set the default gene coverage to 1
	coverage_gene = 1
	
	#the following parameters refer to the generation of the PPI list of interest

	#hard coded PPI infile - can be specified, otherwise this default option will be applied
	PPI_infile = ("/home/kfheil/Dropbox/THESIS/PPI_chapter/2017_04_12_9606_2017_03_PPI_direct_no_mirrored.csv")
	print "The PPI file is hard-coded!"

	#we define a filter for PPI interaction methods
	int_method_interest = []
	int_method_interest.append("-")

	#the taxID of interest is defined
	#if not differently specified the default is human
	taxID = "9606"

	#we define possible interaction type filters
	#the following list contains all direct interaction type MI-IDs (updated 2017-05-19 - includes obsolete terms)
	# this should not make any change to the numbers since the filtering was already done in a previous step, generating the PPI list
	interaction_type_interest = ["MI_1027", "MI_0556", "MI_0557", "MI_0558", "MI_0559", "MI_0570", "MI_0971", "MI_0214", "MI_0844", "MI_1127", "MI_1139", "MI_0210", "MI_1327", "MI_0414", "MI_0882", "MI_0883", "MI_0881", "MI_0945", "MI_0192", "MI_0213", "MI_0212", "MI_0211", "MI_0193", "MI_0217", "MI_0216", "MI_0197", "MI_0198", "MI_0194", "MI_0987", "MI_0986", "MI_0985", "MI_0195", "MI_1250", "MI_1237", "MI_1126", "MI_0571", "MI_0569", "MI_0568", "MI_1230", "MI_0567", "MI_0566", "MI_0701", "MI_0408", "MI_1148", "MI_1251", "MI_0915", "MI_0407", "MI_0199", "MI_1140", "MI_1310", "MI_1143", "MI_1146", "MI_0910", "MI_0902", "MI_0572", "MI_0871", "MI_0200", "MI_0201", "MI_0202", "MI_0203", "MI_0204", "MI_0206", "MI_0207", "MI_0209", "MI_0220", "MI_0914", "MI_0218"]

	#coverage threshold to keep a PPI
	#the default value is set to 1
	coverage_interaction = 1

	#we define the databases we are interested in 
	DB_types = []
	DB_types.append("-")

	#the following parameters describe the name of the output folder to store the PPIs of interest (internal interactions only)
	PPI_int_output = ""

	#for values that might not be supplied we add a "-" for future command line input

	try:
		opts, args = getopt.getopt(sys.argv[1:],"hi:d:g:p:m:t:n:o:e:",["incompartment=", "datasets=", "genecoverage=", "ppiinfile=", "intmethod=", "taxid=", "inttype=", "intcoverage=", "dbtypes=", "PPIoutput=", "home="])
		print opts, args, "opts, args"
	#except appears in case the flag is given, but no argument supplied
	except getopt.GetoptError:
		print 'test_command_line_input.py \n-i <--incompartment> takes the following arguments: pre, post, synapt, - input in form of comma separated list \n-d <--datasets> - specify datasets to be included - comma separated, but no spaces \n-g <--genecoverage> - number of studies in which the genes need to appear to be included in the gene list\n -p <--ppiinfile> - specification of an altered PPI infile - otherwise this is hard coded\n -m <--intmethod> - comma separated list (no spaces) of interaction method filtering - if not given no filters are applied\n -t <--taxid> - if not specified differently the taxid is set to human for both interacting genes\n -n <--inttype> - list with interaction types of interest - comma separated - no spaces\n -c <--intcoverage> - coverage/number of databases listing a PPI - default set to 1\n -b <--dbtypes> - list of databases we filter PPIs for - comma separated, no spaces\n -o <--PPIoutput> - folder name for the output file with PPIs of interest\n -e <--home> - home folder - where folder to output folders and data folder as well as CDM suits are stored (no dash at the end of path)'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'test_command_line_input.py \n-i <--incompartment> takes the following arguments: pre, post, synapt, - input in form of comma separated list \n-d <--datasets> - specify datasets to be included - comma separated, but no spaces \n-g <--genecoverage> - number of studies in which the genes need to appear to be included in the gene list\n -p <--ppiinfile> - specification of an altered PPI infile - otherwise this is hard coded\n -m <--intmethod> - comma separated list (no spaces) of interaction method filtering - if not given no filters are applied\n -t <--taxid> - if not specified differently the taxid is set to human for both interacting genes\n -n <--inttype> - list with interaction types of interest - comma separated - no spaces\n -c <--intcoverage> - coverage/number of databases listing a PPI - default set to 1\n -b <--dbtypes> - list of databases we filter PPIs for - comma separated, no spaces\n -o <--PPIoutput> - folder name for the output file with PPIs of interest\n -e <--home> - home folder - where folder to output folders and data folder as well as CDM suits are stored (no dash at the end of path)'
			sys.exit()
		elif opt in ("-i", "--incompartment"):
			in_genes = arg.split(",")
		elif opt in ("-d", "--datasets"):
			datasets = arg.split(",")
		elif opt in ("-g" , "--genecoverage"):
			coverage_gene = int(arg)
		elif opt in ("-p", "--ppiinfile"):
			PPI_infile = arg
		elif opt in ("-m", "--intmethod"):
			int_method_interest = arg.split(",")
		elif opt in ("-t", "--taxid"):
			taxID = arg
		elif opt in ("-n", "--inttype"):
			interaction_type_interest = arg.split(",")
		elif opt in ("-c", "--intcoverage"):
			coverage_interaction = int(arg)
		elif opt in ("-b", "--dbtypes"):
			DB_types = arg.split(",")
		elif opt in ("-o", "--PPIoutput"):
			PPI_int_output = arg
		elif opt in ("-e", "--home"):
			home = arg

	#we check if a home directory was supplied
	#if not the script stops
	if len(home) == 0:
		print "WARNING\nno home directory specified - use -e // --home to do so (no dash at end of path - home folder - where folder to output folders and data folder as well as CDM suits are stored (no dash at the end of path)"
		sys.exit()

	#we check if an output folder is supplied
	#if not the script is stopped
	if len(PPI_int_output) == 0:
		print "WARNING\nno output folder specified - use -o // --PPIoutput to do so"
		sys.exit()
	#we generate the output folder to safe a log-file to safe command line input options and output of interest
	#we check if the output folder exist - if so a warning is given and the files are not overwritten - the script is stopped
	#if it does not exist the folder is created and the output file can be written
	#we generate the path to the outfolder
	outfolder= os.path.join(home, "002_network_out", PPI_int_output)

	if not os.path.exists(outfolder):
		os.makedirs(outfolder)
		print "outfolder folder generated"
	else: 
		print "WARNING: The outfolder already exists - files are NOT overwritten"
		print "please specify an alternative outfolder!"
		sys.exit()
	log_file = open(os.path.join(outfolder, "log_file_general.txt"), "w")

	print "log_file", log_file

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
	print >> log_file, "\n the following parameter describes the name of the output folder for PPIs of interest (internal) and the log file, as well as output from following scripts"
	print >> log_file, PPI_int_output, "PPI_int_output"
	print >> log_file, home, "path to home directory"
	print >> log_file, "###################\n"
	print "log_file printed - check output folder"

	######################################################################################################################
	################
	#we define argument_1 to run script 1
	argument_1 = "".join("python " + home + "/777_scripts/1_generate_filtered_PPI_list_interest.py -i " + ",".join(in_genes) + " -d " + ",".join(datasets) + " -g " + str(coverage_gene) + " -p " + PPI_infile + " -m " + ",".join(int_method_interest) + " -t " + str(taxID) + " -n " + ",".join(interaction_type_interest) + " -c " + str(coverage_interaction) + " -b " + ",".join(DB_types) + " -o " + PPI_int_output + " -e " + home)

	print "we are now moving on to read the first script!!!!"
	print argument_1, "argument_1"
	print >> log_file, "we are now running the first script"
	print >> log_file, "argument_1", argument_1

	test = subprocess.Popen(argument_1, shell=True)
	test.communicate()
	print "argument/script 1 successfully completed (calling-script)\n"
	print >> log_file, "argument/script 1 successfully completed (calling-script)\n"

	################
	#we define the second argument to run the 2nd script
	argument_2 = "".join("python " + home + "/777_scripts/2_CDM_clustering.py -o " + PPI_int_output + " -e " + home)
	print "we are now moving on to read the second script!!!!\n"
	print argument_2, "argument_2"
	print >> log_file, "we are now running the second script"
	print >> log_file, "argument_2", argument_2
	# subprocess.call(argument_2, shell=True)
	test_2 = subprocess.Popen(argument_2, shell=True)

	test_2.communicate()
	print "argument/script 2 successfully completed (calling-script)"
	print >> log_file, "argument/script 2 successfully completed (calling-script)\n"

	################
	#we define the third argument to run the 3rd script
	argument_3 = "".join("python " + home + "/777_scripts/3_SNP_data_node_list_generation.py -o " + PPI_int_output + " -e " + home)
	print "we are now moving on to read the third script!!!!\n"
	print argument_3, "argument_3"
	print >> log_file, "we are now running the third script"
	print >> log_file, "argument_3", argument_3
	# subprocess.call(argument_2, shell=True)
	test_3 = subprocess.Popen(argument_3, shell=True)

	test_3.communicate()
	print "argument/script 3 successfully completed (calling-script)"
	print >> log_file, "argument/script 3 successfully completed (calling-script)\n"

	################
	#we define the fourth argument to run the 4th script
	argument_4 = "".join("python " + home + "/777_scripts/4_CDM_enrichment.py -o " + PPI_int_output + " -e " + home)
	print "we are now moving on to read the fourth script!!!!\n"
	print argument_4, "argument_4"
	print >> log_file, "we are now running the fourth script"
	print >> log_file, "argument_4", argument_4
	# subprocess.call(argument_2, shell=True)
	test_4 = subprocess.Popen(argument_4, shell=True)

	test_4.communicate()
	print "argument/script 4 successfully completed (calling-script)"
	print >> log_file, "argument/script 4 successfully completed (calling-script)\n"


	################
	#we define the fifth argument to run the 5th script
	argument_5 = "".join("Rscript " + home + "/777_scripts/5_generate_graphs.R -o " + PPI_int_output + " -e " + home)
	print "we are now moving on to read the fifth script!!!!\n"
	print argument_5, "argument_5"
	print >> log_file, "we are now running the fifth script"
	print >> log_file, "argument_5", argument_5
	# subprocess.call(argument_2, shell=True)
	test_5 = subprocess.Popen(argument_5, shell=True)

	test_5.communicate()
	print "argument/script 5 successfully completed (calling-script)"
	print >> log_file, "argument/script 5 successfully completed (calling-script)\n"
	################
	#we define the sixth argument to run the 6th script
	#this is a number of different runs - we iterate over the different algorithms and compute the enrichment
	#we define the list of algorithms 
	algorithms = ["fast_greedy", "infomap", "louvain", "Newman_Girvan_edge_betweenness", "spinglas", "walktrap"]
	print "algorithms", algorithms
	print >> log_file, "we set the list of algorithms", algorithms
	print >> log_file, "we iterate over the algorithms and run enrichment analysis individually"

	#we iterate over the algorithms and run the argument
	for alg in algorithms:
		print >> log_file, "the algorithm is: ", alg
		argument_6 = "".join("python " + home + "/777_scripts/4_CDM_enrichment_R_algorithms.py -o " + PPI_int_output + " -a " + alg + " -e " + home)
		print "we are running the algorithm specific enrichment!!!!\n"
		print argument_6, "argument_6 - algorithm specific"
		print >> log_file, "we are now running the sixth script"
		print >> log_file, "argument_6", argument_6
		# subprocess.call(argument_2, shell=True)
		test_6 = subprocess.Popen(argument_6, shell=True)

		test_6.communicate()
		print "argument/script 6 - algorithm specific - successfully completed (calling-script)"
		print >> log_file, "argument/script 6 - algorithm specific - successfully completed (calling-script)\n"

	#add code to move spectral algorithm files to spectral folder 
	#we generate the path to the folder with the generated output information - and generate the folder
	# /home/kfheil/Dropbox/THESIS/synaptic_PPIN_chapter
	spectral_file_communities = os.path.join(home, "002_network_out", PPI_int_output, "communities.csv")
	store_data_communities = os.path.join(home, "002_network_out", PPI_int_output, "spectral/communities_all_PPIs.csv")

	#we move the spectral algorithm output file communities to the spectral output folder - we rename it since it contains information regarding all nodes undergoing PPIs - and not only the ones in the biggest connected component
	#we generate the command to move the files
	command_move_communities = ''.join("mv " + spectral_file_communities + " " + store_data_communities)
	os.system(command_move_communities)

	#we generate a list of all the other files that need to be moved, iterate over it and move the respective files
	move_files =  ["communities_cytoscape.txt", "communityout.txt", "consensusout.txt", "enrichment_PD.csv", "permute_p_values_PD.csv", "p_sig_tests_summary_PD.csv", "p_values_PD.csv", "removededges.txt", "p_sig_tests_PD.csv"]
	for file in move_files:
		spectral_file_1 = os.path.join(home, "002_network_out", PPI_int_output, file)
		store_data_out = os.path.join(home, "002_network_out", PPI_int_output, "spectral")
		command_move = ''.join("mv " + spectral_file_1 + " " + store_data_out)
		os.system(command_move)

	log_file.close()

if __name__ == "__main__":
	main(sys.argv[1:])
