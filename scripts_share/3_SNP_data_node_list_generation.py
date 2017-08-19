#!/usr/bin/python

##file calling functions to generate tables with geneIDs of interest together with information about their association to a trait (e.g. PD-GWAS) and their grouping to communities

# home folder as well as PD list are hard coded

from enrichment_table_functions import *
import getopt
import os

#we initiate the option to read command line options
def main(argv):

	#we obtain the current working directory and safe this as the home-folder
	home = ""
	
	#the following parameters describe the name of the output folder to store the PPIs of interest (internal interactions only)
	PPI_int_output = ""

	try:
		opts, args = getopt.getopt(sys.argv[1:],"ho:e:",["help=", "PPIoutput=", "home="])

	#except appears in case the flag is given, but no argument supplied
	except getopt.GetoptError:
		print 'CDM_clustering.py -h <--help> \n -o <--PPIoutput> - folder name for the output file with PPIs of interest\n -e <--home> - home folder (where all other folders are'
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h', "--help"):
			print'CDM_clustering.py -h <--help> \n -o <--PPIoutput> - folder name for the output file with PPIs of interest\n -e <--home> - home folder (where all other folders are'
			sys.exit()
		#we obtain the information of the in-/output folder name
		elif opt in ("-o", "--PPIoutput"):
			PPI_int_output = arg
		elif opt in ("-e", "--home"):
			home = arg

	##file listing the genes associated to PD
	#we define the datafile of interest
	PD_trait_infile = os.path.join(home, "111_data/PD_final_list.csv")
	#we define the trait name
	trait = "PD"
	#we call the function to generate a dictionary with values of interest
	PD_hits_dict = read_PD_interest(PD_trait_infile)
	# print PD_hits_dict, "PD_hits_dict"

	#we change the home directory to the folder with output folders
	home = os.path.join(home, "002_network_out")
	print "home", home

	#we re-open the log_file and print information of interest to this file!
	log_file = open(os.path.join(home, PPI_int_output, "log_file.txt"), "a")
	print log_file
	print >> log_file, "\nWe are running the third script \"3_SNP_data_node_list_generation.py\""

	###########
	# here we could add some code to consider the different PD sources independently
	#we filter the GWAS_hits_dict for only overlapping SNPs
	#we generate a list of those EntrezIDs
	# genes_SNPs_overlap = set()
	# for key, value in GWAS_hits_dict.items():
	# 	if "o" in value:
	# 		print value, "value"
	# 		genes_SNPs_overlap.add(key)

	# print "genes_SNPs_overlap", genes_SNPs_overlap
	# print len(genes_SNPs_overlap), "len(genes_SNPs_overlap)"

	#we define the input edge list
	edges_in = os.path.join(home, PPI_int_output, "PPIs_interest.csv")

	nodes = generate_node_list(edges_in)

	# print nodes[:15], "nodes[:15]"
	print >> log_file, "The number of genes in the PPI file is:", len(nodes), "\n"

	#we call the function to generate a first dictionary with all nodes and the respective community they belong to 
	#to add community information to the nodes in the network
	communities_info_in = os.path.join(home, PPI_int_output, "communities_cytoscape.txt")
	dict_0 = nodes_communities(communities_info_in)

	#we call the function that adds trait information to the dictionary
	#the function takes three arguments - the list of geneIDs of interest, the dictionary and a name of the trait

	dict_1 = add_trait(PD_hits_dict.keys(), dict_0, trait)
	# for key, value in dict_1.items():
	# 	for key1, value1 in value.items():
	# 		if key1 == "GWAS_overlap" and value1 == "1":
	# 			print "value1", value1, "key1", key1
	# 			# for value1 in dict_1[value]:
					# if value1 == "1":
						# print key, value, value1, "key", "value", "value1"

	# print len(dict_1.keys()), "len(dict_1.keys())"
	###
	#the following functions print the output 

	# outfile = os.path.join()
	outfolder = os.path.join(home, PPI_int_output)

	print_trait_flatfile(dict_1, trait, outfolder)
	#we print the trait file for the second trait
	# print_trait_flatfile(dict_1, trait1, outfolder)

	#the following function prints the outfile including geneID and community number
	print_community_outfile(dict_1, outfolder)

	#the following function 
	print_igraph_node_table(dict_1, outfolder, trait)

	print >> log_file, "A geneID dictionary was generated including information about gene presence in the trait dataset:", trait, " and the community they belong to\n"

if __name__ == "__main__":
	main(sys.argv[1:])
