#script to run the second round of the CMD suite - enrichment analysis
#an advanced CDM package is required - available upon request

import os, sys
import shutil
import getopt

#we initiate the option to read command line options
def main(argv):

	#we obtain the current working directory and safe this as the home-folder
	home = ""

	#the following parameters describe the name of tohe output folder to store the PPIs of interest (internal interactions only)
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

	#we set the path to the CDMSuite_cpp foler
	CDM_enrichment = os.path.join(home, "CDM_enrichment")

	#we change the home directory to the folder with output folders
	home = os.path.join(home, "002_network_out")
	print "home", home

	#we re-open the log_file and print information of interest to this file!
	log_file = open(os.path.join(home, PPI_int_output, "log_file.txt"), "a")
	print log_file
	print >> log_file, "\nWe are running the fourth file: \"4_CDM_enrichment.py\""
	print >> log_file, "CDM_enrichment", CDM_enrichment

	#we set a random variable as the extension for output files
	extension = "PD"

	#we define the input file with community and trait information
	input_file_communities = os.path.join(home, PPI_int_output, "communities.csv")
	input_file_trait = os.path.join(home, PPI_int_output, "PD_flatfile.csv")

	#we generate the command to run Colins spectral enrichment tool from the command line

	#we empty the DCM output folder for the algorithm output by deleting all files in it
	for the_file in os.listdir(os.path.join(CDM_enrichment, "OUT")):
		file_path = os.path.join(CDM_enrichment, "OUT", the_file)
		try:
			if os.path.isfile(file_path):
				os.unlink(file_path)
		except Exception as e:
			print(e)

	command = "".join(CDM_enrichment + "/run " + input_file_communities + " " + input_file_trait + " " + extension)
	print >> log_file, "command", command

	#the command specifies ./run needs to be run from the CDM_enrichment folder and requires the following command line input
		#tab separated file with GeneID community
		#tab separated file with an ID for the trait of interest - the trait of interest - GeneID
		#data-name/ID - will be added to output file name

	#we change to the CDMSuite_cpp working directory
	os.chdir(CDM_enrichment)
	#we run the command
	os.system(command)

	#######################################
	#######################################
	#after running the algorithm we need to move the files of interest to a "safe" location

	#we generate the path to the folder with the generated output information
	enrichment_out = os.path.join(CDM_enrichment, "OUT/")
	store_data_out = os.path.join(home, PPI_int_output)

	#we move the specific output files to the folder we store all the output information in
	#we generate the command to move the files
	command_1 = ''.join("mv " + enrichment_out + "* " + store_data_out)
	os.system(command_1)
	print >> log_file, "command_1", command_1
	print >> log_file, "command 1 completed - files moved to test output folder: ", os.path.join(home, PPI_int_output)

if __name__ == "__main__":
	main(sys.argv[1:])
