#script to run the first round of the CMD suite - clustering algorithm
#this can be downloaded from
#https://sourceforge.net/projects/cdmsuite/files/CDMSuite_cpp_v1r1/
#the home folder as well as the folder CDMSuite_cpp are hard coded

import os, sys
import shutil
import getopt

#we initiate the option to read command line arguments

def main(argv):
	print "running script 2"
	#we obtain the current working directory and safe this as the home-folder
	home = ""

	#the following parameters describe the name of the output folder to store the PPIs of interest (internal interactions only)
	PPI_int_output = ""

	try:
		opts, args = getopt.getopt(sys.argv[1:],"ho:e:",["help=", "PPIoutput=", "home="])
	#except appears in case the flag is given, but no argument supplied
	except getopt.GetoptError:
		print 'CDM_clustering.py -h <--help> \n -o <--PPIoutput> - folder name for the output file with PPIs of interest \n -e <--home> - home folder (where all other folders are'
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h', "--help"):
			print'CDM_clustering.py -h <--help> \n -o <--PPIoutput> - folder name for the output file with PPIs of interest \n -e <--home> - home folder (where all other folders are'
			sys.exit()
		#we obtain the information of the in-/output folder name
		elif opt in ("-o", "--PPIoutput"):
			PPI_int_output = arg
		elif opt in ("-e", "--home"):
			home = arg

	#we set the path to the CDMSuite_cpp folder
	CDMSuite_cpp = os.path.join(home, "CDMSuite_cpp")

	#we change the home directory to the folder with output folders
	home = os.path.join(home, "002_network_out")
	print "home", home

	#we re-open the log_file and print information of interest to this file!
	log_file = open(os.path.join(home, PPI_int_output, "log_file.txt"), "a")
	print log_file
	print >> log_file, "\nOutput from the second script \"2_CDM_clustering\""
	print >> log_file, "\nCDMSuite_cpp", CDMSuite_cpp
	#we will add input information regarding the algorithm of interest 
	#-a - algorithms to execute: 1 = Geodesic egde Betweenness, 2 = Random Walk edge Betweenness, 3 = Spectral Modularity
	#it is set to the spectral algorithm as the default option
	algorithm = "3"
	#we add the option to specify if edge lists are weighted or not
	#this is indicated by the number of columns in the input file
	#-cols - 2 if only edges; 3 if weighted edges
	#the default option is set to 2 - unweighted
	columns = "2"

	#we define the input file
	input_file = os.path.join(home, PPI_int_output, "PPIs_interest.csv")
		#test.csv is the command line input as the outfile name during the edge-list generation

	#we generate the command to run Colins spectral algorithm from the command line
	#we empty the output folder for the algorithm output by deleting all files in it
	for the_file in os.listdir(os.path.join(CDMSuite_cpp, "OUT")):
		file_path = os.path.join(CDMSuite_cpp, "OUT", the_file)
		try:
			if os.path.isfile(file_path):
				os.unlink(file_path)
		except Exception as e:
			print(e)

	command = "".join(CDMSuite_cpp + "/run -a " + algorithm + " -cols " + columns + " -file " + input_file)
	print >> log_file, "command", command

	#the command specifies ./run needs to be run from the CDMSuite_cpp folder and requires the following command line input
		#-file - edge file
		#-a - algorithms to execute: 1 = Geodesic egde Betweenness, 2 = Random Walk edge Betweenness, 3 = Spectral Modularity
		#-cols - 2 if only edges; 3 if weighted edges

	#we change to the CDMSuite_cpp working directory
	os.chdir(CDMSuite_cpp)
	#we run the command
	os.system(command)
	# print "data:", data
	# print "command", command
	# print "command run"

	#######################################
	#######################################
	#after running the algorithm we need to move the files of interest to a "safe" location

	#we generate the path to the folder with the generated output information
	communities_out = os.path.join(CDMSuite_cpp, "OUT/")
	store_data_out = os.path.join(home, PPI_int_output)

	#we move the specific output files to the folder we store all the output information in
	#we generate the command to move the files
	command_1 = ''.join("mv " + communities_out + "* " + store_data_out)
	os.system(command_1)
	print >> log_file, "command_1", command_1
	print >> log_file, "command 1 completed - files moved to test output folder: ", store_data_out, "\n"

	print "script 2 finished\n"
if __name__ == "__main__":
	main(sys.argv[1:])
