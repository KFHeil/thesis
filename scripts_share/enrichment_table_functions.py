#file with functions to generate
#node table including gene links to a disease/diseases/traits of interest
#two input files to compute enrichment with CDM suite
#flatfile with: traitID - trait name - geneID
#clustering file with: geneID - community
	#both tab separated

import os
from string import strip
import sys

###
#the following functions extracts all PD-GWAS from the input data
#the geneIDs are keys of a dictionary - values are lists with "u", "o", "d" respectively if the SNP is upstream, overlapping or downstream of a gene.

def read_PD_interest(infile):
	#we generate a dictionary to safe trait hits of interest
	PD_interest = {}
	for line in open(infile, "ro"):
		row = map(strip, line.split("\t"))
		#we check the start of the first row - we jump it if it is not numeric (this jumps the header line) 
		if not row[0].isdigit():
			continue
		#we are interested in saving the geneID as the key and keep the information source as a value
		PD_interest.setdefault(row[0],row[2])
	# print "dictionary of interest generated"
	return PD_interest


def generate_node_list(PPI_infile):
	#we generate a list of all the geneIDs in the PPI list of interest
	nodes = set()
	for line in open(PPI_infile, "ro"):
		row = map(strip, line.split("\t"))
		# nodes.add(row[0])
		nodes.update([row[0], row[1]])
		# print nodes, "nodes in function"
	return nodes

def nodes_communities(communities_infile):
	#gene_nodes, 
	#we read the community infile and generate
	#we generate a dictionary to safe the geneID - community pairs - the community is stored in a subdictionary called "community"
	dictionary = {}
	for line in open(communities_infile, "ro"):
		#we split each line at the empty space - elements 0 (geneID) and 2 (community) are of interest
		row = map(strip, line.split(" "))
		#we check that we are not reading the header - if so we ignore it
		if not row[0].isdigit():
			continue
		dictionary.setdefault(row[0], {})
		dictionary[row[0]].setdefault("community", row[2])
	return dictionary

# TWO traits

# def add_trait(trait_list, trait_list_1, dictionary, trait, trait_1):
# 	#we copy the input dictionary and add trait information to it
# 	# counter = 0
# 	dictionary_extended = dictionary
# 	for gene in dictionary.keys():
# 		if gene in trait_list:
# 			dictionary_extended[gene].setdefault(trait,"1")
# 			# counter += 1
# 		else: 
# 			dictionary_extended[gene].setdefault(trait, "0")
# 		#we add the second trait
# 		if gene in trait_list_1:
# 			dictionary_extended[gene].setdefault(trait_1,"1")
# 		else:
# 			dictionary_extended[gene].setdefault(trait_1, "0")
# 	# print "counter", counter
# 	return dictionary_extended

# #function to print the trait_flatfile
# def print_trait_flatfile(dictionary, trait, trait1, output):
# 	#we iter over the dictionary
# 	#and add values to an outrow
# 	#we generate and open the outfile
# 	# print "trait", trait
# 	out = open(os.path.join(output + "/" + trait + "_flatfile.csv"), "wo")
# 	print >> out, ("%s") % ("#flatfile")
# 	for key in dictionary:
# 		#we set up the outrow to add values
# 		outrow = []
# 		#we fill the output row with respective information
# 		if not dictionary[key][trait] == "0":
# 			outrow.append(dictionary[key][trait])
# 			outrow.append(trait)
# 			outrow.append(key)
# 			# print outrow, "outrow"
# 			print >> out, ("%s") % "\t".join(outrow[0:])
# 			outrow = []
# 		if not dictionary[key][trait1] == "0":
# 			outrow.append(str(int(dictionary[key][trait1]) + 1))
# 			outrow.append(trait1)
# 			outrow.append(key)
# 			# print outrow, "outrow"
# 			print >> out, ("%s") % "\t".join(outrow[0:])
# 	out.close()


############
# code as above but for only one trait

def add_trait(trait_list, dictionary, trait):
	#we copy the input dictionary and add trait information to it
	# counter = 0
	dictionary_extended = dictionary
	for gene in dictionary.keys():
		if gene in trait_list:
			dictionary_extended[gene].setdefault(trait,"1")
			# counter += 1
		else: 
			dictionary_extended[gene].setdefault(trait, "0")
	# print "counter", counter
	return dictionary_extended

# function to print the trait_flatfile
def print_trait_flatfile(dictionary, trait, output):
	#we iter over the dictionary
	#and add values to an outrow
	#we generate and open the outfile
	out = open("".join(output + "/" + trait + "_flatfile.csv"), "wo")
	print >> out, ("%s") % ("#flatfile")
	for key in dictionary:
		#we set up the outrow to add values
		outrow = []
		#we fill the output row with respective information
		if not dictionary[key][trait] == "0":
			outrow.append(dictionary[key][trait])
			outrow.append(trait)
			outrow.append(key)
			print >> out, ("%s") % "\t".join(outrow[0:])
	out.close()



def print_community_outfile(dictionary, output):
	#we generate the outfile
	out = open(os.path.join(output, "communities.csv"), "wo")
	print >> out, ("%s") % ("#clustering")
	#we iter over the dictionary and print to the outfile
	for key in dictionary:
		print >> out, ("%s\t%s") % (key, dictionary[key]["community"])
	out.close()


def print_igraph_node_table(dictionary, output, trait):
	#we iter over the dictionary, generate the output row and print it to an output file
	out = open(os.path.join(output, "node_table_igraph.csv"), "wo")
	print >> out, ("%s\t%s\t%s") % ("GeneID", "community", trait)
	for key in dictionary:
		print >> out, ("%s\t%s\t%s") % (key, dictionary[key]["community"], dictionary[key][trait])

	out.close()

# TWO traits
# def print_igraph_node_table(dictionary, output, trait, trait1):
# 	#we iter over the dictionary, generate the output row and print it to an output file
# 	out = open(os.path.join(output + "/node_table_igraph.csv"), "wo")
# 	print >> out, ("%s\t%s\t%s\t%s") % ("GeneID", "community", trait, trait1)
# 	for key in dictionary:
# 		print >> out, ("%s\t%s\t%s\t%s") % (key, dictionary[key]["community"], dictionary[key][trait], dictionary[key][trait1])

# 	out.close()
