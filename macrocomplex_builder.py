import Bio.PDB
import sys
import string
import os
import argparse
import timeit
import logging
import re
from macrocomplex_functions import *
# alias chimera="~/.local/UCSF-Chimera64-1.13.1/bin/chimera"

start = timeit.default_timer()

parser = argparse.ArgumentParser(description = 
	"This program is able to reconstruct biological macrocomplexes of protein-protein interactions as well as protein-DNA/RNA interactions given a set of binary interactions and the desired number of chains of the target complex.")

requiredNamed = parser.add_argument_group('required arguments')

requiredNamed.add_argument('-i', '--indir',			#INPUT FOLDER argument
							dest = "indir",
							action = "store",
							required=True,
							help = "Input folder (or path) containing all PDB files with the protein binary interactions. It is a compulsory argument.")

requiredNamed.add_argument('-nc', '--number_chains',		#NUMBER OF CHAINS argument
							dest = "number_chains",
							action = "store",
							required=True,
							type = int,
							help = "Number of chains desired for the target complex. This is a compulsory argument.")

parser.add_argument('-o', '--outdir',			#OUTPUT FOLDER argument
					dest = "outdir",
					action = "store",
					default = None,
					help = "If set, all the models generated in each iteration, the final macrocomplex structure in PDB format and the log file will be saved in this folder. By default, the output folder will be named as the input folder + \"_output\".")

parser.add_argument('-v', '--verbose',			#VERBOSE argument
					dest = "verbose",
					action = "store_true",
					default = False,
					help = "If set, the progression log printed in standard output file.")

parser.add_argument('-pi', '--pdb_iterations',		#PDB FILES ITERATIONS argument
					dest = "pdb_iterations",
					action = "store_true",
					default = False,
					help = "If set, each time a chain is added to the complex, a new PDB file will be saved.")

parser.add_argument('-rmsd', '--rmsd_threshold',		#RMSD THRESHOLD argument
					dest = "rmsd_threshold",
					action = "store",
					default = 0.3,
					type = float,
					help = "If set, the RMSD threshold for considering a superimposition as correct will take this value. If not, it will be 0.3 by default. The output of the program is very sensitive to this value, we advise to be careful when modifying it.")

parser.add_argument('-cl', '--clashes_threshold',		#CLASHES argument
					dest = "clashes",
					action = "store",
					default = 30,
					type = int,
					help = "If set, the threshold of the number of clashes will take this value. If not, it will be 30 by default. The output of the program is very sensitive to this value, we advise to be careful when modifying it.")

parser.add_argument('-it', '--iterations',		#ITERATIONS argument
					dest = "it",
					action = "store",
					default = 100,
					type = int,
					help = "If set, the maximum number of iterations will run for will take this value. If not, it will be 100 by default. This parameter will only be important if the complex is big and needs lots of iterations to finish.")

### Saving and checking the command-line arguments ###

arguments = parser.parse_args()

## MANDATORY arguments: INPUT and NUMBER OF CHAINS ##

if not arguments.indir:		#checking if an INPUT has been provided
	raise NameError("ERROR! The input directory has not been provided! Please, use the help flag, --help, h to see the program instructions!")
else:		#INPUT has been provided
	if (os.path.isdir(arguments.indir)):		#checking that INPUT is a real directory
		files = sorted(list(filter(lambda x: x.endswith(".pdb"), os.listdir(arguments.indir))))	#Keep the files from the directory in a list of files, but only the ones ending with .pdb
		arguments.indir = os.path.abspath(arguments.indir)
	else:		#provided INPUT is not a directory
		raise NameError("ERROR! Incorrect input folder name!")
if not arguments.number_chains:		#checking if the NUMBER OF CHAINS has been provided
	raise NameError("ERROR! The number of chains of the target complex has not been provided! Please, use the help flag, --help, h to see the program instructions!")

os.chdir(arguments.indir + "/../")

## OPTIONAL arguments ##

if arguments.outdir == None:		# Checking if an OUTPUT directory has been provided
	arguments.outdir = arguments.indir + "_output"		# If not, by default it is the INPUT name + _output
else:
	pass
if not os.path.exists(arguments.outdir):		# Checking if the OUTPUT directory (created by us or provided by the user) already exists
	os.mkdir(arguments.outdir)		# If not, create it
	arguments.outdir = os.path.abspath(arguments.outdir)
else:
	arguments.outdir = os.path.abspath(arguments.outdir)

os.chdir(arguments.outdir)		##changes the current directory to OUTPUT directory

### Initializing the LOG system ###

logging.basicConfig(format = '%(levelname)s:%(message)s', filename = arguments.outdir + '/macrocomplex.log', level = logging.DEBUG)
logging.debug('...STARTING...')		# The LOG file is "macrocomplex.log" by default

if arguments.verbose:		# Checking if VERBOSE argument is set
	logging.getLogger().addHandler(logging.StreamHandler())		# If it is set, the LOG file is also printed in STDOUT

if re.search("\d",files[0]):		# Checking if the files
	files = sorted(files, key=lambda x: str("".join([i for i in x if i.isdigit()])))
else:
	files = sorted(files)
#print(files)
#if re.search("\d",files[0]):
#	popped_last = files.pop(-1)
#	files.insert(0,popped_last)
#print(files)
#sys.exit()

logging.info("Parameters used are:\n  - Number of chains: %d\n  - RMSD threshold: %.4f\n  - Clashes threshold %d" % (arguments.number_chains, arguments.rmsd_threshold, arguments.clashes))

### Using the first file as the reference for the macrocomplex ###	

file_path = arguments.indir + "/" + files[0]		#creates path of the first file, which is going to be the reference structure
pdb_parser = Bio.PDB.PDBParser(QUIET = True)	#creation of PDBParser object
ref_structure = pdb_parser.get_structure("reference", file_path)	#creation of the reference structure with the first file, necessary to call the funtion
logging.info("The initial complex has %d chains and are the following:" % (ref_structure[0].__len__()))
for ID in [chain.get_id() for chain in ref_structure[0].get_chains()]:		#loops through all chains of ref_structure
	logging.info("Chain %s", ID)		#prints the ID

# Calling the RECURSIVE FUNCTION for the first time. See DOC for its parameters #
MacrocomplexBuilder(ref_structure = ref_structure, files_list = files, it = 0, not_added = 0, command_arguments = arguments)	#calling the iterative function

### MACROCOMPLEX BUILDING PROCESS FINISHED ###
if len(list(ref_structure[0].get_atoms())) > 99999:
	io = Bio.PDB.MMCIFIO()		#Creates the PDBIO object, which writes a structure object as a PDB file
	io.set_structure(ref_structure[0])		#Sets the reference structure object to be written in a PDB file
	io.save("macrocomplex.cif")	
else: 
	io = Bio.PDB.PDBIO()
	io.set_structure(ref_structure[0])		#Sets the reference structure object to be written in a PDB file
	io.save("macrocomplex.pdb")	
	#the whole macrocomplex gets saved in "macrocomplex.pdb"
stop = timeit.default_timer()
logging.info("The program has finished running! It took %f seconds" % (stop - start))
logging.info("Output files %s saved in %s" %("macrocomplex.pdb and macrocomplex.log",os.path.abspath(arguments.outdir)))