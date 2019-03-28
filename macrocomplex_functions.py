import Bio.PDB
import sys
import string
import os
import argparse
import timeit
import logging
import re

def Key_atom_retriever(chain):
	"""This function retrieves the key atom, CA in case of proteins and C4' in case of nucleic acids, to do the superimposition and also returns a
	variable indicating the kind of molecule that that chain is: either DNA, RNA or PROTEIN

	Arguments:

	chain (Bio.PDB.Chain.Chain): an instance of class chain

	Returns:
	
	atoms (list): contains all key atoms (CA/C4') instances

	molecule (str): contains type of molecule of the chain

	"""
	### Declaring and creating new variables ###
	nucleic_acids = ['DA','DT','DC','DG','DI','A','U','C','G','I']		#creating a list with all possible nucleic acids letters
	RNA = ['A','U','C','G','I']											#creating a list with all possible RNA letters
	DNA = ['DA','DT','DC','DG','DI']									#creating a list with all possible DNA letters
	atoms = []
	### Loops through all residues of the chain ###
	for res in chain:		
		res_name = res.get_resname()[0:3].strip()		#get the name of the residue (with no spaces)
		## Appends the CA atoms and sets the molecule type of protein ##
		if res.get_id()[0] == " " and res_name not in nucleic_acids:		#checks whether the residue is not a HETATM or nucleic acid
			if 'CA' not in res:		#checks whether the residue has CA atoms
				logging.warning("This protein residue %d %s does not have CA atom" % (res.get_id()[1], res_name))
			else:
				atoms.append(res['CA'])		#append CA atoms to the list of sample atoms
				molecule = 'PROTEIN'		#set the molecule type to protein
		## Append the C4 atoms and sets the molecule type of DNA or RNA ##
		elif res.get_id()[0] == " " and res_name in nucleic_acids:		#checks whether the residue is a nucleic acid and not HETATM
			if res_name in DNA:			#checks whether the residue is a DNA nucleotide
				molecule = 'DNA'		#set the molecule type to DNA
			elif res_name in RNA:		#checks whether the residue is a RNA nucleotide
				molecule = 'RNA'		#set the molecule type to RNA
			atoms.append(res['C4\''])	#append C4' atoms to the list of atoms
	return(atoms, molecule)		# Return all key atoms list and the type of molecule to which they belong

def ID_creator(IDs, ID):
	"""This function creates IDs"""
	UP = list(string.ascii_uppercase)	
	LOW = list(string.ascii_lowercase)
	DIG = list(string.digits)
	alphabet = UP + LOW + DIG		#creates an alphabet containing all the possible characters that can be used as chain IDs

	if len(IDs) < 62:
		if ID not in IDs:
			return ID
		if ID in IDs:
			for i in range(0, len(alphabet)):
				if alphabet[i] not in IDs:
					return alphabet[i]
				else:
					continue
	elif len(IDs) >= 62:
		for char1 in alphabet:
			for char2 in alphabet:
				ID = char1 + char2
				if ID not in IDs:
					return ID
				else:
					continue

def superimposition(ref_structure, sample_structure, rmsd_threshold):
	"""This function, given a reference and a sample structure does the superimposition of every combination of pairs of chains and calculates the RMSD.
	It returns a dictionary with a tuple of the reference and sample chain as a tuple and the superimposer instance resulting from those two chains, as
	well as two variables, the ID of the chain with the smallest RMSD when superimposing with the reference structure and the RMSD itself

	Arguments:

	ref_structure (Bio.PDB.Structure): is the structure on which the macrocomplex is gonna get build on every iteration of the function

	sample_structure (Bio.PDB.Structure): is the structure that is gonna be added on every iteration of the function

	Returns:

	all_superimpositions (dict): dictionary of tuples of chain identifiers as key and Superimposer instances as values

	superimposed_chains (boolean): set to True if there has been at least one superimposition, otherwise is False.

	best_RMSD (float): RMSD of the best superimposition (the lowest RMSD value)

	"""
	### Saving arguments passed on to the function ###
	ref_model = ref_structure[0]			#retrieves the first and only available model of reference structure
	sample_model = sample_structure[0]		#retrieves the first and only available model of the sample structure 
	### Initializing and declaring variables ###
	best_sample_chain_ID = best_ref_chain_ID = ""
	best_RMSD = 0					#variable for the lowest RMSD
	prev_RMSD = True				#variable to know we are in the first combination of pairs of chains					
	superimposed_chains = False		#variable that indicates the presence of a superimposed chain (True if there is superimposed chain)
	all_superimpositions = {}		#start the dictionary that will contain all superimposition instances
	### Superimposition of every combination of pairs of chains between the reference and the sample structures ###
	## loops through all chains in the reference model ##
	for ref_chain in ref_model:		
		logging.info("Processing reference chain %s", ref_chain.id)
		ref_atoms, ref_molecule = Key_atom_retriever(ref_chain)		#Retrieves all key atoms (CA or C4') and molecule type of the sample						
		## loops through all chains in the sample model ##
		for sample_chain in sample_model:							
			logging.info("Processing sample chain %s", sample_chain.id)
			sample_atoms, sample_molecule = Key_atom_retriever(sample_chain)		#Retrieves all key atoms (CA or C4') and molecule type of the sample
			if ref_molecule != sample_molecule:				#checks that the molecular types of ref chain and sample chain are the same
				logging.warning("Cannot superimpose. Reference chain %s is %s and sample chain %s is %s" %(ref_chain.get_id(), ref_molecule, sample_chain.get_id(), sample_molecule))
			elif len(ref_atoms) != len(sample_atoms):		#checks that the length of ref_atoms and sample_atoms is the same
				logging.warning("Cannot superimpose. The number of atoms of the reference chain %s is %d and the number of atoms of the sample chain %s is %d", ref_chain.get_id(), len(ref_atoms), sample_chain.get_id(), len(sample_atoms))
			## Make the superimposition between reference and sample chain ##
			else:		#everything is fine, same type of molecule, same length of atom lists
				super_imposer = Bio.PDB.Superimposer()				#creates superimposer instance
				super_imposer.set_atoms(ref_atoms, sample_atoms)	#creates ROTATION and TRANSLATION matrices from lists of atoms to align
				RMSD = super_imposer.rms 							#retrieves RMSD
				if RMSD > rmsd_threshold:
					logging.info("The RMSD between chain %s of the reference and chain %s of the sample is %f", ref_chain.id, sample_chain.id, RMSD)
					continue
				if prev_RMSD is True or RMSD < prev_RMSD:			#checks that the RMSD of this combination is smaller than the previous one
					best_sample_chain_ID = sample_chain.id 		
					best_ref_chain_ID = ref_chain.id 				#with this condition, the superimposer instance and other important
					best_RMSD = RMSD 								#information pertaining to the superimposition with the smallest
					prev_RMSD = RMSD 								#RMSD will be saved
				all_superimpositions[(ref_chain.id, sample_chain.id)] = super_imposer		#saving ALL superimposer instances in a dictionary
				superimposed_chains = True 							# The superimposition has been made
				logging.info("The RMSD between chain %s of the reference and chain %s of the sample is %f", ref_chain.id, sample_chain.id, RMSD)
	### checks that there has been, at least, one superimposition ###
	if superimposed_chains is True:									
		all_superimpositions = sorted(all_superimpositions.items(), key=lambda k:k[1].rms)		#sorting by the lowest RMSD and saving to a list
		logging.info("The combination of chains with the lowest RMSD is ref chain %s and sample chain %s with an RMSD of %f", best_ref_chain_ID, best_sample_chain_ID, best_RMSD)
	return(all_superimpositions, superimposed_chains, best_RMSD)

def MacrocomplexBuilder(ref_structure, files_list, it, not_added, command_arguments):
	"""This recursive function superimposes the most similar chain of a binary interaction PDB file with a reference structure and adds the transformed chain to the building complex

	Arguments:

	ref_structure (Bio.PDB.Structure): is the structure on which the macrocomplex is gonna get build on every iteration of the function

	files_list (list): a list containing all the pdb files of binary interactions between the different subunits or chains that form the complex

	command_arguments(argparse object): is the object containing all the command-line arguments. Contains:

			RMSD (float): this is the RMSD threshold. If the RMSD of a superimposition between reference and sample structures is greater than this value, it will be
			considered as a wrong superimposition and will not be used to build the complex

			clashes (int): this is the clashes or contacts threshold. If the number of contacts between two chains exceeds this value, the superimposition will not be
			taken into account for the rotated chain is either present in the complex already or clashing with other chains because it should not be there. The chain
			in question will be dismissed

			number_chains (int): this is the numbers of chains that the complex must have in order to stop running. However, if these value is never reached, the program
			will stop after a certain number of iterations
			
			it (int): this is a counter that keeps track of the interaction of the iterative function
			
			indir(str): this is the input directory relative path

			outdir(str): this is the output directory relative path

			iterations(boolean): this is set True if the user wants a pdb file for each iteration of the complex. Otherwise is False


	It is an iterative function, it calls itself until certain condition is met, then:

	Returns:

	ref_structure (Bio.PDB.Structure): pdb structure instance containing all chains of the final macrocomplex.

	"""

	### Saving arguments passed on to the function ###
	i = it 															#number of iterations
	n = not_added													#number of files that have been parsed but no chain has been added
	nc = command_arguments.number_chains							#number of chains		
	clashes_threshold = command_arguments.clashes 					#clashes threshold
	RMSD_threshold = command_arguments.rmsd_threshold 				#RMSD threshold
	indir = command_arguments.indir 								#input directory relative path
	outdir = command_arguments.outdir 								#output directory relative path
	pdb_iterations = command_arguments.pdb_iterations 				#if True, each iteration is stored in a pdb file
	iterations = command_arguments.it 								#maximum number of iterations

	alphabet = list(string.ascii_uppercase)	+ list(string.ascii_lowercase) + list(string.digits)		#creates an alphabet containing all the possible characters that can be used as chain IDs 

	chains = ref_structure[0].__len__()
	### Prints the current iteration and number of chains of the current complex ###
	logging.info("This is the iteration #%d of the recursive function" % i )
	logging.info("The complex has %d chains at this point" % chains)
	print ("The complex has %d chains at this point" % chains )

	### Checks if the current macrocomplex satisfies the desired number of chains or just stops at iteration 150 ### 
	if chains == nc: 
		logging.info("The whole macrocomplex has been successfully build with the desired number of chains")
		logging.info("The final complex has %d chains" % chains)
		print ("The end, we reached %d chains" % nc)
		return 	ref_structure			#END OF THE RECURSIVE FUNCTION
	elif n > len(files_list):
		logging.info("The whole macrocomplex has been build")
		logging.info("The final complex has %d chains, not %d, as requested" % (chains, nc))
		logging.info("We have arrived to iteration %d" %(i))
		print("The end, we cannot add any more chains, it has %d chains" % chains)
		return ref_structure		#END OF THE RECURSIVE FUNCTION

	### Selects the file to analyze in this iteration. It is always the first element of the list of files because once analyzed it is substracted and appended at the end of the list ###
	sample = files_list[0]		#saves the first file name of the list of files as the sample
	logging.info("We are processing the file %s" % (sample))
	file_path = indir + "/" + sample 			#takes the path of the sample file
	pdb_parser = Bio.PDB.PDBParser(QUIET = True)		#parses the sample PDB file and creates a sample PDBParser object
	sample_structure = pdb_parser.get_structure("sample", file_path)		#saves the Structure object of the sample PDBParser object
	sample_model = sample_structure[0]		#obtains the first and only available model of the sample structure 

	### Calling the superimposition function to obtain the superimposition of every combination of pairs of chains between the reference and sample structures
	all_superimpositions, superimposed_chains, best_RMSD = superimposition(ref_structure, sample_structure, RMSD_threshold)

	### There are no superimposed chains or RMSD is above the threshold --> Call again the recursive function ###
	if superimposed_chains is False or best_RMSD > RMSD_threshold:		#if condition is met, there are no superimposed chains, or the RMSD is not small enough to be considered
		file = files_list.pop(0)		#substracts the current file
		files_list.append(file)			#and adds it at the end of the list of files
		i += 1							#calling again the recursive function to analyze the next file
		n += 1
		return MacrocomplexBuilder(ref_structure = ref_structure, files_list = files_list, it = i, not_added = n, command_arguments = command_arguments)	#call again the iterative function, j does not change
	### There are superimposed chains ###
	else:
		## Loops through the superimposition dictionary, obtaining the superimposition instances and the reference and sample IDs ##
		for chains, sup in all_superimpositions:
			logging.info("We are processing the superimposition of ref chain %s with sample chain %s with an RMSD of %f" % (chains[0],chains[1], sup.rms))
			if sup.rms > RMSD_threshold:			#Checks that the superimposition has an RMSD above the threshold
				logging.info("This superimposition of ref chain %s with sample chain %s has an RMSD bigger than the threshold, therefore it is skipped" % (chains[0],chains[1]))
				continue							#if not, skip that superimposition
			sup.apply(sample_model.get_atoms())		#applies ROTATION and TRANSLATION matrices to all the atoms in the sample model
			## Gets the sample chain that was not superimposed with the reference chain --> putative chain to add ##
			chain_to_add = [chain for chain in sample_model.get_chains() if chain.get_id() != chains[0]][0]		
			present_chain = False		#this variable indicates whether the chain to add is present on the building complex or not: False => not present, True => present
			sample_atoms, sample_molecule = Key_atom_retriever(chain_to_add)		#retrieves all key atoms (CA or C4') and molecule type of chain_to_add
			logging.info("Putative chain to add is %s" % chain_to_add.id)
			## Loops through all the chains from the reference structure ##
			all_atoms = []
			for chain in ref_structure[0].get_chains():					
				ref_atoms, ref_molecule = Key_atom_retriever(chain)		#retrieves all key atoms (CA or C4') and molecule type of the reference present chain
				## Makes a Neighbor Search to look for clashes between the chain to add and the chains from the reference structure ##
				all_atoms.extend(ref_atoms)
				Neighbor = Bio.PDB.NeighborSearch(ref_atoms)			#creates an instance of class NeighborSearch, given a list of reference atoms 
				clashes = []		#declares a list that will contain all the atoms that clash between the reference and sample chains
				for atom in sample_atoms:								#loops through the list of atoms of chain_to_add
					atoms_clashed = Neighbor.search(atom.coord,5)		#produces a Neighbor search that returns all atoms/residues/chains/models/structures that have at least one atom within radius of center. 
					if len(atoms_clashed) > 0:				#if there are clashes
						clashes.extend(atoms_clashed)		#adds the atoms list to the list of clashes
				if len(clashes) > clashes_threshold:		#checks that the number of total clashes is above the threshold
					present_chain = True					#then, chain_to_add is considered a chain already present in the complex
					logging.info("The number of clashes between the chain to add %s and reference chain %s is %d, therefore the chain is the same and it is skipped" % (chain_to_add.id, chain.id,len(clashes)))
					break 									#skips continuing through the loop, as it already clashes with one reference chain
				## Checks that the number of total clashes is under the threshold ##
				elif len(clashes) <= clashes_threshold:		
					logging.info("The number of clashes between the chain to add %s and reference chain %s is %d, it is under the threshold" % (chain_to_add.id, chain.id,len(clashes)))
					continue								#continue the loops, as we must ensure that chain_to_add does not clash with ANY reference chain
			#Neighbor2 = Bio.PDB.NeighborSearch(all_atoms)
			#clashes2 = []
			#for atom2 in sample_atoms:								#loops through the list of atoms of chain_to_add
			#	atoms_clashed2 = Neighbor2.search(atom.coord, 65)		#produces a Neighbor search that returns all atoms/residues/chains/models/structures that have at least one atom within radius of center. 
			#	if len(atoms_clashed2) > 0:				#if there are clashes
			#		clashes2.extend(atoms_clashed2)		#adds the atoms list to the list of clashes
			#if len(clashes2) > 100:		#checks that the number of total clashes is above the threshold
			#	logging.info("The number of clashes between the chain to add %s and rhe reference structure is %d, therefore the chain is NOT too far away and is added" % (chain_to_add.id,len(clashes2)))
				#logging.info("The number of clashes between the chain to add %s and rhe reference structure is %d, therefore the chain is too far away and is not added" % (chain_to_add.id,len(clashes2)))
				#continue
			#elif len(clashes2) < 100:		#checks that the number of total clashes is above the threshold
			#	present_chain = True					#then, chain_to_add is considered a chain already present in the complex
			#	logging.info("The number of clashes between the chain to add %s and rhe reference structure is %d, therefore the chain is too far away and is not added" % (chain_to_add.id,len(clashes2)))
				#logging.info("The number of clashes between the chain to add %s and rhe reference structure is %d, therefore the chain is too far away and is not added" % (chain_to_add.id,len(clashes2)))
			#	continue 	
			## Rotated chain to add is not a chain already in the building macrocomplex structure, then adds it, with its original ID or with a new one ##
			if present_chain is False:						
				logging.info("Chain %s superimposed with chain %s yields rotated chain %s which is not in the complex" %(chains[0],chains[1],chain_to_add.id))
				chain_ids = [chain.id for chain in ref_structure[0].get_chains()]	#list containing IDs of all chains present in reference structure
				ID = ID_creator(chain_ids, chain_to_add.id)
				chain_to_add.id = ID
				ref_structure[0].add(chain_to_add)	#adds chain_to_add to the building macrocomplex structure
				logging.info("Added Chain %s" % ID)
				print("Added Chain %s" % ID)
				## Checks whether the user provided the iterations argument, then save each iteration of the current complex in a PDB file ##
				if pdb_iterations:		
					io = Bio.PDB.MMCIFIO()							#creates a PDBIO instance
					io.set_structure(ref_structure[0])				#assigns the reference structure to the IO object
					io.save("macrocomplex_chains_%d.pdb" %(ref_structure[0].__len__()))	#saves the structure on a file
					logging.info("saving macrocomplex_chains_%d.pdb in %s" %(ref_structure[0].__len__(),outdir))
				#break							#breaks loop
				file = files_list.pop(0)		#substracts the first file of the files list
				files_list.append(file)			#adds the file at the end of the files list
				i += 1							#adds one to the iteration variable
				n = 0
				#this is what makes the function recursive, it calls itself on the return, executing the whole function again and again until certain condition is met
				return MacrocomplexBuilder(ref_structure = ref_structure, files_list = files_list, it = i, not_added = n, command_arguments = command_arguments)
	### Once the current file has been analyzed it is substracted and appended at the end of the files list ###
	file = files_list.pop(0)		#substracts the first file of the files list
	files_list.append(file)			#adds the file at the end of the files list
	i += 1							#adds one to the iteration variable
	n += 1
	#this is what makes the function recursive, it calls itself on the return, executing the whole function again and again until certain condition is met
	return MacrocomplexBuilder(ref_structure = ref_structure, files_list = files_list, it = i, not_added = n, command_arguments = command_arguments)		
