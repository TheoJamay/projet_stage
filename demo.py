import os
import sys
from Bio.PDB import *
import wget

def extraction (pdb_filename) :
	parser = PDBParser()
	io = PDBIO()
	structure = parser.get_structure(pdb_filename, pdb_filename + ".pdb")
	model = structure[0]
	dssp = DSSP(model, pdb_filename + ".pdb", acc_array='Wilke')
	residues = list(dssp)
	with open(pdb_filename + ".txt", "w") as pdb_output :
		for i in range(0, len(list(dssp))) : 
			for j in range(1,3) :
				if j == 1 :
					pdb_output.write(str(residues[i][j] + " "))
				if j == 2 :
					pdb_output.write(str(residues[i][j] + " "))
					if residues[i][j] == "H" or residues[i][j] == "G" \
					or residues[i][j] == "I" :
						pdb_output.write(str("H" + "\n"))
					if residues[i][j] == "E" :
						pdb_output.write(str("E" + "\n"))
					if residues[i][j] == "-" or residues[i][j] == "T" \
					or residues[i][j] == "S" or residues[i][j] == "B" :
						pdb_output.write(str("C" + "\n"))


def dssp_converter (filename_dssp) :
	aa_dssp = []
	s2_8_dssp = []
	s2_3_dssp = []
	with open(filename_dssp+".txt", "r") as output_dssp :
		for line in output_dssp :
			words = line.split()
			aa_dssp.append(words[0])
			s2_8_dssp.append(words[1])
			s2_3_dssp.append(words[2])
	return aa_dssp, s2_8_dssp, s2_3_dssp


def psipred_converter (filename_psipred) :
	aa_psipred = []
	s2_psipred = []
	with open(filename_psipred, "r") as output_psipred :
		for line in output_psipred :
			words = line.split()
			if "AA:" in line :
				aa_psipred.append(words[1])
				aa_psipred = list("".join(aa_psipred))
				
			if "Pred:" in line :
				s2_psipred.append(words[1])
				s2_psipred = list("".join(s2_psipred))
	return aa_psipred, s2_psipred


def q3_calculator(aa_psipred, s2_psipred, aa_dssp, s2_8_dssp, s2_3_dssp) :
	H = 0; I = 0; G = 0; B = 0; E = 0; T = 0; S = 0; C = 0
	pred_ok = len(aa_dssp)
	if len(aa_psipred) != len(aa_dssp) :
		print("Erreur nombre AA")
	for aa in range(0, len(aa_dssp)) :
		if s2_psipred[aa] != s2_3_dssp[aa] :
			if s2_8_dssp[aa] == "H" :
				H += 1
			if s2_8_dssp[aa] == "I" :
				I += 1
			if s2_8_dssp[aa] == "G" :
				G += 1
			if s2_8_dssp[aa] == "B" :
				B += 1
			if s2_8_dssp[aa] == "E" :
				E += 1
			if s2_8_dssp[aa] == "T" :
				T += 1
			if s2_8_dssp[aa] == "S" :
				S += 1
			if s2_8_dssp[aa] == "C" :
				C += 1
			pred_ok -= 1
	print("Q3 = {}\n".format(pred_ok/len(aa_dssp)))
	print("Erreur de prédiction :\nH = {}, I = {}, G = {}, B = {}, E = {},\
 T = {}, S = {}, C = {}\n".format(H, I, G, B, E, T, S, C))


liste = os.listdir('.') 
for i in range(0, len(liste)) :
	if ".pdb" in liste[i] : 
		nom = liste[i]
		str(nom)
		simple = nom[0:4] 
		str(simple)
		
		extraction(simple)

		#url = "https://www.rcsb.org/fasta/entry/"+simple+"/download"
		#filename = wget.download(url)
		#os.system("wget https://www.rcsb.org/fasta/entry/"+simple+"/download \
			#-O "+simple+".fasta")

		aa_dssp, s2_8_dssp, s2_3_dssp = dssp_converter(simple)
		aa_psipred, s2_psipred = psipred_converter("rcsb_pdb_1FKQ.horiz")

		q3_calculator(aa_psipred, s2_psipred, aa_dssp, s2_8_dssp, s2_3_dssp)

		#os.system("../../../../psipred/./runpsipred_single rcsb_pdb_"+simple+  \
		#".fasta")
