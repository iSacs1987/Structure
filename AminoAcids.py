from pymol import cmd,stored

def AminoAcids(userSelection):
	AA = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PRO','PHE','SER','THR','TRP','TYR','VAL']
	for Amino in AA:
		stored.carbons = []
		# userSelection "1kim"
		newSelection = userSelection + " and n. CA and resn " + Amino
		# newSelection "1kim and n. CA and resn ALA"
		# selector.process hands over to PyMol console "select 1kim and n. CA and resn ALA"
		cmd.iterate_state(1, selector.process(newSelection), "stored.carbons.append([x,y,z])")
		print(f"There are {len(stored.carbons)} {Amino} residues in this structure")

cmd.extend("AminoAcids",AminoAcids)