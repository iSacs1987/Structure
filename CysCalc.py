import numpy as np
from pymol import stored, cmd
from scipy.spatial.distance import pdist

def Bond(xyz1, xyz2):
    return np.float64(pdist([xyz1, xyz2], "euclidean"))


def Cys_calc(userSelection):
	stored.SG = []
	userSelection = userSelection + " and n. SG"
	cmd.iterate_state(1, selector.process(userSelection), "stored.SG.append([x,y,z])")
	for elem1 in stored.SG:
		for elem2 in stored.SG:
			if not (elem1[0] == elem2[0] and elem1[1] == elem2[1] and elem1[2] == elem2[2]):
				print(Bond(elem1,elem2))

cmd.extend("Cys_calc",Cys_calc)