import numpy as np
from pymol import stored, cmd
from scipy.spatial.distance import pdist

def Bond(xyz1, xyz2):
    return np.float64(pdist([xyz1, xyz2], "euclidean"))


def CaCalc(userSelection):
	stored.carbons = []
	# userSelection = "1kim"
	userSelection = userSelection + " and n. CA"
	# user Selection = "1kim and n. CA"
	cmd.iterate_state(1, selector.process(userSelection), "stored.carbons.append([x,y,z])")
	for i in range(len(stored.carbons)-1):
		print(Bond(stored.carbons[i],stored.carbons[i+1]))

cmd.extend("CaCalc",CaCalc)