import math
import numpy as np
import Bio.PDB as BP
from biopandas.pdb import PandasPdb
from scipy.spatial.distance import pdist


def Bond(xyz1, xyz2):
    return np.float64(
        pdist([list(xyz1.iloc[0, :]), list(xyz2.iloc[0, :])], "euclidean")
    )

def Angle(xyz1,xyz2,xyz3):
    return (360/(2*math.pi))*BP.calc_angle(BP.Vector(xyz1.iloc[0,:]),BP.Vector(xyz2.iloc[0,:]),BP.Vector(xyz3.iloc[0,:]))
	
def Dihedral(xyz1,xyz2,xyz3,xyz4):
    return (360/(2*math.pi))*BP.calc_dihedral(BP.Vector(xyz1.iloc[0,:]),BP.Vector(xyz2.iloc[0,:]),BP.Vector(xyz3.iloc[0,:]),BP.Vector(xyz4.iloc[0,:]))

Bonds = [[], [], []]
Angles = []

ppdb = PandasPdb()
Struc = ppdb.read_pdb("1kim.pdb")
print(Struc.df["ATOM"])
print(Struc.df["OTHERS"])
print(Struc.df["ANISOU"])
print(Struc.df["HETATM"])

Chains = list(set(Struc.df["ATOM"]["chain_id"]))
Chains.sort()
print(Chains)
Numbers = list(
    set(Struc.df["ATOM"][Struc.df["ATOM"]["chain_id"] == Chains[0]]["residue_number"])
)
Numbers.sort()
print(Numbers)

for Nr in Numbers:
    CA = Struc.df["ATOM"][
        (Struc.df["ATOM"]["residue_number"] == Nr)
        & (Struc.df["ATOM"]["atom_name"] == "CA")
        & (Struc.df["ATOM"]["chain_id"] == Chains[0])
    ].loc[:][["x_coord", "y_coord", "z_coord"]]
    print(CA)
    CB = Struc.df["ATOM"][
        (Struc.df["ATOM"]["residue_number"] == Nr)
        & (Struc.df["ATOM"]["atom_name"] == "CB")
        & (Struc.df["ATOM"]["chain_id"] == Chains[0])
    ].loc[:][["x_coord", "y_coord", "z_coord"]]
    C = Struc.df["ATOM"][
        (Struc.df["ATOM"]["residue_number"] == Nr)
        & (Struc.df["ATOM"]["atom_name"] == "C")
        & (Struc.df["ATOM"]["chain_id"] == Chains[0])
    ].loc[:][["x_coord", "y_coord", "z_coord"]]
    O = Struc.df["ATOM"][
        (Struc.df["ATOM"]["residue_number"] == Nr)
        & (Struc.df["ATOM"]["atom_name"] == "O")
        & (Struc.df["ATOM"]["chain_id"] == Chains[0])
    ].loc[:][["x_coord", "y_coord", "z_coord"]]
    if CA.size == 3 and CB.size == 3:
        Distance = Bond(CA, CB)
        Bonds[0].append(Distance)
    if CA.size == 3 and C.size == 3:
        Distance = Bond(CA, C)
        Bonds[1].append(Distance)
    if C.size == 3 and O.size == 3:
        Distance = Bond(C, O)
        Bonds[2].append(Distance)
	if CA.size > 1 and C.size > 1 and O.size > 1:
		Angles.append(Angle(CA,C,O))

for i in range(3):
	print(len(Bonds[i]))
    print(np.mean(Bonds[i]))
    print(np.std(Bonds[i]))
print(len(Angles))
print(np.mean(Angles))
print(np.std(Angles))

NumbersCys = list(set(Struc.df['ATOM'][(Struc.df['ATOM']['residue_name'] == 'CYS') & (Struc.df['ATOM']['atom_name'] == 'SG') & (Struc.df['ATOM']['chain_id'] == Chains[0])]['residue_number']))
print(NumbersCys)