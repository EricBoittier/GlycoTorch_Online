import sys
import os

test_set = os.listdir(r"C:\Users\Eric\Documents\GitHub\GlycoTorch_Online\docking\pdb_redo")

sys.path.insert(1, r'C:\Users\Eric\Documents\GitHub\GlycoTorch_Online')

from Carbohydrate import *

from Protein import *

for pdb_code in test_set:
    p = ProteinPDB(r"C:\Users\Eric\Documents\GitHub\GlycoTorch_Online\docking\pdb_redo\{}".format(pdb_code))
    p.saveLigands(r"C:\Users\Eric\Documents\GitHub\GlycoTorch_Online\docking\ligand_pdb")
