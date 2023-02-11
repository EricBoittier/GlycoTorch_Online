from Carbohydrate import *
from Carbohydrate_to_PDBQT import *
from Vina_Input import *
import os
import shutil

for ligand in os.listdir("C:\Users\localadmin\Documents\GitHub\GlycoTorch_Online\docking\pdb_ligands_redo"):
    if not ligand.__contains__(".DS"):
        print(ligand)
        l = ligand.split(".")[0]
        if not os.path.exists("glycotorch_benchmark/"+l):
            os.makedirs("glycotorch_benchmark/"+l)
        c = Carbohydrate("./data/ligands/pdb/"+ligand)
        pdbqt = Carbohydrate_to_PDBQT(c)
        VinaInput(c).makeConfig(path="./glycotorch_benchmark/"+l)
        pdbqt.save_rigid(path="./glycotorch_benchmark/"+l)

        for pdbqt in os.listdir("data/apoproteins/pdbqt/"):
            if pdbqt.split("_")[0] == ligand.split("_")[0]:
                shutil.copy("data/apoproteins/pdbqt/"+pdbqt, "glycotorch_benchmark/"+l)
