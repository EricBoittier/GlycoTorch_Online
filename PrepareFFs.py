from Carbohydrate import *
from Carbohydrate_to_PDBQT import *

class PrepareFFs(object):
    """docstring for PrepareFFs"""

    def __init__(self, protein, ligand, conf, scale=False):
        super(PrepareFFs, self).__init__()
        self.protein, self.ligand, self.conf = protein, ligand, conf

        self.center_x = None
        self.center_y = None
        self.center_z = None
        self.size_x = None
        self.size_y = None
        self.size_z = None

        self.scale = 0.005
        if scale:
            self.scale = scale

        self.conf_lines = open(self.conf).readlines()

        for line in self.conf_lines:
            if line.__contains__("center_x"):
                self.center_x = float(line.split("=")[-1])
            if line.__contains__("center_y"):
                self.center_y = float(line.split("=")[-1])
            if line.__contains__("center_z"):
                self.center_z = float(line.split("=")[-1])
            if line.__contains__("size_x"):
                self.size_x = float(line.split("=")[-1])
            if line.__contains__("size_y"):
                self.size_y = float(line.split("=")[-1])
            if line.__contains__("size_z"):
                self.size_z = float(line.split("=")[-1])

        self.x_min = self.center_x - self.size_x
        self.y_min = self.center_y - self.size_y
        self.z_min = self.center_z - self.size_z
        self.x_max = self.center_x + self.size_x
        self.y_max = self.center_y + self.size_y
        self.z_max = self.center_z + self.size_z

        self.protein_lines = open(self.protein).readlines()
        self.protein_atom_numbers = []

        self.ligand_lines = open(self.ligand).readlines()
        self.ligand_lines = [line for line in self.ligand_lines if not line.startswith("FF")]
        self.ligand_atom_numbers = []

        for line in self.protein_lines:
            if line.startswith("ATOM"):
                split = line.split()
                resname = split[3]
                atomtype = split[2]
                x = float(split[6])
                y = float(split[7])
                z = float(split[8])

                if self.valid_protein_atom(resname, atomtype, x, y, z):
                    self.protein_atom_numbers.append([resname, split[1]])

        for line in self.ligand_lines:
            if line.startswith("ATOM"):
                split = line.split()
                resname = split[3]
                atomtype = split[2]
                x = split[6]
                y = split[7]
                z = split[8]

                if atomtype.__contains__("@"):
                    self.ligand_atom_numbers.append([atomtype, split[1]])

        print(self.protein_atom_numbers)
        print(self.ligand_atom_numbers)



    def add_FFs(self):
        ffs = []
        for p in self.protein_atom_numbers:
            for l in self.ligand_atom_numbers:
                ffs.append(self.determine_FF(p, l))

        print(ffs)

        outfile = open(self.ligand, "w")
        for ff in ffs:
            outfile.write(ff)
        for line in self.ligand_lines:
            outfile.write(line)

        outfile.close()



    def valid_protein_atom(self, resname, atomtype, x, y, z):
        if (resname == "LYS" and atomtype == "NZ") or \
                (resname == "ARG" and atomtype == "CZ") or \
                (resname == "HIS" and atomtype == "CE1"):
            if (self.x_max > x > self.x_min) and (self.y_max > y > self.y_min) and (self.z_max > z > self.z_min):
                return True

    def determine_FF(self, protein, ligand):
        P_type = protein[0]
        L_type = ligand[0]
        P = protein[1]
        L = ligand[1]
        a = None
        b = None

        if P_type == "ARG":
            if L_type == "@CA":
                a = 9.441
                b = 4.032
            if L_type == "@S":
                a = 5.812
                b = 4.112
            if L_type == "@NS":
                a = 6.253
                b = 4.041
        if P_type == "LYS":
            if L_type == "@CA":
                a = 9.736
                b = 3.210
            if L_type == "@S":
                a = 5.961
                b = 3.441
            if L_type == "@NS":
                a = 10.12
                b = 3.332
        if P_type == "HIS":
            if L_type == "@CA":
                a = 8.723
                b = 3.485
            if L_type == "@S":
                a = 5.534
                b = 3.533
            if L_type == "@NS":
                a = 5.410
                b = 4.098

        return "FFDIST PL LJ {} {} {} {}\n".format(P, L, self.scale * a, b)
