torsion_string = "REMARK       {}    between atoms: {}  and  {} \n"

import os

class PDBQT_Editor(object):
    """docstring for PDBQT_Editor"""

    def __init__(self, path, filename, type_of_file):
        super(PDBQT_Editor, self).__init__()
        self.path = path
        self.filename = filename
        self.file = open(self.path + self.filename, "r")
        self.filelines = self.file.readlines()
        self.torsions = {}
        self.atoms = []
        self.setTorsions()
        if type_of_file == "G":
            self.activateGlycosidicTorsions()
            self.writeNewFile("Glycosidic")

        if type_of_file == "A":
            self.activate_all()
            self.writeNewFile("All")

    def setTorsions(self):
        for line in self.filelines:
            if line.__contains__("REMARK") and line.__contains__("between atoms:"):
                split = line.split()
                atom = split[4]
                neighbour = split[6]
                status = split[1]
                if atom not in self.torsions.keys():
                    self.torsions[atom] = {}
                if neighbour not in self.torsions.keys():
                    self.torsions[neighbour] = {}

                self.torsions[atom][neighbour] = status
                self.torsions[neighbour][atom] = status

    def activateGlycosidicTorsions(self):
        for atom, neighbours in self.torsions.items():
            if atom.__contains__("O") and len(neighbours) == 2:
                keys = list(neighbours.keys())
                if keys[0].__contains__("C") and keys[1].__contains__("C"):
                    for n in self.torsions[atom].items():
                        self.torsions[atom][n[0]] = "A"
                        self.torsions[n[0]][atom] = "A"

    def deactivate_all(self):
        for atom, neighbours in self.torsions.items():
            for n in self.torsions[atom].items():
                self.torsions[atom][n[0]] = "I"
                self.torsions[n[0]][atom] = "I"

    def activate_all(self):
        for atom, neighbours in self.torsions.items():
            for n in self.torsions[atom].items():
                self.torsions[atom][n[0]] = "A"
                self.torsions[n[0]][atom] = "A"


    def writeNewFile(self, type):
        newfile = open(self.path+type+"_"+self.filename, "w")
        for line in self.filelines:
            if line.__contains__("REMARK") and line.__contains__("between atoms:"):
                split = line.split()
                atom = split[4]
                neighbour = split[6]
                status = split[1]

                oldatom = split[4]

                while len(atom) < 5:
                    atom += " "

                linetowrite = torsion_string.format(self.torsions[oldatom][neighbour], atom, neighbour)
                if linetowrite[40] == " ":
                    linetowrite = linetowrite[:40] + linetowrite[(40+1):]

                newfile.write(linetowrite)

            else:
                newfile.write(line)


