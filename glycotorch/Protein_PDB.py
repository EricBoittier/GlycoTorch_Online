from glycotorch.Atom import *
import networkx as nx
from glycotorch.PDBQT_to_PDB import *
import os



class Protein_PDB(object):
    """docstring for LigandPDB"""

    def __init__(self, filename):
        self.filepath = filename
        self.filename = filename

        # if self.filename.endswith(".pdbqt"):
        #     PDBQT_to_PDB(filename)
        #     self.filename = self.filename + ".pdb"

        self.lines = []
        with open(filename, "r") as f:
            self.lines = f.readlines()

        self.filename = self.filename.split("/")[-1]
        self.filename = self.filename[0:-4]

        self.atoms = {}

        self.ring_names = []

        self.graph = nx.Graph()

        self.setAtomsAndConnections()

        self.connections = nx.to_dict_of_lists(self.graph)
        self.connections = {k: v for k,v in self.connections.items() if v}

        self.rings = nx.cycle_basis(self.graph)

        self.hash_to_ring_atoms = {}
        self.setHashToRingAtoms()

        self.nameRings()

    def toString(self):
        string = ""
        for line in self.lines:
            string += line.replace("\n", r" &#160; ")
        return string

    def setHashToRingAtoms(self):
        for ring in self.rings:
            self.hash_to_ring_atoms[self.ringHash(ring)] = ring

    def getHashToRingAtoms(self):
        return self.hash_to_ring_atoms

    def ringHash(self, ring):
        hash = 0
        for atom in ring:
            hash += int(atom) * 3
        return hash


    def getRings(self):
        return self.rings

    def getRingNames(self):
        return self.ring_names


    def getNeighbours(self, atomID):
        self.graph.__getitem__(atomID)

    def setAtomsAndConnections(self):

        ignore = []
        atom_numbers = []
        Hs = []

        for line in self.lines:
            if line.split()[0].startswith("ATOM") and len(line.strip()) > 20:
                atom = Atom(line)
                atom_numbers.append(int(atom.id))

                self.graph.add_node(atom.getID())
                a = Atom(line)

                self.atoms[int(a.getID())] = a

                if a.atomname.__contains__("H"):
                    self.atoms[int(a.getID())].atomtype = "HD"
                    Hs.append(a)


        for key in self.atoms.keys():
            atom = self.atoms[key]
            if atom.atomname == "O":
                for H in Hs:
                    if distance(atom.getXYZ(), H.getXYZ()) < 1.5:
                        print("or here")
                        self.atoms[int(H.getID())].atomname = "HD"
                        self.atoms[int(atom.getID())].atomname = "OA"

            if atom.atomname == "N":
                for H in Hs:
                    if distance(atom.getXYZ(), H.getXYZ()) < 1.5:
                        print("here")
                        self.atoms[int(H.getID())].atomname = "HD"


    def distance(self, point_1, point_2):
        return ((point_1[0] - point_2[0]) ** 2 + (point_1[1] - point_2[1]) ** 2 + (
                    point_1[2] - point_2[2]) ** 2)**0.5

    def getAllConnections(self):
        return self.connections

    def getConnections(self, atomID):
        return self.connections[atomID]

    def getFilename(self):
        return self.filename

    def getAtoms(self):
        return self.atoms

    def getLines(self):
        return self.lines

    def getGraph(self):
        return self.graph

    def nameRings(self):
        for ring in self.rings:
            try:
                self.ring_names.append(
                    (self.ringHash(ring), self.atoms[ring[0]].getLigandType()))
            except KeyError:
                pass

    def getFilepath(self):
        return self.filepath

    def getFilename(self):
        return self.filename

    def getRingName(self, ring):
        for ring_from_list, name in self.ring_names:
            if ring == ring_from_list:
                return name

    def getRingNamefromAtom(self, atom):
        return self.getRingName(self.getRing(atom))

    def getRing(self, atom):  # TODO: give this a sensible name
        for ring in self.rings:
            if atom in ring:
                return self.ringHash(ring)

    def rename(self, new_name):
        with open(new_name, "w") as file:
            for line in self.lines:
                file.write(line)

    def write_pdbqt_file(self, path, filename):
        file = open(os.path.join(path, filename), "w")
        file.write(self.get_pdbqt_string())
        file.close()


    def get_pdbqt_string(self):

        pdbqt_string = ""
        for atom in self.atoms.values():
            line = self.get_pdbqt_line(atom)
            if line:
                pdbqt_string += line
        return pdbqt_string



    def get_pdbqt_line(self, atom):

        write_atom = True

        self.c = int(atom.id)
        id = self.pad_before(self.c, 7)
        name = self.pad_before(atom.atom_type, 4) + " "
        ligand_type = self.pad_before(atom.ligand_type, 4)
        chain = self.pad_before(atom.chain, 2)
        res_id = self.pad_before(atom.ligandID, 4)
        x = "{:7.3f}".format(atom.x)
        x = self.pad_before(x, 12)
        y = "{:7.3f}".format(atom.y)
        y = self.pad_after(y, 8)
        z = "{:7.3f}".format(atom.z)


        ADT = atom.atomname

        if ADT == "H":
            return False

        # if atom.atomname == "O" and len(atom.connections) == 2:
        #     if not self.atoms[atom.connections[0]].atomname.__contains__("H") \
        #             or not self.atoms[atom.connections[0]].atomname.__contains__("H"):
        #         ADT = "O"
        #     else:
        #         ADT = "OA"
        #
        # if atom.atom_type.__contains__("N"):
        #     if len(atom.connections) == 2:
        #         add_H = True
        #         ADT = "N"
        #     elif len(atom.connections) < 4:
        #         ADT = "NA"
        #     else:
        #         ADT = "N"
        #
        #
        # if atom.atomname == "H" and len(atom.connections) == 1:
        #     if self.atoms[atom.connections[0]].atomname.__contains__("O") \
        #             or self.atoms[atom.connections[0]].atomname.__contains__("N"):
        #         ADT = "HD"
        #     else:
        #         write_atom = False

        line = "ATOM{}{}{}{}{}{} {} {}  0.00  0.00     0.000 {}\n".format(id, name, ligand_type, chain, res_id,
                                                                          x, y, z, ADT)

        return line


    def pad_after(self, string, pad):
        string = str(string)
        while len(string) < pad:
            string = string + " "
        return string

    def pad_before(self, string, pad):
        string = str(string)
        while len(string) < pad:
            string = " " + string
        return string

# P_PDB = Protein_PDB(r"C:\Users\localadmin\Desktop\2axm.pdb")
# P_PDB.write_pdbqt_file(r"C:\Users\localadmin\Desktop\\", "test_protein.pdbqt")