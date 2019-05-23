from Atom import *
import networkx as nx
from PDBQT_to_PDB import *
import os


class PDB(object):
    """docstring for LigandPDB"""

    def __init__(self, filename):
        super(PDB, self).__init__()
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

        for line in self.lines:
            if (line.split()[0].startswith("ATOM") or line.split()[0].startswith("HETATM")) and len(line.strip()) > 20:
                atom = Atom(line)
                if not atom.atom_type.__contains__("H"):
                    self.graph.add_node(atom.getID())
                    a = Atom(line)
                    self.atoms[int(a.getID())] = a
                else:
                    ignore.append(int(atom.getID()))
            elif line.startswith("CONECT"):
                split = line.split()
                for connection in split[2:]:
                    if int(connection) not in ignore and int(split[1]) not in ignore:
                        self.atoms[int(line.split()[1])].add_connection(int(connection))
                        self.graph.add_edge(int(split[1]), int(connection))

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

