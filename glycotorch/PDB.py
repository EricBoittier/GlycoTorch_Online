from glycotorch.Atom import *
import networkx as nx
from glycotorch.PDBQT_to_PDB import *
import os
import numpy as np

class PDB(object):
    """docstring for LigandPDB"""

    def __init__(self, filename):
        super(PDB, self).__init__()
        self.filepath = filename
        self.filename = filename
        self.hash_to_ring_atoms = {}
        self.lines = []
        with open(filename, "r") as f:
            self.lines = f.readlines()
        self.filename = self.filename.split("/")[-1]
        self.filename = self.filename[0:-4]

        self.atoms = {}
        self.ring_names = []
        self.graph = nx.Graph()
        self.setAtomsAndConnections()

        self.connections = None
        self.set_connections()

        self.atoms_np = np.array([[a.x, a.y, a.z] for a
                                  in self.get_atoms()])

        self.find_rings()
        self.setHashToRingAtoms()
        self.nameRings()

    def set_connections(self):
        self.connections = nx.to_dict_of_lists(self.graph)
        self.connections = {k: v for k, v
                            in self.connections.items() if v}

    def find_rings(self):
        self.rings = list(nx.cycle_basis(self.graph))

    def add_edges(self, edges):
        self.graph.add_edges_from(edges)

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

        for line in self.lines:
            # if the line is an atom and has a valid length
            if (line.split()[0].startswith("ATOM") or
                line.split()[0].startswith("HETATM")) \
                    and len(line.strip()) > 20:
                atom = Atom(line)
                atom_numbers.append(atom.id)
                self.graph.add_node(atom.getID())
                a = Atom(line)
                self.atoms[int(a.getID())] = a
                if atom.atom_type == "H":
                   ignore.append(int(atom.getID()))

            #  if the line is a connection
            elif line.startswith("CONECT"):
                split = line.split()
                atom_number = int(split[1])
                if atom_number in atom_numbers:
                    for connection in split[2:]:
                        if int(connection) not in ignore \
                                and atom_number not in ignore:
                            #  add the connection to the atom
                            self.atoms[atom_number].\
                                add_connection(int(connection))
                            self.graph.add_edge(
                                atom_number,
                                int(connection))

    def get_atoms(self) -> list:
        keys = list(self.atoms.keys())
        keys.sort()
        return [self.atoms[k] for k in keys]

    def get_atom_types(self) -> list:
        return [a.atom_type for a in self.get_atoms()]

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
                    (self.ringHash(ring),
                     self.atoms[ring[0]].getLigandType()))
            except KeyError:
                pass

    def get_file_path(self):
        return self.filepath

    def get_file_name(self):
        return self.filename

    def get_ring_name(self, ring):
        for ring_from_list, name in self.ring_names:
            if ring == ring_from_list:
                return name

    def get_ring_name_from_atom(self, atom):
        return self.get_ring_name(self.get_ring(atom))

    def get_ring(self, atom):
        for ring in self.rings:
            if atom in ring:
                return self.ringHash(ring)

    def rename(self, new_name):
        with open(new_name, "w") as file:
            for line in self.lines:
                file.write(line)

