import numpy as np
import networkx as nx
"""
Script to calculate the RMSD of two identical ligands in .pdbqt format after docking with Vina-Carb.
"""
class PDBQT(object):
    def __init__(self, filename):
        """
        Takes a path to a .pdbqt file and creates an object that can be used in RMSD calculations
        :param filename:
        """
        super(PDBQT, self).__init__()
        self.file = open(filename, "r")
        self.file_lines = self.file.readlines()
        self.atoms = {}
        self.energy_result = ""
        self.set_atoms()
        self.graph = nx.Graph()
        self.set_connections()
        self.rings = nx.cycle_basis(self.graph)
        self.ring_hashes_dictionary = {self.ringHash(v): v for v in self.rings}
        self.ring_hashes = self.ring_hashes_dictionary.keys()

        self.negative_functional_groups = {}
        self.negative_functional_group_types = ["S", "S2", "S3", "S4", "S5", "S6", "O62", "O61"]
        self.set_negative_functional_groups()

    def set_negative_functional_groups(self):
        for key, value in self.atoms.items():
            if value[0] in self.negative_functional_group_types:
                self.negative_functional_groups[key] = value[1:]

    def get_energy_result(self):
        return self.energy_result

    def set_atoms(self):
        for line in self.file_lines:
            if line.__contains__("REMARK VINA RESULT: "):
                self.energy_result = float(line.split()[3])

            if line.startswith("ATOM"):
                split = line.split()
                xzy = [x for x in line[30:].split() if x != ""]
                self.atoms[int(split[1])] = [split[2], float(xzy[0]), float(xzy[1]), float(xzy[2])]

    def set_connections(self):
        for atom1 in self.atoms.keys():
            for atom2 in self.atoms.keys():
                if atom1 != atom2 and distance(self.atoms[atom1][1:],
                                               self.atoms[atom2][1:]) < 2 and \
                        not self.atoms[atom1][0].__contains__("H") and \
                        not self.atoms[atom2][0].__contains__("H"):

                    self.graph.add_node(atom1)
                    self.graph.add_node(atom2)
                    self.graph.add_edge(atom1, atom2)

    def get_atoms(self):
        return self.atoms

    def ringHash(self, ring):
        hash = 0
        for atom in ring:
            hash += int(atom) * 3
        return hash

    def get_atom_xyz(self, atom):
        return self.atoms[atom][1:]

    def get_atom_type(self, atom):
        return self.atoms[atom][0]


class RMSD:
    def __init__(self, PDBQT1, PDBQT2):
        self.PDBQT1 = PDBQT1
        self.PDBQT2 = PDBQT2
        self.ring_rmsd = self.do_ring_rmsd
        self.negative_functional_group_rmsd = self.do_negative_functional_group_rmsd()
        self.ring_rmsd = self.do_ring_rmsd()

    def get_ring_rmsd(self):
        return self.ring_rmsd

    def get_negative_functional_group_rmsd(self):
        return self.negative_functional_group_rmsd

    def do_ring_rmsd(self):
        RMSD = 0
        count = 0
        for hash in self.PDBQT1.ring_hashes:
            for atom in self.PDBQT1.ring_hashes_dictionary[hash]:
                RMSD += difSqr(self.PDBQT1.get_atom_xyz(atom), self.PDBQT2.get_atom_xyz(atom))
                count += 1
        if count > 0:
            return np.sqrt(RMSD/count)

    def do_negative_functional_group_rmsd(self):
        RMSD = 0
        count = 0
        for atom in self.PDBQT1.negative_functional_groups.keys():
            RMSD += difSqr(self.PDBQT1.negative_functional_groups[atom], self.PDBQT2.negative_functional_groups[atom])
            count += 1

        if count > 0:
            return np.sqrt(RMSD/count)


def difSqr(predicted, actual):
    y = 0
    for x in range(len(predicted)):
        y += (predicted[x] - actual[x])**2
    return y


def distance(point_1, point_2):
    return np.sqrt((point_1[0] - point_2[0])**2 + (point_1[1] - point_2[1])**2 + (point_1[2] - point_2[2])**2)



