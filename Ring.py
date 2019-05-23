from Cremer_Pople import *
from Dihedral import *
from PDB_Abbreviation_to_Sugar_Name import *
from Get_IUPAC_Ring_Conformation import *

class Ring(object):
    """
    A class to handle ring operations

    """

    def __init__(self, c1, c2, c3, c4, c5, o5, atoms):
        self.atoms = atoms
        self.residue = c1.getLigandType()
        self.ID = c1.getLigandID()

        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4
        self.c5 = c5
        self.o5 = o5

        self.ring = [c1.id, c2.id, c3.id, c4.id, c5.id, o5.id]

        self.ring_hash = self.get_ring_hash()

        self.c1_xyz = c1.getXYZ()
        self.c2_xyz = c2.getXYZ()
        self.c3_xyz = c3.getXYZ()
        self.c4_xyz = c4.getXYZ()
        self.c5_xyz = c5.getXYZ()
        self.o5_xyz = o5.getXYZ()

        self.iupac_name = get_iupac_ring(self.c1_xyz, self.c2_xyz, self.c3_xyz,
                                            self.c4_xyz, self.c5_xyz, self.o5_xyz)

        self.ring_name = abbreviationToSugarName(self.residue)

        self.cremer_pople= cremer_pople([self.o5_xyz, self.c1_xyz, self.c2_xyz, self.c3_xyz, self.c4_xyz, self.c5_xyz])

        self.phi = self.cremer_pople["phi"]
        self.theta = self.cremer_pople["theta"]
        self.Q = self.cremer_pople["Q"]
        self.q2 = self.cremer_pople["q_2"]
        self.q3 = self.cremer_pople["q_3"]
        self.phi_radians = self.cremer_pople["phi_radians"]
        self.theta_radians = self.cremer_pople["theta_radians"]

        self.c1_functional_group = self.c2_functional_group = \
            self.c3_functional_group = self.c4_functional_group = self.c5_functional_group = []

    def get_functional_groups(self, glycosidic_atoms):
        self.c1_functional_group = self.get_functional_group(self.c1, glycosidic_atoms)
        self.c2_functional_group = self.get_functional_group(self.c2, glycosidic_atoms)
        self.c3_functional_group = self.get_functional_group(self.c3, glycosidic_atoms)
        self.c4_functional_group = self.get_functional_group(self.c4, glycosidic_atoms)
        self.c5_functional_group = self.get_functional_group(self.c5, glycosidic_atoms)

        # print(self.c2_functional_group)


    def get_functional_group(self, atom, glycosidic_atoms):
        functional_group = []

        for x in atom.connections:
            if x not in self.ring and x not in glycosidic_atoms:
                functional_group.append(x)



        if len(functional_group) > 0:

            visited, stack = set(), [functional_group[0]]

            while stack:
                vertex = stack.pop()
                if vertex not in visited:
                    visited.add(vertex)

                    fg = set()
                    for x in self.atoms[vertex].connections:
                        if x not in self.ring and x not in glycosidic_atoms and \
                                self.atoms[x].ligandID == atom.ligandID:
                            fg.add(x)

                    stack.extend(fg - visited)

            functional_group = visited

        return list(functional_group)


    def get_cremer_pople(self):
        return self.cremer_pople

    def get_ring_name(self):
        return self.ring_name

    def get_ID(self):
        return self.ID

    def get_phi(self):
        return self.phi

    def get_theta(self):
        return self.theta

    def get_ring_hash(self):
        hash = self.c1.getID()*3 + \
               self.c2.getID()*3 + self.c3.getID()*3 + self.c4.getID()*3\
               + self.c5.getID()*3 + self.o5.getID()*3
        return hash

    def get_pdb_lines(self):
        print(self.c1.line)
        print(self.c2.line)
        print(self.c3.line)
        print(self.c4.line)
        print(self.c5.line)
        print(self.o5.line)
