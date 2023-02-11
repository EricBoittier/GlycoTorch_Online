from PDB import *
from PDB_Abbreviation_to_Sugar_Name import *
import networkx as nx
from Dihedral import dihedral
from Ring import *
from LinkageError import *
from Linkage import *
from Geometry import *
from RingError import *

class Carbohydrate(PDB):
    def __init__(self, filename):
        super().__init__(filename)

        print(self.atoms_np)
        print(self.graph)
        #  check if the graph has edges and connections
        if len(self.graph.edges) == 0 or \
                (True in [0 == len(v.connections)
                          for v in self.atoms.values()]):
            print("No edges found in PDB file, "
                  "adding edges from distance matrix")
            connections = self.connections_from_atoms()
            connections = [(self.get_atoms()[k].id,
                            self.get_atoms()[v].id)
                           for k, v in connections]
            for a, b in connections:
                self.atoms[a].add_connection(b)
                self.atoms[b].add_connection(a)

            self.add_edges(connections)
            self.set_connections()
            self.find_rings()

        if len(self.rings) == 0:
            raise RingError("No rings found in PDB file")

        self.carbohydrate_name = ""
        self.SugarGraph = nx.DiGraph()
        self.Linkages = {}
        self.Rings = {}
        self.glycosidic_oxygens = []
        self.ordered_ring_hashes = []
        self.ordered_linkages = []

        self.assign_ring_positions()
        self.assign_glycosidic_atoms()
        self.make_sugar_graph()
        self.set_name()

    def connections_from_atoms(self):
        edges = get_edges_from_distance_matrix(
            self.atoms_np,
            self.get_atom_types())
        return edges

    def assign_ring_positions(self):
        """Assign the ring positions to the atoms in the ring based on
        the usual Carbohydrate nomenclature.

        :return: None
        """
        for ring in self.rings:
            # should have 5 carbons and 1 oxygen
            c1 = c2 = c3 = c4 = c5 = o5 = False
            for atom in ring:
                # the ring oxygen is not a glycosidic oxygen
                # so we can assign it to the ring
                # and set the axial property for glycosidic oxygens
                if not self.isO5(atom):
                    self.atoms[atom].set_axial(self.is_axial(atom))
                #  check if the atom is a c1
                if not c1 and self.isC1(atom):
                    self.atoms[atom].set_sugar_position("C1")
                    c1 = self.atoms[atom]
                elif not c2 and self.isC2(atom):
                    self.atoms[atom].set_sugar_position("C2")
                    c2 = self.atoms[atom]
                elif not c3 and self.isC3(atom):
                    self.atoms[atom].set_sugar_position("C3")
                    c3 = self.atoms[atom]
                elif not c4 and self.isC4(atom):
                    self.atoms[atom].set_sugar_position("C4")
                    c4 = self.atoms[atom]
                elif not c5 and self.isC5(atom):
                    self.atoms[atom].set_sugar_position("C5")
                    c5 = self.atoms[atom]
                elif not o5 and self.isO5(atom):
                    self.atoms[atom].set_sugar_position("O5")
                    o5 = self.atoms[atom]
                else:
                    raise RingError("Something went wrong with" + \
                                    " assigning ring positions for\n "
                                    + f"{atom, self.atoms[atom]}\n" + \
                                    " in ring" + f"{ring}, it was not " + \
                                    f"able to be assigned to a ring position.")

            if c1 and c2 and c3 and c4 and c5 and o5:
                self.Rings[self.ringHash(ring)] = \
                    Ring(c1, c2, c3, c4, c5, o5, self.atoms)
            else:
                raise RingError("Not all atoms in ring")

    def get_ordered_hashes(self):
        return self.ordered_ring_hashes

    def assign_glycosidic_atoms(self):
        for atom in self.atoms.values():
            if self.isGylcosidic(atom.getID()):
                self.atoms[atom.getID()].set_sugar_position = "GO"
                self.glycosidic_oxygens.append(atom)

    def make_sugar_graph(self):
        glycosidic_atoms = [x.id for x in self.glycosidic_oxygens]

        for number, atom in enumerate(self.glycosidic_oxygens):
            neighbours = atom.get_connections()

            a1 = self.atoms[neighbours[0]]
            a2 = self.atoms[neighbours[1]]

            if a1.isC1():
                C1 = a1
                CX = a2
            elif a2.isC1():
                C1 = a2
                CX = a1
            else:
                raise LinkageError("C1 not found")

            O5 = False
            for neighbour in C1.get_connections():
                if self.atoms[neighbour].isO5():
                    O5 = self.atoms[neighbour]
            if not O5:
                raise LinkageError("O5 not found")

            if CX.isC2():
                CX_minus_1 = False
                for neighbour in CX.get_connections():
                    if self.atoms[neighbour].isC1():
                        CX_minus_1 = self.atoms[neighbour]
                if not CX_minus_1:
                    raise LinkageError("C1 not found when finding CX-1")
            elif CX.isC3():
                CX_minus_1 = False
                for neighbour in CX.get_connections():
                    if self.atoms[neighbour].isC2():
                        CX_minus_1 = self.atoms[neighbour]
                if not CX_minus_1:
                    raise LinkageError("C2 not found when finding CX-1")
            elif CX.isC4():
                CX_minus_1 = False
                for neighbour in CX.get_connections():
                    if self.atoms[neighbour].isC3():
                        CX_minus_1 = self.atoms[neighbour]
                if not CX_minus_1:
                    raise LinkageError("C3 not found when finding CX-1")
            else:  # must be C6 by order of elimination
                CX_minus_1 = False
                for neighbour in CX.get_connections():
                    if self.atoms[neighbour].isC5():
                        CX_minus_1 = self.atoms[neighbour]
                if not CX_minus_1:
                    raise LinkageError("C5 not found when finding CX-1")

            ring1 = self.get_ring(C1.getID())
            ring2 = self.get_ring(CX.getID())

            self.Rings[ring1].set_functional_groups(glycosidic_atoms)
            self.Rings[ring2].set_functional_groups(glycosidic_atoms)

            self.SugarGraph.add_node(ring1)
            self.SugarGraph.add_node(ring2)
            self.SugarGraph.add_edges_from([(ring1, ring2)])

            self.atoms[atom.id].atom_type = "O" + CX.sugar_position[-1]

            temp_linkage = Linkage(O5, C1, atom, CX, CX_minus_1)
            temp_linkage.set_rings(self.Rings[ring1], self.Rings[ring2])

            self.Linkages[ring1] = temp_linkage

    def get_ring_object(self, ringHash):
        return self.Rings[ringHash]

    def get_file_name(self):
        return self.filename

    def is_axial(self, atom):
        test = self.atoms[atom].getXYZ()
        neighbour_not_in_ring = ""
        in_ring = []
        in_ring_label = []
        in_rings_neighbour = []

        if len(self.connections[atom]) < 2:
            return None
        try:
            for neighbour in self.connections[atom]:
                if not self.isAtomInRing(neighbour):
                    neighbour_not_in_ring = self.atoms[neighbour].getXYZ()

                elif self.isAtomInRing(neighbour):
                    in_ring_label.append(neighbour)
                    in_ring.append(self.atoms[neighbour].getXYZ())
            for ring_label in in_ring_label:
                for neighbour in self.connections[ring_label]:
                    if self.isAtomInRing(neighbour) and neighbour != atom:
                        in_rings_neighbour.append(self.atoms[neighbour].getXYZ())
        except KeyError:
            pass
        #  try to calculate dihedrals, if not possible return None
        try:
            dih1 = abs(dihedral(neighbour_not_in_ring,
                                test, in_ring[0], in_rings_neighbour[0])) % 180
            if abs(dih1 - 180) < 15 or abs(dih1 - 180) > 155:
                dih1 = 180
            dih2 = abs(dihedral(neighbour_not_in_ring,
                                test, in_ring[-1], in_rings_neighbour[-1])) % 180
            if abs(dih2 - 180) < 15 or abs(dih2 - 180) > 155:
                dih2 = 180
            average_dih = (dih1 + dih2) / 2
        except TypeError:
            return None  # to pass, not all ring atoms are axial/equitorial
            # (i.e. unsaturated uronic acid)
        #  if dihedral is close to 180, return True
        return average_dih <= 135

    def XYZ(self, atomID):
        return self.atoms[atomID].getXYZ()

    def isConnected(self, atomID1, atomID2):
        return self.graph.has_edge(atomID1, atomID2)

    def isC1(self, atomID):
        """Returns True if atom is a C1"""
        #  CALLS isC5 and isO5
        is_carbon5 = self.isC5(atomID)
        inRing = self.isAtomInRing(atomID)
        #  if atom is a carbon and in a ring,
        #  check if it is connected to an O5
        if not is_carbon5 and inRing:
            for neighbours in self.connections[atomID]:
                if self.isO5(neighbours):
                    return True
        return False

    def isC2(self, atomID):
        '''Returns True if atom is a C2'''
        #  CALLS isC1
        if self.atoms[atomID].isAtom("C") \
                and self.isAtomInRing(atomID):
            for connections in self.connections[atomID]:
                if self.isC1(connections):
                    return True
        return False

    def isC3(self, atomID):
        # CALLS C5 and C2
        if self.atoms[atomID].isAtom("C") \
                and self.isAtomInRing(atomID) \
                and not self.isC5(atomID):
            for connections in self.connections[atomID]:
                if self.isC2(connections):
                    return True
        return False

    def isC4(self, atomID):
        """Returns True if atom is a C4"""
        #  CALLS C3
        if self.atoms[atomID].isAtom("C") and self.isAtomInRing(atomID):
            for neighbours in self.getConnections(atomID):
                if self.isC3(neighbours):
                    return True
        return False

    def isO5(self, atom_id):
        #  IF atom is an oxygen and is connected to two atoms
        if len(self.getConnections(atom_id)) == 2 \
                and self.atoms[atom_id].isAtom("O") \
                and self.isAtomInRing(atom_id):
            return True
        return False  # not an O5

    def isGylcosidic(self, atomID):
        if len(self.connections[atomID]) == 2 \
                and self.atoms[atomID].isAtom("O") \
                and not self.isAtomInRing(atomID):
            if self.isAtomInRing(self.connections[atomID][0]) \
                    and self.isAtomInRing(self.connections[atomID][1]):
                return True
        return False  # not a glycosidic bond

    def isC5(self, atomID):
        if self.atoms[atomID].isAtom("C") and self.isAtomInRing(atomID):
            T = True
            for connections in self.connections[atomID]:
                if self.atoms[connections].isAtom("O") \
                        and not self.isAtomInRing(connections):
                    T = False
            if T:
                if len([_ for _ in self.atoms[atomID].connections
                        if self.atoms[_].isAtom("H") >= 3]):
                    return False
                else:
                    return True
        else:
            return False

    def isAtomInRing(self, atomID):
        for ring in self.rings:
            for atom in ring:
                if atom == atomID:
                    return True
        return False

    def get_name(self):
        return self.carbohydrate_name

    def set_name(self):
        name = ""
        start = ""
        end = ""
        for node in self.SugarGraph.nodes:
            if len(list(self.SugarGraph.successors(node))) == 0:
                end = node
            elif len(list(self.SugarGraph.predecessors(node))) == 0:
                start = node

        self.ordered_ring_hashes = shortestpath = \
            nx.shortest_path(self.SugarGraph, start, end)

        for c in shortestpath:
            name += self.Rings[c].get_ring_name()
            try:
                self.ordered_linkages.append(
                    self.Linkages[self.get_ring(self.Rings[c].c1.id)])
            except KeyError:
                pass  # Occurs when final ring, which is not a linkage key, is found
            try:
                name += "-{}-".format(self.Linkages[c].get_linkage_type())
            except KeyError:
                pass

        self.carbohydrate_name = catchNameExceptions(name)

    def getAnomericName(self, ring):
        for ring_from_list, name, a in self.ring_names:
            if ring == ring_from_list:
                return [name, a]

    def getSNFGname(self):
        for linkage in self.psi_angles:
            self.SugarGraph.add_node(self.get_ring(linkage[3][0]))
            self.SugarGraph.add_node(self.get_ring(linkage[3][1]))
            self.SugarGraph.add_edge(self.get_ring(linkage[3][0]), self.get_ring(linkage[3][1]))
        name = ""
        start = ""
        end = ""
        for node in self.SugarGraph.nodes:
            if len(list(self.SugarGraph.successors(node))) == 0:
                end = node
            elif len(list(self.SugarGraph.predecessors(node))) == 0:
                start = node
        shortestpath = nx.shortest_path(self.SugarGraph, start, end)
        for c in shortestpath:
            name += getSNFGname(self.get_ring_name(c))
        return name

    def linkage_hash(self, ring_hash_1, ring_hash_2):
        return ring_hash_1 * ring_hash_2




if __name__ == "__main__":
    from sys import argv
    args = argv[1:]
    if len(args) != 1:
        print("Usage: python3 carbohydrate.py <pdb file>")
        exit(1)
    c = Carbohydrate(args[0])
    print(c.ordered_linkages[0].__dict__)
    for _ in c.Rings.values():
        print(_.print_functional_groups())