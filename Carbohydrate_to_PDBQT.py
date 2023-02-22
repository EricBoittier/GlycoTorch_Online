import Angle
from Angle import *
from Carbohydrate import *
import os
import sys
import argparse


class TMP:
    def __init__(self):
        self.tmp = []
#        self.tmp_ = open("tmp.pdb", "w")
        self.root = None
        self.branches = {}
        self.leaves = []

    def write(self, string):
        if string is not None:
            self.tmp.append(string)
#            self.tmp_.write(string)
#            print(string)

    def close(self):
        #self.tmp_.close()
        pass

    def readlines(self):
        return self.tmp


class Carbohydrate_to_PDBQT(object):
    """docstring for Carbohydrate_to_PDBQT"""

    def __init__(self, carbohydrate):
        super(Carbohydrate_to_PDBQT, self).__init__()
        self.branch_count = 0
        self.carbohydrate = carbohydrate
        self.atoms = self.carbohydrate.atoms
        self.resnum = 1
        self.c = 1
        self.atom_id_to_c = {}
        self.functional_groups = []
        self.GO_ids = []
        self.branches = []
        self.branch_counts = []
        self.end_branches = []
        self.branches_to_write_after_root = []
        self.tmp = None
        self.branch_str = "BRANCH {} {}\n"
        self.end_branch_str = "ENDBRANCH {} {}\n"

        self.root = None
        self.branches = {}
        self.leaves = []

    def process_link(self, link_index):
        linkage = self.carbohydrate.ordered_linkages[link_index]
        self.branches_to_write_after_root = []
        # keep track of glycosidic oxygen ids
        self.GO_ids.append(linkage.GO.id)

        # if link_index != 0:
        #     l = self.carbohydrate.ordered_linkages[link_index - 1]
        #     self.tmp.write("BRANCH {} {}\n".format(l.GO.id, l.CX.id))
        #     self.branch_count += 1

        # #  write the branches
        self.tmp.write("BRANCH {} {}\n".format(
            linkage.C1.id,
            linkage.GO.id))
        self.tmp.write(self.get_pdbqt_line(linkage.GO))
        self.branch_count += 1
        self.end_branches.append("ENDBRANCH {} {}\n"
                                 .format(linkage.C1.id, linkage.GO.id))

        #

    def process_ring(self, ring, link_index):
        for atom in ring.ring_cs:
            self.tmp.write(self.get_pdbqt_line(self.atoms[atom.id]))

        #  loop through functional_groups in ring
        for ci, func_group in enumerate(ring.get_functional_groups()):
            #  write functional groups for ring
            ci = ring.ring_cs[ci]
            print("ci: ", ci, func_group)
            functional_group_atoms = []
            connected_atom = None
            # for each atom in the functional group, find the connected atoms
            for a in func_group:
                #  if not the glycosidic oxygen
                if self.atoms[a].id not in self.GO_ids:
                    #  found the connected atom
                    if self.atoms[a].id in ci.connections:
                        if connected_atom is None:
                            connected_atom = self.atoms[a]
                    if self.atoms[a] not in ring.ring:
                        functional_group_atoms.append(self.atoms[a])

            if len(functional_group_atoms) > 3:
                # if rotatable group, write to the root/branch
                # will assume some double bonds are rotatable...  TODO
                # write the branch name
                print("writing branch ", ci.id, connected_atom.id)
                self.branches_to_write_after_root.append(
                    self.branch_str.format(
                        ci.id,
                        connected_atom.id)
                )
                # write the branch atoms
                self.branches_to_write_after_root. \
                    append(connected_atom)
                #  increment branch count
                self.branch_count += 1
            else:
                #  just a non-rotatable group, write to the root/branch
                for a in functional_group_atoms:
                    self.tmp.write(self.get_pdbqt_line(a))

            # #  write the functional group to the branch
            if len(functional_group_atoms) > 3:
                for a in functional_group_atoms:
                    self.branches_to_write_after_root.append(a)

            #  write the functional group to the root
            if len(functional_group_atoms) > 3:
                self.branches_to_write_after_root.append(
                    self.end_branch_str.format(
                        ci.id,
                        connected_atom.id)
                )

        #  write end of ring1
        if link_index == 0:
            self.tmp.write("ENDROOT\n")

        #  write branches after root
        for line in self.branches_to_write_after_root:
            #  branch id
            if isinstance(line, str):
                self.tmp.write(line)
            else:  # an atom, writes in first/out first
                self.tmp.write(self.get_pdbqt_line(line))

        # #  append end branch to list if not in the first linkage
        # if link_index != 0:
        #     l = self.carbohydrate.ordered_linkages[link_index - 1]
        #     self.end_branches.append(self.end_branch_str.format(
        #         l.GO.id,
        #         l.CX.id))

        #  increment resnum
        self.resnum += 1



    def save_flex(self, path: str = None):
        """ Save a flexible carbohydrate to a pdbqt file

        :param path:
        :return:
        """
        if path:
            self.carbohydrate.filepath \
                = os.path.join(path, self.carbohydrate.filename)
            print(self.carbohydrate.filepath)

        self.tmp = TMP()  # open("tmp.pdb", "w")
        #  glycosidic oxygen ids (i.e. the main branch points)

        # write root
        self.tmp.write("ROOT\n")

        #  loop through the linkages in order
        for link_index, linkage in enumerate(
                self.carbohydrate.ordered_linkages):
            #  process ring1
            print(f"L{link_index}" "ring1")
            self.process_ring(linkage.ring1, link_index)
            self.process_link(link_index)

            #  if the last linkage
            if link_index == len(self.carbohydrate.ordered_linkages) - 1:
                # self.tmp.write("BRANCH {} {}\n".
                #                format(linkage.GO.id,
                #                       linkage.CX.id))
                # self.branch_count += 1
                # self.end_branches.append("ENDBRANCH {} {}\n".
                #                          format(linkage.GO.id,
                #                                 linkage.CX.id))
                #  process ring2
                print(f"L{link_index} ring2")
                self.process_ring(linkage.ring2, link_index)

        self.end_branches.reverse()
        for line in self.end_branches:
            self.tmp.write(line)
        self.tmp.write("TORSDOF {}".format(self.branch_count))
        self.tmp.close()

        #  open the temporary file and write the pdbqt file
        # tmp = open("tmp.pdb", "r")
        print(self.carbohydrate.filepath + ".pdbqt")
        test = open(self.carbohydrate.filepath + ".pdbqt", "w")
        #  some files added endroot after it had already been added, checking for consistency
        end_root_coming = True
        for line in self.tmp.readlines():
            if line.__contains__("BRANCH"):
                branch_list = line.split()
                try:
                    test.write(branch_list[0]
                               + " {} {}\n".format(
                        self.atom_id_to_c[int(branch_list[1])],
                        self.atom_id_to_c[int(branch_list[2])]
                    ))
                except:
                    pass
            else:
                if line.__contains__("ENDROOT"):
                    #  only write endroot once
                    if end_root_coming:
                        test.write(line)
                        end_root_coming = False
                else:
                    test.write(line)

        test.close()

    def get_pdbqt_line(self, atom):
        id = self.pad_before(self.c, 7)
        print(atom.id)
        if int(atom.id) not in self.atom_id_to_c.keys():
            self.atom_id_to_c[int(atom.id)] = int(self.c)
            self.c += 1
        else:
            return None

        name = self.pad_before(atom.atom_type, 4) + " "
        ligand_type = self.pad_before(atom.ligand_type, 4)
        chain = self.pad_before(atom.chain, 2)
        res_id = self.pad_before(self.resnum, 4)
        x = "{:7.3f}".format(atom.x)
        x = self.pad_before(x, 12)
        y = "{:7.3f}".format(atom.y)
        y = self.pad_after(y, 8)
        z = "{:7.3f}".format(atom.z)

        add_H = False
        ADT = atom.atomname
        if atom.atomname == "O" and len(atom.connections) == 1:
            if not self.atoms[atom.connections[0]].atomname.__contains__("S") \
                    and not self.is_next_to_trigonal_planar_carbon(atom) and not self.is_carboxylate(atom):
                add_H = True
                ADT = "O"
            else:
                ADT = "OA"

        if atom.atom_type.__contains__("N"):
            if len(atom.connections) == 2:
                add_H = True
                ADT = "N"
            elif len(atom.connections) < 4:
                ADT = "NA"
            else:
                ADT = "N"

        line = "ATOM{}{}{}{}{}{} {} {}  0.00  0.00     0.000 {}\n".format(id, name, ligand_type, chain, res_id,
                                                                          x, y, z, ADT)
        # if add_H:
        #     id = self.pad_before(self.c, 7)
        #     self.c += 1
        #     # self.atom_id_to_c[int(atom.id)] = self.c
        #     x = "{:7.3f}".format(atom.x + 0.5)
        #     x = self.pad_before(x, 12)
        #     y = "{:7.3f}".format(atom.y + 0.5)
        #     y = self.pad_after(y, 8)
        #     z = "{:7.3f}".format(atom.z + 0.5)
        #     name = self.pad_before("HD", 4) + " "
        #     line += "ATOM{}{}{}{}{}{} {} {}  0.00  0.00     0.000 {}\n".format(id, name, ligand_type, chain,
        #                                                                        res_id, x, y, z, "HD")
        return line

    def is_carboxylate(self, atom):
        neighbour = self.atoms[atom.id].connections[0]
        neighbour_connections = self.atoms[neighbour].connections
        count = 0
        for id in neighbour_connections:
            if self.atoms[id].atomname.__contains__("O"):
                count += 1
        if count == 2:
            return True
        return False

    def is_carboxylate_carbon(self, atom):
        if self.is_trigonal_planar(atom):
            connections = atom.connections
            count = 0
            for id in connections:
                if self.atoms[id].atomname.__contains__("O"):
                    count += 1
            if count == 2:
                return True
            return False

    def is_sulfate(self, atom):
        if atom.atom_type == "S":
            neighbour_connections = atom.connections
            count = 0
            for id in neighbour_connections:
                if self.atoms[id].atomname.__contains__("O"):
                    count += 1
            if count == 4:
                return True
        return False

    def is_Nsulfate(self, atom):
        if atom.atom_type == "S":
            neighbour_connections = atom.connections
            n_count = 0
            count = 0
            for id in neighbour_connections:
                if self.atoms[id].atomname.__contains__("O"):
                    count += 1
                if self.atoms[id].atomname.__contains__("N"):
                    n_count += 1
            if count == 3 and n_count == 1:
                return True
        return False

    def is_trigonal_planar(self, atom):
        connections_xyz = []
        for connection in atom.connections:
            a = self.atoms[connection]
            xyz = a.getXYZ()
            connections_xyz.append(xyz)

        if len(connections_xyz) != 3:
            return False

        atom_xyz = atom.getXYZ()
        angle_1 = Angle.angle_between(connections_xyz[0], atom_xyz, connections_xyz[1])
        angle_2 = Angle.angle_between(connections_xyz[1], atom_xyz, connections_xyz[2])
        angle_3 = Angle.angle_between(connections_xyz[2], atom_xyz, connections_xyz[0])
        sum_of_angles = abs(angle_1) + abs(angle_2) + abs(angle_3)

        return abs(360 - sum_of_angles) < 5

    def is_next_to_trigonal_planar_carbon(self, atom):
        neighbour = atom.connections[0]
        neighbour = self.atoms[neighbour]
        neighbour_xyz = neighbour.getXYZ()
        connections_xyz = []
        for connection in neighbour.connections:
            a = self.atoms[connection]
            xyz = a.getXYZ()
            connections_xyz.append(xyz)

        if len(connections_xyz) < 3:
            return False

        angle_1 = Angle.angle_between(connections_xyz[0], neighbour_xyz, connections_xyz[1])
        angle_2 = Angle.angle_between(connections_xyz[1], neighbour_xyz, connections_xyz[2])
        angle_3 = Angle.angle_between(connections_xyz[2], neighbour_xyz, connections_xyz[0])
        sum_of_angles = abs(angle_1) + abs(angle_2) + abs(angle_3)

        return abs(360 - sum_of_angles) < 5

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


# c = Carbohydrate("/home/eric/Documents/github/GlycoTorch-Vina/Tutorial/2axm_ligand.pdb")
# c = Carbohydrate("/home/eric/Documents/github/GlycoTorch-Vina/Tutorial/glycam_test1.pdb")
# c = Carbohydrate("/home/EricBoittier/GlycoTorch-Vina/Tutorial/2axm_ligand.pdb")
# c = Carbohydrate("/Users/ericboittier/Downloads/Unsulphated_HS_tetramer_glycam.pdb")
# pdbqt = Carbohydrate_to_PDBQT(c)
# pdbqt.save_flex()

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description='Convert a GAG ligand from pdb to pdbqt format')
#     parser.add_argument("--file", type=str, default=None, help="full path to file")
#     args = parser.parse_args()
#     c = Carbohydrate(args.file)
#     pdbqt = Carbohydrate_to_PDBQT(c)
#     pdbqt.save_flex()
