import Angle
from Angle import *
from Carbohydrate import *
import os
import sys
import argparse



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

    def process_ring(self, ring, link_index):
        linkage = self.carbohydrate.ordered_linkages[link_index]
        # keep track of glycosidic oxygen ids
        self.GO_ids.append(linkage.GO.id)

        if link_index != 0:
            l = self.carbohydrate.ordered_linkages[link_index - 1]
            self.tmp.write("BRANCH {} {}\n".format(l.GO.id, l.CX.id))
            self.branch_count += 1

        #  write ring atoms for ring1
        for atom in ring.ring_cs:
            self.tmp.write(self.get_pdbqt_line(self.atoms[atom.id]))

        #  loop through functional_groups in ring1
        for ci, func_group in enumerate(ring.get_functional_groups()):
            #  write functional groups for ring1
            ci = ring.ring_cs[ci]
            print("func group:", func_group)
            functional_group_atoms = []
            connected_atom = False
            for a in func_group:
                if self.atoms[a].id not in self.GO_ids:
                    if self.atoms[a].id in ci.connections:
                        connected_atom = self.atoms[a]
                    else:
                        functional_group_atoms.append(self.atoms[a])
            print("fga:", functional_group_atoms)
            if connected_atom is not False:
                if len(functional_group_atoms) > 1:
                    self.branches_to_write_after_root.append(
                        self.branch_str.format(
                            ci.id,
                            connected_atom.id)
                    )
                    self.branches_to_write_after_root.append(connected_atom)
                    self.branch_count += 1
                else:
                    for a in functional_group_atoms:
                        self.tmp.write(self.get_pdbqt_line(a))

            if len(functional_group_atoms) > 2:
                for a in functional_group_atoms:
                    self.branches_to_write_after_root.append(a)

            if connected_atom is not False:
                if len(func_group) > 2:
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
            if isinstance(line, str):
                self.tmp.write(line)
            else:
                self.tmp.write(self.get_pdbqt_line(line))

        #  append end branch to list if not in the first linkage
        if link_index != 0:
            l = self.carbohydrate.ordered_linkages[link_index - 1]
            self.end_branches.append(self.end_branch_str.format(
                l.GO.id,
                l.CX.id))

        #  increment resnum
        self.resnum += 1
        #  write the branches
        self.tmp.write("BRANCH {} {}\n".format(linkage.C1.id, linkage.GO.id))
        self.tmp.write(self.get_pdbqt_line(linkage.GO))
        self.branch_count += 1
        self.end_branches.append("ENDBRANCH {} {}\n".format(linkage.C1.id, linkage.GO.id))


    def save_rigid(self, path=None):
        if path:
            self.carbohydrate.filepath = path + "/" + self.carbohydrate.filename
            print(self.carbohydrate.filepath)

        tmp = open("tmp.pdb", "w")
        GO_ids = []
        end_branches = []

        #  write header
        tmp.write("ROOT\n")
        #  loop through linkages
        for index in range(len(self.carbohydrate.ordered_linkages)):
            linkage = self.carbohydrate.ordered_linkages[index]
            #  keep track of glycosidic oxygen ids
            GO_ids.append(linkage.GO.id)
            if index != 0:
                link = self.carbohydrate.ordered_linkages[index - 1]

            #  write ring atoms for ring1
            for atom in linkage.ring1.ring:
                tmp.write(self.get_pdbqt_line(self.atoms[atom]))

            #  write functional groups for ring1
            for a in linkage.ring1.c1_functional_group:
                if self.atoms[a].id not in GO_ids:
                    tmp.write(self.get_pdbqt_line(self.atoms[a]))
            for a in linkage.ring1.c2_functional_group:
                if self.atoms[a].id not in GO_ids:
                    tmp.write(self.get_pdbqt_line(self.atoms[a]))
            for a in linkage.ring1.c3_functional_group:
                if self.atoms[a].id not in GO_ids:
                    tmp.write(self.get_pdbqt_line(self.atoms[a]))
            for a in linkage.ring1.c4_functional_group:
                if self.atoms[a].id not in GO_ids:
                    tmp.write(self.get_pdbqt_line(self.atoms[a]))
            for a in linkage.ring1.c5_functional_group:
                if self.atoms[a].id not in GO_ids:
                    tmp.write(self.get_pdbqt_line(self.atoms[a]))
            #  end writing functional groups for ring1

            if index == 0:
                pass
            else:
                link = self.carbohydrate.ordered_linkages[index - 1]

            self.resnum += 1

            tmp.write(self.get_pdbqt_line(linkage.GO))

            if index == len(self.carbohydrate.ordered_linkages) - 1:
                for atom in linkage.ring2.ring:
                    tmp.write(self.get_pdbqt_line(self.atoms[atom]))

                for a in linkage.ring2.c1_functional_group:
                    if self.atoms[a].id not in GO_ids:
                        tmp.write(self.get_pdbqt_line(self.atoms[a]))

                for a in linkage.ring2.c2_functional_group:
                    if self.atoms[a].id not in GO_ids:
                        tmp.write(self.get_pdbqt_line(self.atoms[a]))

                for a in linkage.ring2.c3_functional_group:
                    if self.atoms[a].id not in GO_ids:
                        tmp.write(self.get_pdbqt_line(self.atoms[a]))

                for a in linkage.ring2.c4_functional_group:
                    if self.atoms[a].id not in GO_ids:
                        tmp.write(self.get_pdbqt_line(self.atoms[a]))

                for a in linkage.ring2.c5_functional_group:
                    if self.atoms[a].id not in GO_ids:
                        tmp.write(self.get_pdbqt_line(self.atoms[a]))

        tmp.write("ENDROOT\n")

        tmp.write("TORSDOF {}".format(2 * len(self.carbohydrate.ordered_linkages)))
        tmp.close()

        tmp = open("tmp.pdb", "r")
        test = open(self.carbohydrate.filepath + ".pdbqt", "w")

        for line in tmp.readlines():
            if line.__contains__("BRANCH"):
                branch_list = line.split()
                test.write(branch_list[0] + " {} {}\n".format(self.atom_id_to_c[int(branch_list[1])],
                                                              self.atom_id_to_c[int(branch_list[2])]))
            else:
                test.write(line)

    def save_flex(self, path: str = None):
        """ Save a flexible carbohydrate to a pdbqt file

        :param path:
        :return:
        """
        if path:
            self.carbohydrate.filepath \
                = os.path.join(path, self.carbohydrate.filename)
            print(self.carbohydrate.filepath)

        self.tmp = open("tmp.pdb", "w")
        #  glycosidic oxygen ids (i.e. the main branch points)

        # write root
        self.tmp.write("ROOT\n")

        #  loop through the linkages in order
        for link_index, linkage in enumerate(
                self.carbohydrate.ordered_linkages):
            #  process ring1
            self.process_ring(linkage.ring1, link_index)

            if link_index != 0:
                l = self.carbohydrate.ordered_linkages[link_index - 1]
                self.tmp.write(self.branch_str.format(l.GO.id, l.CX.id))
                self.branch_count += 1

            #  if not in the last linkage
            if link_index == len(self.carbohydrate.ordered_linkages) - 1:
                self.tmp.write("BRANCH {} {}\n".
                               format(linkage.GO.id,
                                      linkage.CX.id))
                self.branch_count += 1
                self.end_branches.append("ENDBRANCH {} {}\n".
                                         format(linkage.GO.id,
                                                linkage.CX.id))
                #  process ring2
                self.process_ring(linkage.ring2, link_index)

        self.end_branches.reverse()
        for line in self.end_branches:
            self.tmp.write(line)
        self.tmp.write("TORSDOF {}".format(self.branch_count))
        self.tmp.close()

        #  open the temporary file and write the pdbqt file
        tmp = open("tmp.pdb", "r")
        print(self.carbohydrate.filepath + ".pdbqt")
        test = open(self.carbohydrate.filepath + ".pdbqt", "w")
        for line in tmp.readlines():
            if line.__contains__("BRANCH"):
                branch_list = line.split()
                test.write(branch_list[0]
                           + " {} {}\n".format(
                    self.atom_id_to_c[int(branch_list[1])],
                    self.atom_id_to_c[int(branch_list[2])]
                )
                           )
            else:
                test.write(line)

    # def save(self, path=None):
    #     if path:
    #         self.carbohydrate.filepath \
    #             = path + "/" + self.carbohydrate.filename
    #         print(self.carbohydrate.filepath)
    #
    #     tmp = open("tmp.pdb", "w")
    #     GO_ids = []
    #     end_branches = []
    #
    #     for index in range(len(self.carbohydrate.ordered_linkages)):
    #         linkage = self.carbohydrate.ordered_linkages[index]
    #         GO_ids.append(linkage.GO.id)
    #
    #         if index == 0:
    #             tmp.write("ROOT\n")
    #         else:
    #             l = self.carbohydrate.ordered_linkages[index - 1]
    #             tmp.write("BRANCH {} {}\n".format(l.GO.id, l.CX.id))
    #
    #         #  write ring atoms for ring1
    #         for atom in linkage.ring1.ring:
    #             tmp.write(self.get_pdbqt_line(self.atoms[atom]))
    #
    #         #  write functional groups for ring1
    #         for a in linkage.ring1.c1_functional_group:
    #             if self.atoms[a].id not in GO_ids:
    #                 tmp.write(self.get_pdbqt_line(self.atoms[a]))
    #
    #         for a in linkage.ring1.c2_functional_group:
    #             if self.atoms[a].id not in GO_ids:
    #                 tmp.write(self.get_pdbqt_line(self.atoms[a]))
    #
    #         for a in linkage.ring1.c3_functional_group:
    #             if self.atoms[a].id not in GO_ids:
    #                 tmp.write(self.get_pdbqt_line(self.atoms[a]))
    #
    #         for a in linkage.ring1.c4_functional_group:
    #             if self.atoms[a].id not in GO_ids:
    #                 tmp.write(self.get_pdbqt_line(self.atoms[a]))
    #
    #         for a in linkage.ring1.c5_functional_group:
    #             if self.atoms[a].id not in GO_ids:
    #                 tmp.write(self.get_pdbqt_line(self.atoms[a]))
    #         #  end writing functional groups for ring1
    #
    #         if index == 0:
    #             tmp.write("ENDROOT\n")
    #         else:
    #             l = self.carbohydrate.ordered_linkages[index - 1]
    #             end_branches.append("ENDBRANCH {} {}\n".format(l.GO.id, l.CX.id))
    #
    #         self.resnum += 1
    #
    #         tmp.write("BRANCH {} {}\n".format(linkage.C1.id, linkage.GO.id))
    #         tmp.write(self.get_pdbqt_line(linkage.GO))
    #         end_branches.append("ENDBRANCH {} {}\n".format(linkage.C1.id, linkage.GO.id))
    #
    #         if index == len(self.carbohydrate.ordered_linkages) - 1:
    #             tmp.write("BRANCH {} {}\n".format(linkage.GO.id, linkage.CX.id))
    #             end_branches.append("ENDBRANCH {} {}\n".format(linkage.GO.id, linkage.CX.id))
    #
    #             for atom in linkage.ring2.ring:
    #                 tmp.write(self.get_pdbqt_line(self.atoms[atom]))
    #
    #             for a in linkage.ring2.c1_functional_group:
    #                 if self.atoms[a].id not in GO_ids:
    #                     tmp.write(self.get_pdbqt_line(self.atoms[a]))
    #
    #             for a in linkage.ring2.c2_functional_group:
    #                 if self.atoms[a].id not in GO_ids:
    #                     tmp.write(self.get_pdbqt_line(self.atoms[a]))
    #
    #             for a in linkage.ring2.c3_functional_group:
    #                 if self.atoms[a].id not in GO_ids:
    #                     tmp.write(self.get_pdbqt_line(self.atoms[a]))
    #
    #             for a in linkage.ring2.c4_functional_group:
    #                 if self.atoms[a].id not in GO_ids:
    #                     tmp.write(self.get_pdbqt_line(self.atoms[a]))
    #
    #             for a in linkage.ring2.c5_functional_group:
    #                 if self.atoms[a].id not in GO_ids:
    #                     tmp.write(self.get_pdbqt_line(self.atoms[a]))
    #
    #     end_branches.reverse()
    #     for line in end_branches:
    #         tmp.write(line)
    #     tmp.write("TORSDOF {}".format(2 * len(self.carbohydrate.ordered_linkages)))
    #     tmp.close()
    #
    #     tmp = open("tmp.pdb", "r")
    #     print(self.carbohydrate.filepath + ".pdbqt")
    #     test = open(self.carbohydrate.filepath + ".pdbqt", "w")
    #
    #     for line in tmp.readlines():
    #         if line.__contains__("BRANCH"):
    #             branch_list = line.split()
    #             test.write(branch_list[0]
    #                        + " {} {}\n".format(
    #                 self.atom_id_to_c[int(branch_list[1])],
    #                 self.atom_id_to_c[int(branch_list[2])]))
    #         else:
    #             test.write(line)

    def get_pdbqt_line(self, atom):
        id = self.pad_before(self.c, 7)
        self.atom_id_to_c[int(atom.id)] = self.c
        self.c += 1
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

        """
        These lines have been commented out because they are not used in the standard program
        """
        # if self.is_carboxylate_carbon(atom):
        #     name = self.pad_before("@CA", 4) + " "
        # if self.is_sulfate(atom):
        #     name = self.pad_before("@S", 4) + " "
        # if self.is_Nsulfate(atom):
        #     name = self.pad_before("@NS", 4) + " "

        line = "ATOM{}{}{}{}{}{} {} {}  0.00  0.00     0.000 {}\n".format(id, name, ligand_type, chain, res_id,
                                                                          x, y, z, ADT)
        if add_H:
            self.c += 1
            id = self.pad_before(self.c, 7)
            self.atom_id_to_c[int(atom.id)] = self.c
            x = "{:7.3f}".format(atom.x + 0.5)
            x = self.pad_before(x, 12)
            y = "{:7.3f}".format(atom.y + 0.5)
            y = self.pad_after(y, 8)
            z = "{:7.3f}".format(atom.z + 0.5)
            name = self.pad_before("HD", 4) + " "
            line += "ATOM{}{}{}{}{}{} {} {}  0.00  0.00     0.000 {}\n".format(id, name, ligand_type, chain,
                                                                               res_id, x, y, z, "HD")
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

c = Carbohydrate("/Users/ericboittier/Downloads/Unsulphated_HS_tetramer_glycam.pdb")
pdbqt = Carbohydrate_to_PDBQT(c)
pdbqt.save_flex()

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description='Convert a GAG ligand from pdb to pdbqt format')
#     parser.add_argument("--file", type=str, default=None, help="full path to file")
#     args = parser.parse_args()
#     c = Carbohydrate(args.file)
#     pdbqt = Carbohydrate_to_PDBQT(c)
#     pdbqt.save_flex()