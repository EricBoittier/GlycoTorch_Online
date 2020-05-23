from Angle import *
from Carbohydrate import *
import os
import sys
import argparse



class Carbohydrate_to_PDBQT(object):
    """docstring for Carbohydrate_to_PDBQT"""

    def __init__(self, Carbohydrate):
        super(Carbohydrate_to_PDBQT, self).__init__()
        self.Carbohydrate = Carbohydrate
        self.atoms = self.Carbohydrate.atoms
        self.resnum = 1
        self.c = 1
        self.atom_id_to_c = {}

    def save_rigid(self, path=None):
        if path:
            self.Carbohydrate.filepath = path + "/" + self.Carbohydrate.filename
            print(self.Carbohydrate.filepath)

        tmp = open("./tmp.pdb", "w")

        GO_ids = []

        end_branches = []

        tmp.write("ROOT\n")

        for index in range(len(self.Carbohydrate.ordered_linkages)):

            linkage = self.Carbohydrate.ordered_linkages[index]

            GO_ids.append(linkage.GO.id)

            if index == 0:
                pass
            else:
                l = self.Carbohydrate.ordered_linkages[index - 1]

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
                l = self.Carbohydrate.ordered_linkages[index - 1]

            self.resnum += 1

            tmp.write(self.get_pdbqt_line(linkage.GO))

            if index == len(self.Carbohydrate.ordered_linkages) - 1:

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

        tmp.write("TORSDOF {}".format(2 * len(self.Carbohydrate.ordered_linkages)))
        tmp.close()

        tmp = open("./tmp.pdb", "r")
        test = open(self.Carbohydrate.filepath + ".pdbqt", "w")

        for line in tmp.readlines():
            if line.__contains__("BRANCH"):
                branch_list = line.split()
                test.write(branch_list[0] + " {} {}\n".format(self.atom_id_to_c[int(branch_list[1])],
                                                              self.atom_id_to_c[int(branch_list[2])]))
            else:
                test.write(line)

    def save_flex(self, path=None):

        if path:
            self.Carbohydrate.filepath = os.path.join(path, self.Carbohydrate.filename)
            print(self.Carbohydrate.filepath)

        tmp = open("./tmp.pdb", "w")
        GO_ids = []
        end_branches = []

        branch_count = 0

        for index in range(len(self.Carbohydrate.ordered_linkages)):

            linkage = self.Carbohydrate.ordered_linkages[index]

            GO_ids.append(linkage.GO.id)

            if index == 0:
                tmp.write("ROOT\n")
            else:
                l = self.Carbohydrate.ordered_linkages[index - 1]
                tmp.write("BRANCH {} {}\n".format(l.GO.id, l.CX.id))
                branch_count += 1

            #  write ring atoms for ring1
            for atom in linkage.ring1.ring:
                tmp.write(self.get_pdbqt_line(self.atoms[atom]))

            branches_to_write_after_root = []

            #  write functional groups for ring1

            # functional group 1
            functional_group_atoms = []
            connected_atom = False
            for a in linkage.ring1.c1_functional_group:
                if self.atoms[a].id not in GO_ids:

                    if self.atoms[a].id in linkage.ring1.c1.connections:
                        connected_atom = self.atoms[a]
                    else:
                        functional_group_atoms.append(self.atoms[a])

            if connected_atom:
                if len(linkage.ring1.c1_functional_group) > 2:
                    branches_to_write_after_root.append("BRANCH {} {}\n".format(linkage.ring1.c1.id, connected_atom.id))
                    branches_to_write_after_root.append(connected_atom)
                    branch_count += 1
                else:
                    for a in linkage.ring1.c1_functional_group:
                        tmp.write(self.get_pdbqt_line(self.atoms[a]))

            if len(linkage.ring1.c1_functional_group) > 2:
                for a in functional_group_atoms:
                    branches_to_write_after_root.append(a)

            if connected_atom:
                if len(linkage.ring1.c1_functional_group) > 2:
                    branches_to_write_after_root.append(
                        "ENDBRANCH {} {}\n".format(linkage.ring1.c1.id, connected_atom.id))

            # functional group 2
            functional_group_atoms = []
            connected_atom = False
            for a in linkage.ring1.c2_functional_group:
                if self.atoms[a].id not in GO_ids:

                    if self.atoms[a].id in linkage.ring1.c2.connections:
                        connected_atom = self.atoms[a]
                    else:
                        functional_group_atoms.append(self.atoms[a])

            if connected_atom:
                if len(linkage.ring1.c2_functional_group) > 2:
                    branches_to_write_after_root.append("BRANCH {} {}\n".format(linkage.ring1.c2.id, connected_atom.id))
                    branches_to_write_after_root.append(connected_atom)
                    branch_count += 1
                else:
                    for a in linkage.ring1.c2_functional_group:
                        tmp.write(self.get_pdbqt_line(self.atoms[a]))

            if len(linkage.ring1.c2_functional_group) > 2:
                for a in functional_group_atoms:
                    branches_to_write_after_root.append(a)

            if connected_atom:
                if len(linkage.ring1.c2_functional_group) > 2:
                    branches_to_write_after_root.append(
                        "ENDBRANCH {} {}\n".format(linkage.ring1.c2.id, connected_atom.id))

            # functional group 3
            functional_group_atoms = []
            connected_atom = False
            for a in linkage.ring1.c3_functional_group:
                if self.atoms[a].id not in GO_ids:

                    if self.atoms[a].id in linkage.ring1.c3.connections:
                        connected_atom = self.atoms[a]
                    else:
                        functional_group_atoms.append(self.atoms[a])

            if connected_atom:
                if len(linkage.ring1.c3_functional_group) > 2:
                    branches_to_write_after_root.append("BRANCH {} {}\n".format(linkage.ring1.c3.id, connected_atom.id))
                    branches_to_write_after_root.append(connected_atom)
                    branch_count += 1
                else:
                    for a in linkage.ring1.c3_functional_group:
                        tmp.write(self.get_pdbqt_line(self.atoms[a]))

            if len(linkage.ring1.c3_functional_group) > 2:
                for a in functional_group_atoms:
                    branches_to_write_after_root.append(a)

            if connected_atom:
                if len(linkage.ring1.c3_functional_group) > 2:
                    branches_to_write_after_root.append(
                        "ENDBRANCH {} {}\n".format(linkage.ring1.c3.id, connected_atom.id))

            # functional group 4
            functional_group_atoms = []
            connected_atom = False
            for a in linkage.ring1.c4_functional_group:
                if self.atoms[a].id not in GO_ids:

                    if self.atoms[a].id in linkage.ring1.c4.connections:
                        connected_atom = self.atoms[a]
                    else:
                        functional_group_atoms.append(self.atoms[a])

            if connected_atom:
                if len(linkage.ring1.c4_functional_group) > 2:
                    branches_to_write_after_root.append("BRANCH {} {}\n".format(linkage.ring1.c4.id, connected_atom.id))
                    branches_to_write_after_root.append(connected_atom)
                    branch_count += 1
                else:
                    for a in linkage.ring1.c4_functional_group:
                        tmp.write(self.get_pdbqt_line(self.atoms[a]))

            if len(linkage.ring1.c4_functional_group) > 2:
                for a in functional_group_atoms:
                    branches_to_write_after_root.append(a)

            if connected_atom:
                if len(linkage.ring1.c4_functional_group) > 2:
                    branches_to_write_after_root.append(
                        "ENDBRANCH {} {}\n".format(linkage.ring1.c4.id, connected_atom.id))

            # functional group 5
            functional_group_atoms = []
            connected_atom = False
            for a in linkage.ring1.c5_functional_group:
                if self.atoms[a].id not in GO_ids:

                    if self.atoms[a].id in linkage.ring1.c5.connections:
                        connected_atom = self.atoms[a]
                    else:
                        functional_group_atoms.append(self.atoms[a])

            if connected_atom:
                if len(linkage.ring1.c5_functional_group) > 2:
                    branches_to_write_after_root.append("BRANCH {} {}\n".format(linkage.ring1.c5.id, connected_atom.id))
                    branches_to_write_after_root.append(connected_atom)
                    branch_count += 1
                else:
                    for a in linkage.ring1.c5_functional_group:
                        tmp.write(self.get_pdbqt_line(self.atoms[a]))

            if len(linkage.ring1.c5_functional_group) > 2:
                for a in functional_group_atoms:
                    branches_to_write_after_root.append(a)

            if connected_atom:
                if len(linkage.ring1.c5_functional_group) > 2:
                    branches_to_write_after_root.append(
                        "ENDBRANCH {} {}\n".format(linkage.ring1.c5.id, connected_atom.id))

            #  end writing functional groups for ring1

            if index == 0:
                tmp.write("ENDROOT\n")

            for line in branches_to_write_after_root:
                if isinstance(line, str):
                    tmp.write(line)
                else:
                    tmp.write(self.get_pdbqt_line(line))

            if index > 0:
                l = self.Carbohydrate.ordered_linkages[index - 1]
                end_branches.append("ENDBRANCH {} {}\n".format(l.GO.id, l.CX.id))

            self.resnum += 1

            tmp.write("BRANCH {} {}\n".format(linkage.C1.id, linkage.GO.id))
            tmp.write(self.get_pdbqt_line(linkage.GO))
            branch_count += 1
            end_branches.append("ENDBRANCH {} {}\n".format(linkage.C1.id, linkage.GO.id))



            if index == len(self.Carbohydrate.ordered_linkages) - 1:
                tmp.write("BRANCH {} {}\n".format(linkage.GO.id, linkage.CX.id))
                branch_count += 1
                end_branches.append("ENDBRANCH {} {}\n".format(linkage.GO.id, linkage.CX.id))

                for atom in linkage.ring2.ring:
                    tmp.write(self.get_pdbqt_line(self.atoms[atom]))

                # functional group 1
                functional_group_atoms = []
                connected_atom = False
                for a in linkage.ring2.c1_functional_group:
                    if self.atoms[a].id not in GO_ids:

                        if self.atoms[a].id in linkage.ring2.c1.connections:
                            connected_atom = self.atoms[a]
                        else:
                            functional_group_atoms.append(self.atoms[a])

                if connected_atom:
                    if len(linkage.ring2.c1_functional_group) > 2:
                        tmp.write("BRANCH {} {}\n".format(linkage.ring2.c1.id, connected_atom.id))
                        tmp.write(self.get_pdbqt_line(connected_atom))
                        branch_count += 1
                    else:
                        for a in functional_group_atoms:
                            tmp.write(self.get_pdbqt_line(a))

                for a in functional_group_atoms:
                    tmp.write(self.get_pdbqt_line(a))

                if connected_atom:
                    if len(linkage.ring2.c1_functional_group) > 2:
                        tmp.write(
                            "ENDBRANCH {} {}\n".format(linkage.ring2.c1.id, connected_atom.id))

                # functional group 2
                functional_group_atoms = []
                connected_atom = False
                for a in linkage.ring2.c2_functional_group:
                    if self.atoms[a].id not in GO_ids:

                        if self.atoms[a].id in linkage.ring2.c2.connections:
                            connected_atom = self.atoms[a]
                        else:
                            functional_group_atoms.append(self.atoms[a])

                if connected_atom:
                    if len(linkage.ring2.c2_functional_group) > 2:
                        tmp.write("BRANCH {} {}\n".format(linkage.ring2.c2.id, connected_atom.id))
                        tmp.write(self.get_pdbqt_line(connected_atom))
                        branch_count += 1
                    else:
                        for a in functional_group_atoms:
                            tmp.write(self.get_pdbqt_line(a))

                for a in functional_group_atoms:
                    tmp.write(self.get_pdbqt_line(a))

                if connected_atom:
                    if len(linkage.ring2.c2_functional_group) > 2:
                        tmp.write(
                            "ENDBRANCH {} {}\n".format(linkage.ring2.c2.id, connected_atom.id))

                # functional group 3
                functional_group_atoms = []
                connected_atom = False
                for a in linkage.ring2.c3_functional_group:
                    if self.atoms[a].id not in GO_ids:

                        if self.atoms[a].id in linkage.ring2.c3.connections:
                            connected_atom = self.atoms[a]
                        else:
                            functional_group_atoms.append(self.atoms[a])

                if connected_atom:
                    if len(linkage.ring2.c3_functional_group) > 2:
                        tmp.write("BRANCH {} {}\n".format(linkage.ring2.c3.id, connected_atom.id))
                        tmp.write(self.get_pdbqt_line(connected_atom))
                        branch_count += 1
                    else:
                        for a in functional_group_atoms:
                            tmp.write(self.get_pdbqt_line(a))

                for a in functional_group_atoms:
                    tmp.write(self.get_pdbqt_line(a))

                if connected_atom:
                    if len(linkage.ring2.c3_functional_group) > 2:
                        tmp.write(
                            "ENDBRANCH {} {}\n".format(linkage.ring2.c3.id, connected_atom.id))

                # functional group 4
                functional_group_atoms = []
                connected_atom = False
                for a in linkage.ring2.c4_functional_group:
                    if self.atoms[a].id not in GO_ids:

                        if self.atoms[a].id in linkage.ring2.c4.connections:
                            connected_atom = self.atoms[a]
                        else:
                            functional_group_atoms.append(self.atoms[a])

                if connected_atom:
                    if len(linkage.ring2.c4_functional_group) > 2:
                        tmp.write("BRANCH {} {}\n".format(linkage.ring2.c4.id, connected_atom.id))
                        tmp.write(self.get_pdbqt_line(connected_atom))
                        branch_count += 1
                    else:
                        for a in functional_group_atoms:
                            tmp.write(self.get_pdbqt_line(a))

                for a in functional_group_atoms:
                    tmp.write(self.get_pdbqt_line(a))

                if connected_atom:
                    if len(linkage.ring2.c4_functional_group) > 2:
                        tmp.write(
                            "ENDBRANCH {} {}\n".format(linkage.ring2.c4.id, connected_atom.id))

                # functional group 5
                functional_group_atoms = []
                connected_atom = False
                for a in linkage.ring2.c5_functional_group:
                    if self.atoms[a].id not in GO_ids:

                        if self.atoms[a].id in linkage.ring2.c5.connections:
                            connected_atom = self.atoms[a]
                        else:
                            functional_group_atoms.append(self.atoms[a])

                if connected_atom:
                    if len(linkage.ring2.c5_functional_group) > 2:
                        tmp.write("BRANCH {} {}\n".format(linkage.ring2.c5.id, connected_atom.id))
                        tmp.write(self.get_pdbqt_line(connected_atom))
                        branch_count += 1
                    else:
                        for a in functional_group_atoms:
                            tmp.write(self.get_pdbqt_line(a))

                for a in functional_group_atoms:
                    tmp.write(self.get_pdbqt_line(a))

                if connected_atom:
                    if len(linkage.ring2.c5_functional_group) > 2:
                        tmp.write(
                            "ENDBRANCH {} {}\n".format(linkage.ring2.c5.id, connected_atom.id))


        end_branches.reverse()
        for line in end_branches:
            tmp.write(line)
        tmp.write("TORSDOF {}".format(branch_count))
        tmp.close()

        tmp = open("./tmp.pdb", "r")
        test = open(self.Carbohydrate.filepath + ".pdbqt", "w")

        for line in tmp.readlines():
            if line.__contains__("BRANCH"):
                branch_list = line.split()
                test.write(branch_list[0] + " {} {}\n".format(self.atom_id_to_c[int(branch_list[1])],
                                                              self.atom_id_to_c[int(branch_list[2])]))
            else:
                test.write(line)

    def save(self, path=None):

        if path:
            self.Carbohydrate.filepath = path + "/" + self.Carbohydrate.filename
            print(self.Carbohydrate.filepath)

        tmp = open("./tmp.pdb", "w")

        GO_ids = []

        end_branches = []

        for index in range(len(self.Carbohydrate.ordered_linkages)):

            linkage = self.Carbohydrate.ordered_linkages[index]

            GO_ids.append(linkage.GO.id)

            if index == 0:
                tmp.write("ROOT\n")
            else:
                l = self.Carbohydrate.ordered_linkages[index - 1]
                tmp.write("BRANCH {} {}\n".format(l.GO.id, l.CX.id))

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
                tmp.write("ENDROOT\n")
            else:
                l = self.Carbohydrate.ordered_linkages[index - 1]
                end_branches.append("ENDBRANCH {} {}\n".format(l.GO.id, l.CX.id))

            self.resnum += 1

            tmp.write("BRANCH {} {}\n".format(linkage.C1.id, linkage.GO.id))
            tmp.write(self.get_pdbqt_line(linkage.GO))
            end_branches.append("ENDBRANCH {} {}\n".format(linkage.C1.id, linkage.GO.id))

            if index == len(self.Carbohydrate.ordered_linkages) - 1:
                tmp.write("BRANCH {} {}\n".format(linkage.GO.id, linkage.CX.id))
                end_branches.append("ENDBRANCH {} {}\n".format(linkage.GO.id, linkage.CX.id))

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

        end_branches.reverse()
        for line in end_branches:
            tmp.write(line)
        tmp.write("TORSDOF {}".format(2 * len(self.Carbohydrate.ordered_linkages)))
        tmp.close()

        tmp = open("./tmp.pdb", "r")
        test = open(self.Carbohydrate.filepath + ".pdbqt", "w")

        for line in tmp.readlines():
            if line.__contains__("BRANCH"):
                branch_list = line.split()
                test.write(branch_list[0] + " {} {}\n".format(self.atom_id_to_c[int(branch_list[1])],
                                                              self.atom_id_to_c[int(branch_list[2])]))
            else:
                test.write(line)

    def get_pdbqt_line(self, atom):
        id = self.pad_before(self.c, 7)
        self.atom_id_to_c[int(atom.id)] = self.c
        self.c += 1
        name = self.pad_before(atom.atom_type, 4) + " "
        ligand_type = self.pad_before(atom.ligand_type, 4)
        chain = self.pad_before(atom.chain, 2)
        # res_id = self.pad_before(atom.ligandID, 4)
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

        angle_1 = angle_between(connections_xyz[0], atom_xyz, connections_xyz[1])
        angle_2 = angle_between(connections_xyz[1], atom_xyz, connections_xyz[2])
        angle_3 = angle_between(connections_xyz[2], atom_xyz, connections_xyz[0])
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

        angle_1 = angle_between(connections_xyz[0], neighbour_xyz, connections_xyz[1])
        angle_2 = angle_between(connections_xyz[1], neighbour_xyz, connections_xyz[2])
        angle_3 = angle_between(connections_xyz[2], neighbour_xyz, connections_xyz[0])
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert a GAG ligand from pdb to pdbqt format')
    parser.add_argument("--file", type=str, default=None, help="full path to file")
    args = parser.parse_args()
    c = Carbohydrate(args.file)
    pdbqt = Carbohydrate_to_PDBQT(c)
    pdbqt.save_flex()