import math

def PDBQT_to_PDB(filepath):
    """"

    """
    file = open(filepath, "r")
    lines = file.readlines()

    atoms = [line for line in lines if line.__contains__("ATOM ") and not line.__contains__("H   ")]

    ordered_atoms = {}
    positions = {}

    for l in atoms:
        id = int(l.split()[1])
        ordered_atoms[id] = l
        positions[id] = (float(l.split()[6]),  float(l.split()[7]),  float(l.split()[8]))


    atom_keys = list(ordered_atoms.keys())
    atom_keys.sort()


    connections = {}
    used_connections = [] # atom1_atom2

    for atom1 in positions.keys():
        for atom2 in positions.keys():
            if atom1 != atom2 and distance(positions[atom1], positions[atom2]) < 2:
                if atom1 not in connections.keys():
                    connections[atom1] = []
                if "{}_{}".format(atom1, atom2) not in used_connections and "{}_{}".format(atom2, atom1):
                    used_connections.append("{}_{}".format(atom1, atom2))
                    connections[atom1].append(atom2)

    new_file = open(filepath+".pdb", "w")
    for key in atom_keys:
        new_file.write(ordered_atoms[key])
    for key in connections.keys():
        if len(connections[key]) > 0:
            new_file.write("CONECT   {}   ".format(key))
            for connection in connections[key]:
                new_file.write(str(connection)+"   ")
            new_file.write("\n")


def distance(point_1, point_2):
    return math.sqrt((point_1[0] - point_2[0])**2 + (point_1[1] - point_2[1])**2 + (point_1[2] - point_2[2])**2)



