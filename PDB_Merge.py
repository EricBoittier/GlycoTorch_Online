import math
def calculate_distance(a, b):
    """
    :param a: vector a
    :param b: vector b
    :return:
    """
    dif1= a[0] - b[0]
    dif2=a[1] - b[1]
    dif3=a[2] - b[2]

    return math.sqrt(dif1*dif1+dif2*dif2+dif3*dif3)


class PDB_Merge(object):
    """

    """
    def __init__(self):
        self.molecules = []
        self.output = []
        self.counter = Counter()
        self.ligand = ""
        self.ligand_center = []
        self.protein = []

    def add_molecule(self, PDB):
        self.molecules.append(PDB)

    def merge(self, filename):
        for pdb in self.molecules:
            file = open(pdb, "r")
            for line in file.readlines():
                if line.__contains__("ATOM") or line.__contains__("HETATM"):
                    self.output.append(line[:7]+"{}".format(self.counter.get_count()) + line[12:])

        new_file = open(filename, "w")
        print(filename)
        for line in self.output:
            new_file.write(line)

    def add_ligand(self, PDB):

        self.ligand = PDB

        f = open(PDB, "r")
        x = []
        y = []
        z = []
        count = 0
        for l in f.readlines():
            if l.__contains__("ATOM") or l.__contains__("HETATOM"):
                x.append(float(l.split()[6]))
                y.append(float(l.split()[7]))
                z.append(float(l.split()[8]))
                count += 1

        self.ligand_center = [sum(x)/count, sum(y)/count, sum(z)/count]

    def add_protein(self, PDB):

        pdb = []
        f = open(PDB+".pdb", "r")
        for l in f.readlines():
            if l.__contains__("ATOM") or l.__contains__("HETATOM"):
                x = (float(l.split()[6]))
                y = (float(l.split()[7]))
                z = (float(l.split()[8]))
                if calculate_distance([x, y, z], self.ligand_center) < 30:
                    pdb.append(l)
            else:
                pdb.append(l)

        self.protein = pdb

    def merge_protein_ligand(self, filename):
        x = self.ligand_center[0]
        y = self.ligand_center[1]
        z = self.ligand_center[2]

        #o = ["HETATM    0  H   IDS F 504     {:7.2} {:7.2} {:7.2}  1.00 43.71           H\n".format(x, y, z)]
        o = []
        # for line in self.protein:
        #     if line.__contains__("ATOM") or line.__contains__("HETATM"):
        #         o.append(line[:7] + "{}".format(self.counter.get_count()) + line[12:14]+ "  "+line[16:])
        #     else:
        #         o.append(line)

        for line in self.protein:
            if line.__contains__("ATOM") or line.__contains__("HETATM"):
                o.append("{} {} {} {}\n".format(line[13:14], line[30:37], line[38:45], line[46:53]))

        new_file = open(filename, "w")
        print(filename)
        for line in o:
            new_file.write(line)

class Counter(object):
    def __init__(self):
        self.count = 1

    def get_count(self):
        return_count = str(self.count)
        self.count += 1
        while len(return_count) < 5:
            return_count = " "+return_count

        return return_count

