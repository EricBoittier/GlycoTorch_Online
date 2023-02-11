out_path = "/home/eric/Projects/GlycoTorch/data/ligands/pdb_fixed/"
class PDB_Out(object):
    def __init__(self, PDB):
        self.lines = PDB.lines
        self.outlines = []
        for line in self.lines:
            if line.__contains__("HETATM") or line.__contains__("ATOM"):
                id = int(line.split()[1])
                atom = PDB.atoms[id]

                s = atom.atom_type[0:4]
                print(len(s))
                c = 0
                while c < 4 - len(s):
                    s += " "
                    c += 1

                line = line[:13] + s + line[16:]
                self.outlines.append(line)
            else:
                self.outlines.append(line)

    def make_output(self):
        with open(out_path+"2NWG_LIGAND_2.pdb", "w") as f:
            f.writelines(self.outlines)
