import os
from Carbohydrate import *
# from FindCarbohydrateName import *

non_GAG_residues = ["NI", "CH2", "SO4", "TDG", 'CAC', 'ASG', 'CA', 'EPE', 'IPA', 'NA', 'PO4', 'GOL', 'A3P', 'CIT',
                    'MPD', 'ACY'
                    , '0G6', 'ZN', 'FMT', 'PEG', 'MG', 'E64', 'NO3', 'ACT', 'JHM', 'TLA', 'DTT', 'PCA', 'CO3', 'MN',
                    'PA5', 'FAD', 'THJ', 'SIA', 'PT', 'K', 'AMP',  'TL',  'MRD', 'DMJ', 'CS', 'ACP', 'CO', 'EDO',
                    'LDA', 'ACE', 'RET', 'PLM', 'C8E', 'BME', 'HG',  'FOR',  'BEK', 'ADE', 'NHE',  'BEZ', 'HEM',
                    'OS', 'NH2',  'CD', 'HEZ', 'ASO',  'VO4',  'FSM', 'PG4',  'BEN',  '6PG', 'NPO', '2PO', 'ACH',
                    'PEP', 'NO2', 'UNX',  '293', 'IFL', 'XX6',  'XX7', 'AGG', 'MPT', 'PGE', 'FE', '1PG', 'GBL','MES',
                    'A46', 'AZI', 'AVE', 'AVF', 'AVD',  'BCD', 'OH', 'CYS']


class Count:
    def __init__(self, start):
        self.count = start

    def forward(self):
        self.count += 1
        return self.count

    def getCount(self):
        return self.count


class ProteinPDB(PDB):
    '''
    A class to handle a .pdb file and extract information from it
    '''

    def __init__(self, filename):
        super().__init__(filename)
        self.resolution = "-"
        self.protein_name = "-"
        self.description = "-"
        self.pH = "-"
        self.heteroatoms = []
        self.DOI = ""
        print(self.graph)
        self.ligands = list(self.graph.subgraph(c) for c in nx.connected_components(self.graph))
        self.all_ligands_list = []
        temp = []

        self.findInfo()

        self.known_residues = []

        self.residue_atoms = {}

        self.findKnownResidues()

        for ligand in self.ligands:
            if len(ligand.nodes) > 1:
                ligand = list(ligand.nodes)
                catch = True
                for atomID in ligand:
                    try:
                        if not self.atoms[atomID].isHeteroAtom() and not self.atoms[atomID].isNonGAG():
                            catch = False
                    except KeyError:
                        print("key error")
                        pass
                if catch:
                    temp.append(ligand)
        self.ligands = temp
        self.makeAllLigandsList()

    def makeAllLigandsList(self):
        for ligand in self.ligands:
            for id in ligand:
                self.all_ligands_list.append(id)

    def findKnownResidues(self):
        for line in self.heteroatoms:
            if not self.known_residues.__contains__(line):
                self.known_residues.append(line)
                self.residue_atoms[line.split()[1]] = line

    def findInfo(self):
        for line in self.lines:
            if line.startswith("COMPND   2 MOLECULE"):
                self.protein_name = line.split(" ", 5)[5].replace(
                    ";", "").strip("\n").capitalize()
            if line.startswith("REMARK 200  PH  "):
                self.pH = line.split(":")[-1]
                if self.pH.__contains__("NULL"):
                    self.pH = "-"
            if line.startswith("JRNL        DOI"):
                self.DOI = line.split()[2].strip("\n")
            if line.startswith("HEADER"):
                self.description = line.split(
                    " ", 1)[1][0:-40].strip(" ").capitalize()
            if line.startswith("REMARK   2 RESOLUTION."):
                self.resolution = line.split()[3]
            if line.startswith("HETATM"):
                if line.split()[3] != "MSE":
                    self.heteroatoms.append(line)

    def getInfo(self):
        return [str(self.filename.strip("PDBs/")), self.resolution, self.protein_name, self.description, self.DOI, self.pH]

    def getLigand(self, int):
        return self.ligands[int]

    def getLigands(self):
        return list(self.ligands)

    def getProteinName(self):
        return self.protein_name

    def getKnownResidues(self):
        return self.known_residues

    def getHeteroAtoms(self):
        return self.heteroatoms

    def saveProteinWithoutLigands(self, save_path):
        print(self.filename)
        with open("{}/{}_WithoutLigands.pdb".format(save_path, self.filename), "w") as file:
            for line in self.lines:
                try:
                    if line.split()[1] not in self.all_ligands_list:
                        file.write(line)
                except IndexError:
                    # This is ok, it just means the line is only one word long, e.g. END, so it's valid
                    file.write(line)

    def saveLigands(self, load_path):
        count = 1
        print("ligands:")
        print(self.ligands)

        for ligand in self.ligands:

            print("saving {}/{}_LIGAND_{}.pdb".format(load_path, self.filename, count))

            with open("temp.pdb".format(load_path, self.filename, count), "w") as file:

                connections = []
                for id in ligand:
                    for line in self.lines:
                        try:
                            if line.split()[1] == id:
                                if line.split()[0] == "HETATM":
                                    file.write(line)
                                else:
                                    connections.append(line)
                        except IndexError:
                            # This is ok, it just means the line is only one word long, e.g. END, so it's valid
                            pass

                for connection in connections:
                    file.write(connection)

            you_shall_pass = True

            try:
                l = Carbohydrate("temp.pdb")
                print("rings ", len(l.getRings()))
            except:
                you_shall_pass = False

            if you_shall_pass and len(l.getRings()) > 1:
                with open("{}/{}_LIGAND_{}.pdb".format(load_path, self.filename, count), "w") as file:

                    connections = []
                    for id in ligand:
                        for line in self.lines:
                            try:
                                if line.split()[1] == id:
                                    if line.split()[0] == "HETATM":
                                        file.write(line)
                                    else:
                                        connections.append(line)
                            except IndexError:
                                # This is ok, it just means the line is only one word long, e.g. END, so it's valid
                                pass

                    for connection in connections:
                        file.write(connection)

                count += 1


