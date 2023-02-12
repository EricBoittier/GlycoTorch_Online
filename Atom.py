non_GAG_residues = ["NI", "CH2", "SO4", "TDG", 'CAC', 'ASG', 'CA', 'EPE', 'IPA', 'NA', 'PO4', 'GOL', 'A3P', 'CIT', 'MPD', 'ACY'
                    , '0G6', 'ZN', 'FMT', 'PEG', 'MG', 'E64', 'NO3', 'ACT', 'JHM', 'TLA', 'DTT', 'PCA', 'CO3', 'MN', 'PA5', 'FAD', 'THJ', 'SIA', 'PT', 'K', 'AMP',  'TL',  'MRD', 'DMJ', 'CS', 'ACP', 'CO', 'EDO', 'LDA', 'ACE', 'RET', 'PLM', 'C8E', 'BME', 'HG',  'FOR',  'BEK', 'ADE', 'NHE',  'BEZ', 'HEM',  'OS', 'NH2',  'CD', 'HEZ', 'ASO',  'VO4',  'FSM', 'PG4',  'BEN',  '6PG', 'NPO', '2PO', 'ACH', 'PEP', 'NO2', 'UNX',  '293', 'IFL', 'XX6',  'XX7', 'AGG', 'MPT', 'PGE', 'FE', '1PG', 'GBL','MES', 'A46', 'AZI', 'AVE', 'AVF', 'AVD',  'BCD', 'OH', 'CYS']


class Atom(object):
    """

    """
    def __init__(self, line):
        super(Atom, self).__init__()
        self.line = line
        self.kind = line[0:5].strip()
        self.id = int(line[6:11].strip())
        self.atom_type = line[12:15].strip()
        self.ligand_type = line[17:20].strip()
        self.chain = line[21].strip()
        self.ligandID = line[22:26].strip()
        self.x = float(line[30:37].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:53].strip())
        self.atomname = line[76:79].strip()

        if len(self.atomname) > 0 and len(self.atomname) < 3:
            self.atom_type = self.atomname

        self.connections = []
        self.is_axial = None  # True for axial, False for equitorial, None for neither or has not been set...
        self.sugar_position = None  # Can be C1, C2, C3, C5,  O5 or GO (for glycosidic oxygen)

    def __repr__(self):
        return repr(self.line)

    def get_connections(self):
        return self.connections

    def add_connection(self, ID):
        self.connections.append(ID)

    def remove_connection(self, ID):
        self.connections.remove(ID)

    def get_is_axial(self):
        return self.is_axial

    def set_axial(self, b):
        """
        :param bool: True for axial, False for equitorial, None for neither or has not been set...
        """
        self.is_axial = b

    def set_sugar_position(self, pos):
        self.sugar_position = pos
        if not pos == "GO":
            self.atom_type = pos

    def get_sugar_position(self):
        return self.sugar_position

    def set_atom_type(self, type):
        self.atom_type = type

    def isC1(self):
        return self.sugar_position == "C1"

    def isC2(self):
        return self.sugar_position == "C2"

    def isC3(self):
        return self.sugar_position == "C3"

    def isC4(self):
        return self.sugar_position == "C4"

    def isC5(self):
        return self.sugar_position == "C5"

    def isO5(self):
        return self.sugar_position == "O5"

    def isGO(self):
        return self.sugar_position == "GO"

    def getID(self):
        return self.id

    def getLigandType(self):
        return self.ligand_type

    def getLigandID(self):
        return self.ligandID

    def getAtomType(self):
        return self.atom_type

    def getAtomename(self):
        return self.atomname

    def getXYZ(self):
        return [self.x, self.y, self.z]

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def getZ(self):
        return self.z

    def isAtom(self, atomtype):
        return self.atom_type.__contains__(atomtype)

    def getKind(self):
        return self.kind

    def isHeteroAtom(self):
        return self.kind.__contains__("HETAT")

    def isNonGAG(self):
        return self.ligand_type in non_GAG_residues
