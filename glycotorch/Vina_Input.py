
class VinaInput(object):
    """docstring for Vina_Input"""
    def __init__(self, Ligand_Object):
        super(VinaInput, self).__init__()
        self.Ligand_Object = Ligand_Object
        self.atoms = Ligand_Object.atoms

        self.xMin = self.yMin = self.zMin = 10000000000
        self.zMax = self.xMax = self.yMax = -10000000000
        self.centreX = self.centreY = self.centreZ = self.sizeX = self.sizeY = self.sizeZ = 0

        self.findADbox()

    def findADbox(self):
        for atom in self.atoms.values():
            atom = atom.getXYZ()

            if atom[0] < self.xMin:
                self.xMin = atom[0]

            if atom[0] > self.xMax:
                self.xMax = atom[0]

            if atom[1] < self.yMin:
                self.yMin = atom[1]

            if atom[1] > self.yMax:
                self.yMax = atom[1]

            if atom[2] < self.zMin:
                self.zMin = atom[2]

            if atom[2] > self.zMax:
                self.zMax = atom[2]

        self.centreX = int((self.xMin + self.xMax) / 2)
        self.centreY = int((self.yMin + self.yMax) / 2)
        self.centreZ = int((self.zMin + self.zMax) / 2)
        self.sizeX = int(abs(self.xMax - self.xMin))
        self.sizeY = int(abs(self.yMax - self.yMin))
        self.sizeZ = int(abs(self.zMax - self.zMin))

        print(self.sizeX+8)
        print(self.sizeY+8)
        print(self.sizeZ+8)

    def getADBox(self):
        """

        the bounds of the ligand plus 5 angstroms (2.5 angstroms either side of the box) to allow for
        flexibility of hydrogen bonding residues.

        :return: array size 6 of ints corresponding to centre x, y, z and corresponding sizes.
        """
        return [self.centreX, self.centreY, self.centreZ, self.sizeX + 5, self.sizeY + 5, self.sizeZ + 5]

    def makeConfig(self, path=None):
        if not path:
            path = "./glycotorch_benchmark"
        ADbox = self.getADBox()
        text = "center_x = {}\ncenter_y = {}\ncenter_z = {}\nsize_x = {}\nsize_y = {}\nsize_z" \
               " = {}\nenergy_range = 12\nexhaustiveness = {}\nnum_modes = 100".format(ADbox[0], ADbox[1], ADbox[2],
                                                                                       ADbox[3], ADbox[4], ADbox[5],
                                                                                       int(5*2*8))

        with open(path + "/" + self.Ligand_Object.filename.split("/")[-1] + ".conf", "w") as f:
            print(path + "/" + self.Ligand_Object.filename + ".conf")
            f.write(text)
