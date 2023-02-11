from glycotorch.Dihedral import *
from glycotorch.PDB_Abbreviation_to_Sugar_Name import *


class Linkage(object):
    """docstring for Linkage"""

    def __init__(self, O5, C1, GO, CX, CX_minus_1):
        super(Linkage, self).__init__()
        self.O5 = O5
        self.C1 = C1
        self.GO = GO
        self.CX = CX
        self.CX_minus_1 = CX_minus_1
        self.phi = None
        self.set_phi()
        self.psi = None
        self.set_psi()
        self.linkage_type = None
        self.linkage_number = None
        self.linkage_name = None
        self.ring1 = None
        self.ring2 = None

    def set_rings(self, ring1, ring2):
        self.ring1 = ring1
        self.ring2 = ring2
        self.set_linkage_type()
        self.set_linkage_name()


    def get_linkage_name(self):
        return self.linkage_name

    def set_linkage_type(self):
        self.linkage_number = int(self.CX.sugar_position[-1])
        link = "1-{}".format(self.CX.sugar_position[-1])

        anomer1 = ""
        if self.ring1.iupac_name.__contains__("C"):
            anomer1 = self.is_alpha_or_beta(self.C1.get_is_axial())

        anomer2 = ""
        if self.ring2.iupac_name.__contains__("C"):
            anomer2 = self.is_axial_or_eq(self.CX.get_is_axial())

        self.linkage_type = "({}{}{})".format(anomer1, link, anomer2)

    def set_linkage_name(self):
        self.linkage_name = "{}-{}-{}".format(abbreviationToSugarName(self.C1.getLigandType()),
                                              self.linkage_type,
                                              abbreviationToSugarName(self.CX.getLigandType()))

    def get_linkage_type(self):
        return self.linkage_type

    def set_psi(self):
        """

        C1 - O - Cx - Cx-1

        Psi is defined as the dihedral angle between C1, the oxygen of the
        glycosidic bond, CX (where x is the ring carbon number connected to the
        glycosidic oxygen), and CX-1 is the carbon numbered 1 less than the CX.

        :return:
        """
        self.psi = float(dihedral(self.C1.getXYZ(), self.GO.getXYZ(), self.CX.getXYZ(), self.CX_minus_1.getXYZ()))

    def set_phi(self):
        """
        O5 - C1 - O - Cx

        Psi is defined as the dihedral angle between the ring oxygen, carbon 1,
        the oxygen of the glycosidic bond,CX (where x is the ring carbon number
        connected to the glycosidic oxygen).

        :return:
        """
        self.phi = float(dihedral(self.O5.getXYZ(), self.C1.getXYZ(), self.GO.getXYZ(), self.CX.getXYZ()))



    def get_phi(self):
        return self.phi

    def get_psi(self):
        return self.psi

    def is_axial_or_eq(self, axial):
        if not axial:
            return "e"
        else:
            return "a"

    def is_alpha_or_beta(self, alpha):
        if not alpha:
            return "β"
        else:
            return "α"
