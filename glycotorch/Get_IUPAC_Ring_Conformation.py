from glycotorch.Dihedral import *
import math

endocyclic_torsional_angles = {"1C4": [60.0, -60.0, 60.0, -60.0, 60.0, -60.0],
                               "4C1": [-60.0, 60.0, -60.0, 60.0, -60.0, 60.0],
                               "1,4B": [0.0, 60.0, -60.0, 0.0, 60.0, -60.0],
                               "B2,5": [60.0, 0.0, -60.0, 60.0, 0.0, -60.0],
                               "O,3B": [60.0, -60.0, 0.0, 60.0, -60.0, 0.0],
                               "B1,4": [0.0, -60.0, 60.0, 0.0, -60.0, 60.0],
                               "2,5B": [-60.0, 0.0, 60.0, -60.0, 0.0, 60.0],
                               "BO_3": [-60.0, 60.0, 0.0, -60.0, 60.0, 0.0],
                               "1S5": [30.0, 30.0, -60.0, 30.0, 30.0, -60.0],
                               "OS2": [60.0, -30.0, -30.0, 60.0, -30.0, -30.0],
                               "S1": [30.0, -60.0, 30.0, 30.0, -60.0, 30.0],
                               "5S1": [-30.0, -30.0, 60.0, -30.0, -30.0, 60.0],
                               "2SO": [-60.0, 30.0, 30.0, -60.0, 30.0, 30.0],
                               "1S3": [-30.0, 60.0, -30.0, -30.0, 60.0, -30.0],
                               "1H2": [45.0, -15.0, 0.0, -15.0, 45.0, -60.0],
                               "3H2": [60.0, -45.0, 15.0, 0.0, 15.0, -45.0],
                               "3H4": [45.0, -60.0, 45.0, -15.0, 0.0, -15.0],
                               "5H4": [15.0, -45.0, 60.0, -45.0, 15.0, 0.0],
                               "5HO": [0.0, -15.0, 45.0, -60.0, 45.0, -15.0],
                               "1HO": [15.0, 0.0, 15.0, -45.0, 60.0, -45.0],
                               "4H5": [-15.0, 45.0, -60.0, 45.0, -15.0, 0.0],
                               "OH5": [0.0, 15.0, -45.0, 60.0, -45.0, 15.0],
                               "OH1": [-15.0, 0.0, -15.0, 45.0, -60.0, 45.0],
                               "2H1": [-45.0, 15.0, 0.0, 15.0, -45.0, 60.0],
                               "2H3": [-60.0, 45.0, -15.0, 0.0, -15.0, 45.0],
                               "4H3": [-45.0, 60.0, -45.0, 15.0, 0.0, 15.0], "1E": [30.0, 0.0, 0.0, -30.0, 60.0, -60.0],
                               "E2": [60.0, -30.0, 0.0, 0.0, 30.0, -60.0], "3E": [60.0, -60.0, 30.0, 0.0, 0.0, -30.0],
                               "E4": [30.0, -60.0, 60.0, -30.0, 0.0, 0.0], "5E": [0.0, -30.0, 60.0, -60.0, 30.0, 0.0],
                               "EO": [0.0, 0.0, 30.0, -60.0, 60.0, -30.0], "4E": [-30.0, 60.0, -60.0, 30.0, 0.0, 0.0],
                               "E5": [0.0, 30.0, -60.0, 60.0, -30.0, 0.0], "OE": [0.0, 0.0, -30.0, 60.0, -60.0, 30.0],
                               "E1": [-30.0, 0.0, 0.0, 30.0, -60.0, 60.0], "2E": [-60.0, 30.0, 0.0, 0.0, -30.0, 60.0],
                               "E3": [-60.0, 60.0, -30.0, 0.0, 0.0, 30.0]}


def get_iupac_ring(C1, C2, C3, C4, C5, O):
    T1 = dihedral(C1, C2, C3, C4, negative=True)
    T2 = dihedral(C2, C3, C4, C5, negative=True)
    T3 = dihedral(C3, C4, C5, O, negative=True)
    T4 = dihedral(C4, C5, O, C1, negative=True)
    T5 = dihedral(C5, O, C1, C2, negative=True)
    T6 = dihedral(O, C1, C2, C3, negative=True)

    angles = [T1, T2, T3, T4, T5, T6]

    values = {}

    for key in endocyclic_torsional_angles.keys():
        dif_squared = 0
        for torsion in range(0, 6):
            dif_squared += (angles[torsion] - endocyclic_torsional_angles[key][torsion]) **2

        answer = math.sqrt(dif_squared)
        values[answer] = key

    keys = list(values.keys())

    return values[min(keys)]
