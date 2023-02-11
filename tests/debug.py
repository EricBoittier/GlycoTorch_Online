from glycotorch.Carbohydrate import Carbohydrate
from glycotorch.Carbohydrate_to_PDBQT import Carbohydrate_to_PDBQT

# c = Carbohydrate("../../Downloads/Unsulphated_HS_tetramer_glycam.pdb")
c = Carbohydrate("../../../Downloads/1-2.pdb")
print(c.ordered_linkages[0].__dict__)
for _ in c.Rings.values():
    print(_.print_functional_groups())
print(c.get_name())
pdbqt = Carbohydrate_to_PDBQT(c)
pdbqt.save_flex()
# print(c.getSNFGname())