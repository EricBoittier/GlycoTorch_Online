path = "/home/eric/Projects/GlycoTorch/data/ligands/pdb"
#path = '/Volumes/Eric/Projects/GlycoTorch/data/ligands/pdb'

import sys
sys.path.insert(0, '/home/eric/Projects/GlycoTorch')
sys.path.insert(0, '/Volumes/Eric/Projects/GlycoTorch')
from Carbohydrate import *

import matplotlib.pyplot as plt
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.image as mpimg


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

files = [x for x in os.listdir(path) if x.__contains__(".pdb")]

unique_linkages = set()

all_linkages = []

df = pd.DataFrame()

def color_point(ring):
    if ring.__contains__("B"):
        return "cornflowerblue"
    if ring.__contains__("E"):
        return "darkslateblue"
    if ring.__contains__("S"):
        return "gold"
    if ring.__contains__("H"):
        return "firebrick"
    if ring.__contains__("4C"):
        return "forestgreen"
    if ring.__contains__("1C"):
        return "grey"
    else:
        return "k"

from matplotlib.patches import Patch
from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], marker='o', color='forestgreen', label=r'$^{4}C_{1}$',
                          markerfacecolor='forestgreen', markersize=5, linestyle='None'),
                   Line2D([0], [0], marker='o', color='grey', label=r'$^{4}C_{1}$',
                          markerfacecolor='grey', markersize=5, linestyle='None'),
                   Line2D([0], [0], marker='o', color='gold', label=r'$^{x}S_{y}$',
                          markerfacecolor='gold', markersize=5, linestyle='None'),
                   Line2D([0], [0], marker='o', color='cornflowerblue', label=r'$^{x}B_{y}$',
                          markerfacecolor='cornflowerblue', markersize=5, linestyle='None'),
                   Line2D([0], [0], marker='o', color='darkslateblue', label=r'$^{x}E_{y}$',
                          markerfacecolor='darkslateblue', markersize=5, linestyle='None'),
                   Line2D([0], [0], marker='o', color='firebrick', label=r'$^{x}H_{y}$',
                          markerfacecolor='firebrick', markersize=5, linestyle='None')]


# plt.xlim(0, 360)
# plt.ylim(0, 180)

for file in files:
    ligand = Carbohydrate(path+"/"+file)
    for r in ligand.Rings.values():

        PHI = r.theta_radians
        THETA = r.phi_radians
        R = 1
        X = R * np.sin(PHI) * np.cos(THETA)
        Y = R * np.sin(PHI) * np.sin(THETA)
        Z = R * np.cos(PHI)

        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_zlim(-1, 1)
        ax.scatter(X, Y, Z,
                    color=color_point(r.iupac_name), alpha=0.5)




plt.gca().invert_yaxis()

tropic_of_cancer = [r"$^{0}E$", r"$^{O}H_{1}$", r"$E_{1}$", r"$^{2}H_{1}$", r"$^{2}E$",
r"$^{2}H_{3}$", r"$E_{3}$", r"$^{4}H_{3}$", r"$^{4}E$",
r"$^{4}H_{5}$", r"$E_{5}$", r"$^{O}H_{5}$", r"$^{O}E$"]

meridian = [r"$^{3,O}B$", r"$^{3}S_{1}$", r"$B_{1,4}$", r"$^{5}S_{1}$", r"$^{2,5}B$",
r"$^{2}S_{O}$",r"$B_{3,O}$", r"$^{1}S_{3}$", r"$^{1,4}B$", r"$^{1}S_{5}$", r"$H_{2,5}$",
r"$^{O}S_{2}$",r"$^{3,O}B$",]

tropic_of_capricorn = [r"$^{3}E$", r"$^{3}H_{4}$", r"$E_{4}$", r"$^{5}H_{4}$", r"$^{5}E$",
r"$^{5}H_{O}$", r"$E_{O}$", r"$^{1}H_{O}$", r"$^{1}E$",
r"$^{1}H_{2}$", r"$E_{2}$", r"$^{1}H_{5}$", r"$^{3}E$"]
#
# ax = plt.gca()
# ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handles=legend_elements)
ax.legend(handles=legend_elements)
# circle_rad = 5
# for number, label in enumerate(tropic_of_cancer):
#     plt.text(number*30, 50, 1, label)
#
# ax.axhline(45, color="grey", alpha=0.5)
#
# for number, label in enumerate(meridian):
#     plt.text(number*30-5, 95, 1, label)
#
# ax.axhline(95, color="grey", alpha=0.5)
#
# for number, label in enumerate(tropic_of_capricorn):
#     plt.text(number*30, 140, 1, label)
#
# ax.axhline(135, color="grey", alpha=0.5)
#
# for x in range(len(tropic_of_cancer)):
#     ax.axvline(x*30, color="grey", alpha=0.5)

#ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

ax.set_aspect("equal")

# draw sphere
u, v = np.mgrid[0:2*np.pi:10j, 0:np.pi:10j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.set_axis_off()
ax.plot_wireframe(x, y, z, edgecolor="k", alpha=0.5, linewidth=0.2)

plt.yticks(range(0, 181, 45), ['{}'.format(x) for x in range(0, 181, 45)])
plt.xticks(range(0, 361, 30), ['{}'.format(x) for x in range(0, 361, 30)])
# plt.xlabel("")
# plt.ylabel("")
# plt.text(172, -1, s= 1, r'$^{4}C_{1}$')
# plt.text(172, 180, 1, r'$^{1}C_{4}$')
plt.savefig('globe.pdf', bbox_inches="tight")
plt.show()

#plt.savefig("cremer_pople_flat.png", dpi=300, bbox_inches="tight")