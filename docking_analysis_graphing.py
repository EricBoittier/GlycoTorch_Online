import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
from statistics import *

files = os.listdir("data/results/rigid/standard/small")

current_palette = sns.color_palette("muted")
sns.palplot(current_palette)

my_pal = {"standard": "#2471A3", "gauss1_0": "#F1948A", "gauss2_0": "#E74C3C", "hbond_0": "#A93226",
          "hydrophobic_0": "#ABEBC6",
          "repulsive_0": "#27AE60", "rot_0": "#117A65"}

print(files)

scores = {"standard": 0, "gauss1_0": 0, "gauss2_0": 0, "hbond_0": 0, "hydrophobic_0": 0,
          "repulsive_0": 0, "rot_0": 0}
scorer = []
all_rmsd_list = {"standard": [], "gauss1_0": [], "gauss2_0": [], "hbond_0": [], "hydrophobic_0": [],
                 "repulsive_0": [], "rot_0": []}

all_rmsd = {"standard": 0, "gauss1_0": 0, "gauss2_0": 0, "hbond_0": 0, "hydrophobic_0": 0,
            "repulsive_0": 0, "rot_0": 0}

for name in files:
    pass_ = True

    if pass_:
        try:
            print(name)
            results = []

            result_to_key = {}

            fig, ax = plt.subplots()

            chi_cutoff_co_ef_values = ["gauss1_0", "gauss2_0", "hbond_0", "hydrophobic_0", "repulsive_0", "rot_0",
                                       "standard"]

            for chi_chi in ["gauss1_0", "gauss2_0", "hbond_0", "hydrophobic_0", "repulsive_0", "rot_0", "standard"]:

                f = open("data/results/rigid/{}/small/{}/RMSD.csv".format(chi_chi, name))

                lines = f.readlines()
                model_rmsd = {}
                for line in lines[1:]:
                    model_rmsd[int(line.split(",")[0])] = float(line.split(",")[1])

                top_5 = []

                for k in range(5):
                    top_5.append(model_rmsd[k + 1])

                results.append(min(top_5))

                result_to_key[min(top_5)] = chi_chi

                all_rmsd[chi_chi] += min(top_5)
                all_rmsd_list[chi_chi].append(min(top_5))

            scorer.append(result_to_key)

            for x in range(len(chi_cutoff_co_ef_values)):
                plt.bar(x, results[x])

            ax = sns.barplot(chi_cutoff_co_ef_values, results, palette=current_palette)

            ax.set_xticklabels(chi_cutoff_co_ef_values)

            plt.xticks(range(7), [x.replace("_", " = ") for x in chi_cutoff_co_ef_values])
            plt.xticks(rotation=45)
            plt.title("{}".format(name))
            ax.set(xlabel='CHI coefficent, CHI cut-off', ylabel='RMSD (Å)')
            plt.clf()
        except Exception as e:
            print(e)

for dictionary in scorer:
    c = 1
    l = list(dictionary.keys())
    l.sort()
    for key in l:
        scores[dictionary[key]] += c
        c += 1

amount_of_results = len(scorer)
print(amount_of_results)

for key in scores:
    scores[key] = scores[key] / amount_of_results // 1

print(scores)

avg_rmsd = []
for key in chi_cutoff_co_ef_values:
    avg_rmsd.append(all_rmsd[key] / amount_of_results)

print(avg_rmsd)

ax = sns.barplot(chi_cutoff_co_ef_values, avg_rmsd, palette=current_palette)
# ax.set_xticklabels(chi_cutoff_co_ef_values)

plt.xticks(range(7), [x.replace("_", " = ") for x in chi_cutoff_co_ef_values])
plt.xticks(rotation=45)
plt.title(r"Average min$_{top 5}$(RMSD) (Å)")
ax.set(xlabel='Scoring function term', ylabel='RMSD (Å)')

