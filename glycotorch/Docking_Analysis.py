from glycotorch.Vina_Output import *
from glycotorch.RMSD_Vina import *
# import markov_clustering as mc
# import networkx as nx
import matplotlib.pyplot as plt
import os
import os.path
from matplotlib.pylab import show, cm, axis
from glycotorch.PDBQT_to_PDB import *

"""
Takes output structures, protein structure and crystal pose
"""
class docking_analysis(object):
    def __init__(self, protein, output, crystal_pose=None, location="."):
        """
        :param protein: protein in .pdbqt format
        :param output: output file in .pdbqt format
        :param crystal_pose: crystal structure ligand (if exists) in .pdbqt format
        """
        self.location = location
        self.protein_path = protein
        self.output_path = output
        self.crystal_pose_path = crystal_pose
        self.rmsd = {}
        self.pairwise_rmsd = nx.Graph()
        self.markov_rmsd = nx.Graph()

        output_file = OutputFile(self.output_path)
        output_file.make_pdbqt_files(location=location)
        self.do_known_rmsd_calculations()

        # PDBQT_to_PDB(self.protein_path)
        # merge = PDB_Merge()
        # merge.add_ligand(self.crystal_pose_path)
        # merge.add_protein(self.protein_path)
        # merge.merge_protein_ligand(os.path.join(self.location, "complex.xyz"))
        #

    def do_known_rmsd_calculations(self):
        """

        :return:
        """
        ref_ligand = self.crystal_pose_path
        reference_ligand = PDBQT(ref_ligand)

        for ligand in os.listdir(self.location):
            if ligand.__contains__("MODEL_") and not ligand.__contains__(".pdbqt.pdb"):
                docked_ligand = PDBQT(self.location + ligand)
                self.rmsd[int(ligand.split("_")[1])] = [ligand, RMSD(reference_ligand, docked_ligand).get_ring_rmsd(),
                                                        docked_ligand.get_energy_result()]




        s = list(self.rmsd.keys())
        s.sort()

        rmsd_data = False
        rmsd_values = [x[1] for x in self.rmsd.values()]

        if None not in rmsd_values:
            rmsd_data = []
            binding_data = []
            new_file = open(os.path.join(self.location, "RMSD.csv"), "w")
            new_file.write("Model,RMSD,Energy\n")
            count = 0
            for key in s:
                if count < 10:
                    rmsd_data.append(self.rmsd[key][1])
                    binding_data.append(self.rmsd[key][0])
                    count += 1
                new_file.write("{0},{1:.2f},{2:.2f}\n".format(key, self.rmsd[key][1], self.rmsd[key][2]))

        if rmsd_data and len(rmsd_data) > 10:
            title = ref_ligand.split("/")[-1].replace("_", " ")

            fig, axs = plt.subplots(1, 2, sharey=True, gridspec_kw={'width_ratios': [1, 3]})
            plt.subplots_adjust(wspace=0.00001)

            fig.suptitle(title)

            axs[1].bar(range(0, 10), rmsd_data[0:11])
            axs[1].set_xticks(range(0, 10))
            axs[1].set_xticklabels(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"])
            box_plot = [float(x) for x in rmsd_data]

            axs[0].boxplot(box_plot[:])

            axs[1].set(xlabel='Docking pose')
            axs[0].set(ylabel=r'RMSD ($\AA$) of ring atoms', xticks=[])

            axs[0].spines["top"].set_visible(False)
            axs[0].spines["right"].set_visible(False)
            axs[0].spines["bottom"].set_visible(False)

            plt.savefig(self.location + "{}_barchart_.png".format(title))

            plt.clf()

            plt.scatter(rmsd_data, binding_data)
            plt.xlabel = "RMSD"
            plt.ylabel = "Binding energy"

            plt.savefig(self.location + "RMSD_vs_Binding.png")


    def do_pairwise_rmsd_calculations(self):
        """

        :return:
        """

        def draw_graph(matrix, clusters, **kwargs):
            """
            Visualize the clustering

            :param matrix: The unprocessed adjacency matrix
            :param clusters: list of tuples containing clusters as returned
                             by 'get_clusters'
            :param kwargs: Additional keyword arguments to be passed to
                           networkx.draw_networkx
            """
            # make a networkx graph from the adjacency matrix
            graph = nx.Graph(matrix)

            # map node to cluster id for colors
            cluster_map = {node: i for i, cluster in enumerate(clusters) for node in cluster}
            colors = [cluster_map[i] for i in range(len(graph.nodes()))]

            # if colormap not specified in kwargs, use a default
            if not kwargs.get("cmap", False):
                kwargs["cmap"] = cm.tab20

            # draw
            nx.draw_networkx(graph, node_color=colors, **kwargs)
            axis("off")
            plt.savefig(self.location + "MarkovClustering.pdf", bbox_inches='tight')

        labels = {}

        int_to_model = {}

        labels["A"] = "A"
        int_to_model["A"] = self.crystal_pose_path

        for ligand1 in os.listdir(self.location):
            for ligand2 in os.listdir(self.location):
                if ligand1.__contains__("MODEL_") and not ligand1.\
                        __contains__(".pdbqt.pdb") and ligand2.__contains__("MODEL_") and not ligand2.\
                        __contains__(".pdbqt.pdb"):

                    int_to_model[int(ligand1.split("_")[1])] = self.location+ligand1

                    try:
                        reference_ligand = PDBQT(self.location + ligand1)
                        docked_ligand = PDBQT(self.location + ligand2)


                        if RMSD(reference_ligand, docked_ligand):
                            self.pairwise_rmsd.add_node(int(ligand1.split("_")[1]))
                            self.pairwise_rmsd.add_node(int(ligand2.split("_")[1]))

                            labels[int(ligand1.split("_")[1])] = int(ligand1.split("_")[1])


                            self.pairwise_rmsd.add_edge(int(ligand1.split("_")[1]), int(ligand2.split("_")[1]),
                                                        weight=RMSD(reference_ligand, docked_ligand))

                        if RMSD(reference_ligand, docked_ligand) < 3:
                            self.markov_rmsd.add_node(int(ligand1.split("_")[1]))
                            self.markov_rmsd.add_node(int(ligand2.split("_")[1]))
                            self.markov_rmsd.add_edge(int(ligand1.split("_")[1]), int(ligand2.split("_")[1]),
                                                      weight=RMSD(reference_ligand, docked_ligand))


                    except ZeroDivisionError or KeyError as e:
                        pass


        ref_ligand = self.crystal_pose_path
        reference_ligand = PDBQT(ref_ligand)
        PDBQT_to_PDB(self.crystal_pose_path)

        for ligand in os.listdir(self.location):
            if ligand.__contains__("MODEL_") and not ligand.\
                        __contains__(".pdbqt.pdb"):
                try:

                    docked_ligand = PDBQT(self.location + ligand)

                    if RMSD(reference_ligand, docked_ligand):
                        self.pairwise_rmsd.add_node("A")
                        self.pairwise_rmsd.add_node(int(ligand.split("_")[1]))
                        self.pairwise_rmsd.add_edge("A", int(ligand.split("_")[1]),
                                                    weight=RMSD(reference_ligand, docked_ligand))

                    if RMSD(reference_ligand, docked_ligand) < 2:
                        self.markov_rmsd.add_node("A")
                        self.markov_rmsd.add_node(int(ligand.split("_")[1]))
                        self.markov_rmsd.add_edge("A", int(ligand.split("_")[1]),
                                                  weight=RMSD(reference_ligand, docked_ligand))

                except ZeroDivisionError or KeyError as e:
                    pass

        labels = {}
        for idx, node in enumerate(self.markov_rmsd.nodes()):
            labels[idx] = node

        pairwise = nx.to_pandas_adjacency(self.pairwise_rmsd)
        pairwise.to_csv(self.location+"pairwise.csv")

        if len(self.markov_rmsd.edges) > 2:
            matrix = nx.to_scipy_sparse_matrix(self.markov_rmsd, weight='weight',)

            markov = nx.to_pandas_adjacency(self.markov_rmsd)
            markov.to_csv(self.location+"markov.csv")

            inflation = 2
            plt.figure(figsize=(22, 22))
            result = mc.run_mcl(matrix, inflation=inflation)
            clusters = mc.get_clusters(result)

            ccc = clusters
            ccc.sort(key = lambda s: -len(s))

            count = 1

            # if not os.path.exists(self.protein_path + '.pdb'):
            #     PDBQT_to_PDB(self.protein_path)

            for cluster in ccc[0:3]:
                f = open(self.location+"/cluster_{}".format(count), "w")
                for c in cluster:
                    if not os.path.exists(int_to_model[labels[c]] + '.pdb') and not int_to_model[labels[c]].__contains__("out.pdbqt.pdb"):
                        PDBQT_to_PDB(int_to_model[labels[c]])
                    f.write(int_to_model[labels[c]].split("/")[-1] + '.pdb' + "\n")
                count += 1

            draw_graph(matrix, clusters, labels=labels, node_size=200, with_labels=True, edge_color="#DCDCDC",
                       linewidths=0.001)



# /home/eric/Projects/GlycoTorch/data/results/_1BFB_WithoutLigands.pdb.pdbqt

#
# da = docking_analysis("/home/eric/Projects/GlycoTorch/data/results/glycosidic/0_0/small/_glycosidic_1BFB_LIGAND_1_uncharged/_1BFB_WithoutLigands.pdb.pdbqt",
#                  "/home/eric/Projects/GlycoTorch/data/results/glycosidic/0_0/small/_glycosidic_1BFB_LIGAND_1_uncharged/_glycosidic_1BFB_LIGAND_1_uncharged.pdb.mol2.pdbqt",
#                  crystal_pose="/home/eric/Projects/GlycoTorch/data/results/glycosidic/0_0/small/_glycosidic_1BFB_LIGAND_1_uncharged/_glycosidic_1BFB_LIGAND_1_uncharged.pdb.mol2.pdbqt",
#                  location="/home/eric/Projects/GlycoTorch/data/results/glycosidic/0_0/small/_glycosidic_1BFB_LIGAND_1_uncharged/")
#
# da.do_known_rmsd_calculations()
