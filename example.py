""" Example of ChroCoS application"""
# AUTHOR : Franck Delaplace
# CREATION DATE: 17/03/2023
# OBJECT:   this script shows the basic use of the ChroCoDe algorithm on different graphs
#           and draw the chromatic community structure solution on the graphs.
#
import networkx as nx
import matplotlib.pyplot as plt

from chrocos import (
    K,
    Kappa,
    Gamma,
    DrawColoredGraph,
    DrawChroCoS,
    RandomColoring,
    GenerateSeeds,
    MonochromeCommunityStructure,
    ChroCoDe
)



def printcommunities(P):
    """print the community structure

    Args:
        P (set[frozenset]): community structure.
    """
    for p in P:
        s = ""
        for v in p:
            s += " " + str(v)
        print("{", s, "}")


def graphexample(G, position, title, transparency=0.0):
    """ "Example of ChroCoDe computation on a graph.

    Args:
        G (Graph): initial graph
        position (dict): node position
        title (str): title of the figure
        transparency (float, optional): probability of transparent nodes. Defaults to 0.0.
    """

    seeds = GenerateSeeds(
        G, r
    )  # generate seeds - they represent the 'corners' of the grid graph.
    RandomColoring(
        G, seeds, density=0.3, transparency=transparency
    )  # color the graph randomly.
    cp = nx.get_node_attributes(G, "color")  # get the color profile.

    P0 = MonochromeCommunityStructure(G)
    plt.subplot(221)
    plt.title(
        title
        + " network: Monochrome community: Kk=%4.2f, Kg=%4.2f"
        % (K(P0, cp, 4, Kappa), K(P0, cp, 4, funK=Gamma))
    )
    DrawColoredGraph(G, pos=position)  # Display the graph

    # RADIUS = 1
    print("radius=1")
    P = ChroCoDe(G, r, radius=1, funK=Gamma)

    # print the community structure.
    print("P=")
    printcommunities(P)

    # Display the result
    plt.subplot(223)
    plt.title(
        "radius=1 - Kg="
        + "{:.2f}".format(K(P, cp, 4, Gamma))
        + " - Kk="
        + "{:.2f}".format(K(P, cp, 4, funK=Kappa))
    )
    DrawChroCoS(G, P, pos=position)  # Display the community structure on the graph

    # RADIUS = 2 (default value)
    print("radius=2")
    P = ChroCoDe(G, r, funK=Gamma)

    #   print the community structure.
    print("P=")
    printcommunities(P)

    # Display the result
    plt.subplot(224)
    plt.title(
        "radius=2 - Kg="
        + "{:.2f}".format(K(P, cp, 4, Gamma))
        + " - Kk="
        + "{:.2f}".format(K(P, cp, 4, funK=Kappa))
    )
    DrawChroCoS(
        G, P, theme="pastel", pos=position
    )  # Display the community structure on the graph

    # Show
    plt.show()


# MAIN ==============================================================================
# NOTE: Close the pyplot windows to pursue.

r = 4  # number of colors

# GRID GRAPH ===========================================================
print("GRID GRAPH")
plt.figure(figsize=(10, 10))  # set size of the output graphic view
n = 8
GD = nx.grid_2d_graph(n, n)
G = nx.convert_node_labels_to_integers(GD)  # rename the vertices as integers

gridposition = dict(zip(G, GD))  # define position as label of the initial graph
graphexample(G, gridposition, "Grid")

# SMALL WORLD  ===========================================================
print("SMALL WORLD")
plt.figure(figsize=(15, 10))  # set size of the output graphic view
n = 40
G = nx.connected_watts_strogatz_graph(n, 2, 0.6)

position = nx.circular_layout(G)
graphexample(G, position, "Small World")

# ERDOS RENY ================================================================
print("ERDOS RENY")
plt.figure(figsize=(15, 10))  # set size of the output graphic view
n = 40
G = nx.erdos_renyi_graph(n, 0.3)

position = nx.spring_layout(G)
graphexample(G, position, "Erdös Reny")

# SCALE FREE ================================================================
print("Scale Free")
plt.figure(figsize=(15, 10))  # set size of the output graphic view
n = 40
SFG = nx.DiGraph(nx.scale_free_graph(n))
SFG.remove_edges_from(nx.selfloop_edges(SFG))
G = SFG.to_undirected()

position = nx.spring_layout(G)
graphexample(G, position, "Scale Free")
# GRID with Transparent color ===============================================
print("GRID GRAPH WITH TRANSPARENCY")
plt.figure(figsize=(10, 10))  # set size of the output graphic view
n = 8
GD = nx.grid_2d_graph(n, n)
G = nx.convert_node_labels_to_integers(GD)  # rename the vertices as integers

gridposition = dict(zip(G, GD))  # define position as label of the initial graph
graphexample(G, gridposition, "Grid with transparency", transparency=0.2)
