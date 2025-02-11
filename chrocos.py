""" Chromatic Community Structure (ChroCos) analysis module.
Module content:
    'ChroCoDe',
    'CommunityColorProfile',
    'DominantSigs',
    'DrawChroCoS',
    'DrawColoredGraph',
    'Gamma',
    'GenerateSeeds',
    'H',
    'Hcore',
    'Kappa',
    'MonochromeCommunityStructure',
    'RandomColoring'"""

# ** CHROMATIC COMMUNITY STRUCTURE
# AUTHOR : Franck Delaplace
# CREATION DATE: 12/03/2023

# ! DATA STRUCTURES
# a community p is a (frozen) set of nodes p=frozenset({n1,n2,...})
# a community structure P is a set of frozen node sets:
# P=[p1,p2,...] where pi is a frozenset.
# a color profile c is defined by a dictionary {node:color ...}.

# Throughout the functions:
#   - r is the number of colors, r>=0
#   - n is the size of the community,
#   - d is the number of nodes with the same color. 0<=d<=n
# The colors are integers from 1 to r and 0 stands for the transparent color

# ** Import functions & packages, typing ====================================
from math import factorial, comb, ceil, exp, inf, log2
from typing import TypeAlias, Callable
from random import choices, random
from collections import Counter
from functools import reduce
from scipy.stats import gmean  # type: ignore

import networkx as nx
import seaborn as sns  # type: ignore

# graph type alias
graph_t: TypeAlias = nx.classes.graph.Graph

# enumeration function type
fun3int2int_t: TypeAlias = Callable[[int, int, int], int]


# ** BASIC FUNCTIONS =====================================================


def CommunityColorProfile(p: frozenset, c: dict, basic=True) -> dict:
    """Define a color profile restricted to a community

    Args:
        p (frozenset): community
        c (dict): color profile
        basic (bool, optional): Boolean setting whether the transparent color node
                                should be removed (False) or not (True). Defaults to True.

    Returns:
        dict: _description_
    """
    return {k: v for (k, v) in c.items() if k in p and (basic or v != 0)}


def DominantSigs(r: int, n: int, d: int) -> list[list[int]]:
    """compute all the d-dominant signatures (list) for community of size n considering r colors.

    Args:
        r (int): number of colors
        n (int): community size
        d (int): dominant color size

    Returns:
        list: all d-dominant signatures
    """
    assert r >= 0
    assert n >= d >= 0

    # sub function of Dominant signatures
    # depth=depth, rest=remaining value for which elements to sum must be found,
    # dbound= max bound value, sig= signature
    def DSS(depth: int, rest: int, dbound: int, sig: list):
        if rest == 0:
            ds = [depth * [0] + sig]
        elif dbound == 0:
            ds = []
        elif depth == 0:
            if rest == 0:
                ds = [sig]
            else:
                ds = []
        else:
            ds = []
            for v in range(ceil(rest / depth), min(dbound, rest) + 1):
                ds = DSS(depth - 1, rest - v, v, [v] + sig) + ds
        return ds

    # main
    domsigset = []  # type:list[list[int]]
    if r == 0:
        if n == 0:
            domsigset = [[]]
    else:
        domsigset = DSS(r - 1, n - d, d, [d])

    return domsigset


# ** CHROMARITIES ===========================================================

# Enumeration --------------------------------------------------


def Kappa(r: int, n: int, d: int) -> int:
    """Count the number of d-coloring profiles of size n considering r colors.

    Args:
        r (int): number of colors
        n (int): size of the community
        d (int): number of nodes of the same color

    Returns:
        int: number of communities complying with these requirements.
    """
    kappa = 0
    for k in range(1, min(r, n // d) + 1):
        kappa += (
                         (-1) ** (k - 1) * comb(r, k) * factorial(n) * (r - k) ** (n - k * d)
                 ) // (factorial(n - k * d) * factorial(d) ** k)
    return kappa


def Gamma(r: int, n: int, d: int) -> int:
    """Count the number of d-dominant coloring profiles of size n considering r colors.

    Args:
        r (int): number of colors
        n (int): size of the community
        d (int): number of nodes of the same color which dominates in cardinality

    Returns:
        int: number of communities complying with these requirements.
    """
    gamma = 0
    factorialnr = factorial(n) * factorial(r)

    for sig in DominantSigs(r, n, d):
        ps = reduce(lambda x, y: x * y, map(factorial, sig), 1)
        pcs = reduce(lambda x, y: x * y, map(factorial, Counter(sig).values()), 1)
        gamma += factorialnr // (ps * pcs)
    return gamma


#  Chromatic entropy -------------------------------------


def Hcore(r: int, n: int, d: int, funenum: fun3int2int_t = Gamma) -> float:
    """Compute the core chromarity

    Args:
        r (int): number of colors,
        n (int): size of a community,
        d (int): number of nodes with the same color,
        funenum (function): enumeration function. (Default Gamma)

    Returns:
        float: chromarity value.
    """
    assert r > 0
    assert n >= d >= 0
    p = funenum(r, n, d) / r ** n

    try:
        h = -p * log2(1 - p)
    except ValueError:  # case if a community is empty  skip it
        h = inf

    return h


def H(P: set, c: dict, r: int, funenum: fun3int2int_t = Gamma) -> float:
    """compute the chromarity

    Args:
        funenum:
        P (set):  community structure,
        c (dict):  color profile,
        r (int): number of colors,

    Returns:
        float: chromatic entropy.
    """

    h = 0.0
    for p in P:
        colors = CommunityColorProfile(p, c, basic=False).values()
        try:
            h += Hcore(r, len(p), max(Counter(colors).values()), funenum)
        except (
                ValueError
        ):  # case if a community is empty or only with transparent nodes - skip it (provisionally )
            pass

    return h


# ** GRAPH =============================================================================

# Display  ----------------------------------------------------------


# Default palette
__CHROCOS_PALETTE__ = {
    0: "gainsboro",
    1: "lightgreen",
    2: "crimson",
    3: "gold",
    4: "steelblue",
    5: "mediumpurple",
    6: "darkorange",
    7: "burlywood",
    8: "salmon",
    9: "orchid",
    10: "teal",
}

# Default font
__FONT__ = "Franklin Gothic Heavy"  # other nice fonts  'Tahoma'  'Impact'


def DrawColoredGraph(G, palette: dict[int, str] = __CHROCOS_PALETTE__, pos=None, color="color"):
    """Display a colored graph
    :param G:  Graph
    :type G:  nx.Graph
    :param palette:  color palette associating a color to integer. Defaults to __CHROCOS_PALETTE__ (up to 10 colors)
    :type palette: dict[int, str]
    :param pos: position where to draw the graph (Defaults: None)
    :type pos: dict|None
    :param  color: attribute denominating the color (Defaults: 'color')
    :type color: str
    :return: Draw graph
    :rtype: None
    """
    color = nx.get_node_attributes(G, name=color)
    nx.draw_networkx(
        G,
        with_labels=True,
        pos=pos,
        node_color=[palette[c] for c in color.values()],
        node_size=500,
        font_size=11,
        font_color="black",
        font_family=__FONT__,
    )


# Display the community structure on graph
def DrawChroCoS(G, P: set[frozenset], theme: str = "Set2", pos=None):
    """Display the chromarity structure on a graph. The nodes of the same community are the same color.
    Args:
        G (Graph): Undirected colored graph_
        P (set[frozenset]): community structure_
        theme (str, optional): theme color of seaborn package. Defaults to 'Set2'.
        pos (dict|None, optional): position of nodes. Defaults to None.
    """

    Pl = list(P)
    palette = sns.color_palette(theme, len(Pl))
    color = {v: palette[i] for i in range(len(Pl)) for v in Pl[i]}  # type: ignore
    nx.draw_networkx(
        G,
        with_labels=True,
        pos=pos,
        node_color=[color[v] for v in G.nodes()],
        node_size=500,
        font_size=11,
        font_color="black",
        font_family=__FONT__,
    )


# Random Graph  --------------------------------------------------


def RandomColoring(
        G: graph_t, seeds: list, density: float = 0.2, transparency: float = 0.0, color: str = "color"
):
    """Attributes colors to nodes of graph G randomly.
    :param G: Graph
    :type G: nx.Graph
    :param seeds: list of seeds
    :type seeds: list
    :param density:  probability parameter:  higher the value less the colors are scattered .
    (Defaults: 0.2).
    :type density: float
    :param transparency: probability of transparent nodes.
    (Defaults: 0).
    :type transparency: float
    :param color: attribute denominating the color (Default :'color')
    :type color: str
    :return: a colored graph
    :rtype: nx.Graph with color attributes
    """
    assert 0.0 <= density <= 1.0
    assert 0.0 <= transparency <= 1.0
    assert nx.is_connected(G)  # require that a path exists for any pair of vertices.

    # sub function selecting the color from seeds w.r.t. an exponential probability
    # the farther from a color seed the node, the lower the probability of the color is.
    def ChooseColorRandomly(seeds: list, v) -> int:
        return choices(
            range(1, len(seeds) + 1),
            weights=[
                exp(-density * nx.shortest_path_length(G, seed, v)) for seed in seeds
            ],
            k=1,
        )[0]

    nx.set_node_attributes(
        G, dict([(v, ChooseColorRandomly(seeds, v)) for v in G.nodes()]), name=color)

    if transparency > 0:
        transparent = [v for v in G.nodes() if random() < transparency]
        nx.set_node_attributes(G, dict.fromkeys(transparent, 0), name=color)


def GenerateSeeds(G: graph_t, r: int) -> list:
    """Generate r color seeds for graph G by maximizing the  geometric mean distance between them.

    Args:
        G (Graph): undirected graph that must be a single component.
        r (int): number of colors with r>1
    Returns:
    list: color seeds
    """
    assert r > 1
    assert nx.is_connected(G)  # require that a path exists for any pair of vertices.
    pathlength = dict(nx.all_pairs_shortest_path_length(G))

    # first select two nodes with a maximal distance
    maxi = 0
    v = w = vmax = wmax = next(iter(pathlength))
    for v, vpathlength in pathlength.items():
        for w, distance in vpathlength.items():
            if distance > maxi:
                maxi = distance
                vmax = v
                wmax = w
    seeds = [vmax, wmax]

    # Next complete the seed list with nodes maximizing the mean geometric distance to each already selected seed
    V = list(G.nodes())
    V.remove(vmax)
    V.remove(wmax)
    for _ in range(r - 2):
        maxi = 0
        vmax = None
        for v in V:
            disttoseeds = [
                pathlength[v][seed] for seed in seeds if seed in pathlength[v]
            ]

            if disttoseeds:
                meandisttoseeds = gmean(disttoseeds)
                if meandisttoseeds > maxi:
                    maxi = meandisttoseeds
                    vmax = v
        seeds.append(vmax)
        V.remove(vmax)
    return seeds


# ** CHROCODE ===========================================================================


def MonochromeCommunityStructure(G: graph_t, color: str) -> set:
    """Compute a monochrome community structure of graph G.
    :param G: graph
    :type G: nx.Graph
    :param color: color demonating the color
    :type color: str
    :return: community structure
    :rtype:  set[frozenset[node]]
    """
    pendingnodes = set(G.nodes())
    P = set()
    while pendingnodes:
        v = next(iter(pendingnodes))  # take the first available node
        referencecolor = G.nodes[v][color]  # get its color used as reference color
        p = set()
        monochromeneighbors = {v}
        while monochromeneighbors:
            v = monochromeneighbors.pop()  # take a node for including its neighbors
            p.add(v)  # add the current node to the community
            for w in G[v]:  # extend the neighbor set with the same color
                if G.nodes[w][color] == referencecolor and w not in p:
                    monochromeneighbors.add(w)

        pendingnodes = (
                pendingnodes - p
        )  # remove the community of the examined vertices.
        P.add(frozenset(p))  # add the community to the structure
    return P


def ChroCoDe(
        G: graph_t, r: int, radius: int = 1, funenum: fun3int2int_t = Gamma, color: str = "color"
) -> set:
    """Find a chromatic community structure from a colored graph.
    :param G: undirected graph
    :type G: nx.Graph
    :param r: number of colors
    :type r: int
    :param radius: radius of community (Default: 1)
    :type radius: int
    :param funenum: function that returns the number of colors
    :type funenum:  function(int,int,int)->int
    :param color: attribute denominating the color (Default: 'color')
    :type color: str
    :return: community structure
    :rtype: set[frozenset[node]]
    """
    assert radius > 0
    assert r > 0

    colorprofile = nx.get_node_attributes(G, name=color)
    QG = nx.quotient_graph(
        G, MonochromeCommunityStructure(G, color)
    )  # Quotient graph of the monochrome community structure.
    P = set(QG.nodes())
    Pscan = P.copy()  # Initialize the running community structure Pscan.

    while Pscan:  # When no improvements of P occur then Pscan will empty progressively.
        hmax = -1.0
        p = set()
        for q in Pscan:  # find the community p in P with the maximal chromatic entropy.
            h = H({q}, colorprofile, r, funenum)
            if h > hmax:
                hmax = h
                p = q
        Pscan.remove(p)  # remove p of the running community structure PScan.

        N = list(nx.ego_graph(QG, p, radius=radius, center=False).nodes())  # compute the neighborhood of size radius but p.

        hmin = H(P, colorprofile, r, funenum)
        minpath = set()
        improved = False
        for q in N:  # find the neighbor q of p minimizing H by merging the path p-q.
            communitypath = set(nx.shortest_path(QG, p, q))
            pmerge = reduce(lambda p, q: p | q, communitypath)  # Union of communities of the path.

            h = H((P - communitypath) | {pmerge}, colorprofile, r, funenum)  # add pmerge to P
            if hmin >= h:
                improved = True
                hmin = h
                minpath = communitypath.copy()

        if (
                improved
        ):  # Update P with the shortest path connecting p to this neighbor by merging their community.
            pmerge = reduce(lambda p, q: p | q, minpath)
            P = (P - minpath) | {pmerge}
            QG = nx.quotient_graph(G, P)  # re-compute the quotient graph QG
            Pscan = P.copy()  # re-initialize Pscan
    return P
