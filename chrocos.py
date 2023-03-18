# CHOMATIC COMMUNITY STRUCTURE
# AUTHOR : Franck Delaplace
# CREATION DATE: 12/03/2023

# DATA STRUCTURES
# a community p is a (frozen) set of nodes p={n1,n2,...} or p=frozenset({n1,n2,...})
# a community structure P is a set of frozen node sets  P=[p1,p2,...] where pi is a frozenset.
# a color profile c is defined by a dictionary {node:color ...}.

        # Throughout the functions:
#   - r is the number of colors, r>=1
#   - n is the size of the community, 
#   - d is the number of nodes with the same color. 0<=d<=n
# The colors are integers from 1 to r and 0 stands for the transparent color

# Import functions & packages ======================================================
from math import dist, factorial,comb,ceil,exp,inf 
from scipy.stats import gmean
from random import choices, gammavariate, random
from collections import Counter
from functools import reduce

import networkx as nx
import seaborn as sns

# Basic Functions ==================================================================

def CommunityColorProfile(p,c):
    """CleanCommunityColorProfile(p,c): Define a color profile restricted to a community

    Args:
        p (set): community
        c (dict): color profil

    Returns:
    dict: comunity color profile.
    """
    return {k:v  for (k,v) in c.items() if k in p}

def CleanCommunityColorProfile(p,c):
    """CleanCommunityColorProfile(p,c): Define a color profile restricted to a community and remove the transparent profile (v:0)

    Args:
        p (set): community
        c (dict): color profil

    Returns:
    dict: comunity color profile.
    """
    return {k:v  for (k,v) in c.items() if k in p and v!=0}
      

def DominantSigs(r:int,n:int,d:int)->list:
    """compute all the d-dominant signatures (list) for community of size n considering r colors.

    Args:
        r (int): number of colors 
        n (int): community size
        d (int): dominant color size

    Returns:
        list: all d-dominant signatures 
    """
    assert r>= 0
    assert n>=d>=0
    # sub function of Dominant signatures (rest=remaining sum, dbound= bound value, seq= sequence of sum)
    def DSS(depth:int,rest:int,dbound:int,seq):
        if rest==0:
            ds=[depth*[0]+seq]
        elif dbound==0:
            ds=[]
        elif depth==0:
            if rest==0:
                ds=[seq]
            else:
                ds=[]
        else:
            ds=[]
            for v in range(ceil(rest/depth),min(dbound,rest)+1):
                ds= DSS(depth-1,rest-v,v,[v]+seq) + ds
        return ds
    # main
    if r==0:
        if n==0:
            dss= [[]]
        else:
            dss= []
    else:
        dss=DSS(r-1,n-d,d,[d])
    return dss

# CHROMARITIES =========================================================================

# Enumeration --------------------------------------------------


def Kappa(r:int,n:int,d:int) -> int:
    """Count the number of d-coloring profiles of size n considering r colors.

    Args:
        r (int): number of colors
        n (int): size of the community 
        d (int): number of nodes of the same color

    Returns:
        int: number of communities complying with these requirements.
    """
    assert r >0
    assert n >= d >=0
    kappa=0
    for k in range(1, min(r,n//d)+1):
        kappa+= ((-1)**(k-1)*comb(r,k)*factorial(n)*(r-k)**(n-k*d)) // (factorial(n-k*d)*factorial(d)**k)
    return kappa


def Gamma(r:int,n:int,d:int)->int:
    """Count the number of d-dominant coloring profiles of size n considering r colors.

    Args:
        r (int): number of colors
        n (int): size of the community 
        d (int): number of nodes of the same color which dominates in cardinality

    Returns:
        int: number of communities complying with these requirements.
    """
    assert r >0
    assert n >= d >=0
    gamma=0
    factorialnr=factorial(n)*factorial(r)
    
    for sig in DominantSigs(r,n,d):
        ps=1
        for s in sig:
            ps*=factorial(s)
        pcs=1
        for cs in Counter(sig).values():
            pcs*=factorial(cs)
        gamma+=factorialnr//(ps*pcs)
    return gamma

# Core chromarities -------------------------------------

def Kck(r:int,n:int,d:int)->float:
    """Kappa core chromarity.

    Args:
        r (int): number of colors
        n (int): community size
        d (int): size of nodes with the same color

    Returns:
        float: chromarity value
    """
    return d/n*(1-Kappa(r,n,d)/r**n)

# Gamma Core Chromarity
def Kcg(r:int,n:int,d:int)->float:
       
       """Gamma core chromarity.

    Args:
        r (int): number of colors
        n (int): community size
        d (int): size of the dominant color nodes

    Returns:
        float: chromarity value
    """
       return d/n*(1-Gamma(r,n,d)/r**n)
    
# Chromarities ------------------------------------------

#Kappa Chromarity
def Kk(P:set,c:dict,r:int)->float:
    """"Compute the Kappa chromarity.

    Args:
        P (set): community structure
        c (dict): color profile
        r (int): number of colors

    Returns:
        float: chromarity value.
    """
    if P:   
        k=0
        for p in P:
            colors=CleanCommunityColorProfile(p,c).values()
            if colors:
                k+=Kck(r,len(p),max(Counter(colors).values()))
        return k/len(P)
    else:
        return 0

# Gamma Chromarity
def Kg(P:set,c:dict,r:int)->float:
    """"Compute the Gamma chromarity.

    Args:
        P (set): community structure
        c (dict): color profile
        r (int): number of colors

    Returns:
        float: chromarity value.
    """
    if P:   
        k=0
        for p in P:
            colors=CleanCommunityColorProfile(p,c).values()
            if colors:
                k+=Kcg(r,len(p),max(Counter(colors).values()))
        return k/len(P)
    else:
        return 0

# GRAPH =============================================================================
# Graph ----------------------------------------------------------

#Default palette
__palette__={0:'lightgray', 1:'crimson', 2:'steelblue', 3:'gold', 4:'lightgreen',5:'mediumpurple',6:'darkorange',7:'burlywood'}

#Default font
__font__='Franklin Gothic Heavy'  #other nice fonts  'Tahoma'  'Impact'

def DrawColoredGraph(G, palette=__palette__,pos=None):
    """Display a colored graph

    Args:
        G (Graph): colored graph
        palette (dict, optional): color palette associating a color to integer. Defaults to __palette__ (6 colors).
        pos (dict|None, optional): node position. Defaults to None.
    """
    color=nx.get_node_attributes(G,"color")
    nx.draw_networkx(G,with_labels=True, pos=pos, node_color=[palette[c] for c in color.values()],  node_size=500, \
                      font_size=11, font_color='black', font_family=__font__)

# Display the community structure on graph 
def DrawChroCoS(G,P, theme='Set2', pos=None):
    """ Display the chromarity structure on a graph. The nodes of the same community are the same color.
    Args:
        G (Graph): _Undirected colored graph_
        P (set): _community structure_
        theme (str, optional): theme color of seaborn package. Defaults to 'Set2'.
        pos (dict|None, optional): position of nodes. Defaults to None.
    """
  
    Pl=list(P)
    palette=sns.color_palette(theme, len(Pl))
    color= {v:palette[i] for i in range(len(Pl)) for v in Pl[i]}
    nx.draw_networkx(G,with_labels=True, pos=pos, node_color=[color[v] for v in G.nodes()],  node_size=500, \
                      font_size=11, font_color='black', font_family=__font__)
    
# Random Graph  --------------------------------------------------

def RandomColoring(G,seeds:list,density:float=0.2,transparency:float=0.):
    """Attributes colors to nodes of graph G randomly.

    Args:
        G (Graph): unidrected graph that must be a single component.
        seeds (list): list of seeds (nodes)
        density (float, optional): probability parameter higher the value less the colors are scattered . Defaults to 0.2.
        transparency (float, optional): probability of transparent nodes. Defaults to 0..
    """
 
     # sub function selecting the color from seeds w.r.t. an exponential probability law. The farther from a color seed the node, the lower the probability of the color is.
    def ChooseColorRandomly(seeds,v):
        return choices(range(1,len(seeds)+1),weights=[exp(-density*nx.shortest_path_length(G,seed,v)) for seed in seeds],k=1)[0]

    assert 0.<= density <=1.0
    assert 0.<= transparency <=1.0
    assert nx.is_connected(G) #require that a path exists for any pair of vertices.
    
    nx.set_node_attributes(G, dict([(v,ChooseColorRandomly(seeds,v)) for v in G.nodes()]),"color")
    
    if transparency >0:
        transparent=[v for v in G.nodes() if random()< transparency]
        nx.set_node_attribute(G,dict.fromkey(transparent,0),"color")

#Generate a list of r seeds.      
def GenerateSeeds(G,r:int):
    """Generate r color seeds for graph G by maximizing the  geometric mean distance between them.

    Args:
        G (Graph): undirected graph that must be a single component.
        r (int): number of colors with r>1
    Returns: 
    list: color seeds
    """
    assert(r>1)
    assert nx.is_connected(G) #require that a path exists for any pair of vertices.
    pathlength=dict(nx.all_pairs_shortest_path_length(G))

    # first select two nodes with a maximal distance
    maxi=0
    v=w=next(iter(pathlength))
    for v,vpathlength in pathlength.items():
        for w,distance in vpathlength.items():
            if distance>maxi:
                maxi=distance
                vmax=v
                wmax=w
    seeds=[vmax,wmax]
    
    # Next complete the seed list with nodes maximizing the mean geometric distance.
    V=list(G.nodes())
    V.remove(vmax)
    V.remove(wmax) 
    for step in range(r-2):
        maxi=0
        vmax=None
        for v in V:
            disttoseeds=[pathlength[v][seed] for seed in seeds if seed in pathlength[v]]
            if disttoseeds:
                meandisttoseeds = gmean(disttoseeds)
                if meandisttoseeds > maxi:
                    maxi=meandisttoseeds
                    vmax=v
        seeds.append(vmax)
        V.remove(vmax)
    return(seeds)
            
# CHROCODE ===========================================================================

def MonochromeCommunityStructure(G):
    """Compute a monochrome community structure of graph G.

    Args:
        G (Graph): Undirected colored graph 

    Returns:
        set of frozensets : community structure
    """
    V=set(G.nodes())
    P=set()
    while V:
        v= next(iter(V))                    # take a node 
        referencecolor=G.nodes[v]['color']  # get its color used as reference colore
        p=set()
        monochromeneighbors={v}
        while monochromeneighbors: 
            v=monochromeneighbors.pop()     # take a node for including its neighbors 
            p.add(v)                        # add the current node to the community 
            for w in G[v]:                  # extend the neighbor set with the same color
                if G.nodes[w]['color']==referencecolor and not w in p:
                    monochromeneighbors.add(w)
                    
        V = V - p                           # remove the community of the examined vertices. 
        P.add(frozenset(p))                 # add the community to the structure
    return P


def ChroCoDe(G,r:int,radius:int=2,K=Kg):
    """Find a chromatic community structure.
    Args:
        G (Graph): Colored undirected graph 
        r (int): number of colors
        radius (int, optional): neigborhood distance. Defaults to 2.
        K (function, optional): Chromarity. Defaults to Kg.

    Returns:
        set of frozensets : community structure.
    """
    assert(radius>=1)
    assert(r>0)
    
    colorprofile=nx.get_node_attributes(G,'color')
    QG=nx.quotient_graph(G,MonochromeCommunityStructure(G)) # Quotient graph of the community structure
    P=set(QG.nodes())
    Pscan=P.copy()
    
    while Pscan:
        kmin=inf
        p=set()
        for q in Pscan:                                     # find the community p in P with the minimal core chromarity
            k=K({q},colorprofile,r)
            if k<kmin:
                k=kmin
                p=q
        Pscan.remove(p)

        N= list(nx.ego_graph(QG,p,radius=radius).nodes())   # compute the neighborhood of size dist but p
        N.remove(p)
        
        kmax=K(P,colorprofile,r)
        improved=False
        for q in N:                                         # find the neighbor q of p maximizing K
            communitypath=set(nx.shortest_path(QG,p,q))
            pmerge = reduce(lambda p,q: p|q,communitypath)

            k=K((P - communitypath) | {pmerge},colorprofile,r)
            if kmax<k:
                improved=True
                kmax=k
                maxpath=communitypath.copy()
                
        if improved:                                        # Update P with the shortest path connecting p to this neighbor 
            pmerge = reduce(lambda p,q: p|q,maxpath)
            P = (P - maxpath) | {pmerge}
            QG=nx.quotient_graph(G,P)                       # Update the quotient graph QG and the running community set Pscan
            Pscan=P.copy()
    return P