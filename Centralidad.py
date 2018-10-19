import random
import copy
import time
import csv
import math
import networkx as nx
import matplotlib.pyplot as plt

G_APMS = nx.read_edgelist('yeast_AP-MS.txt')
G_LIT = nx.read_edgelist('yeast_LIT.txt')
G_Y2H = nx.read_edgelist('yeast_Y2H.txt')

from funciones import ldata
LITR = ldata("yeast_LIT_Reguly.txt")
e_LITR = [[x[0],x[1]]  for x in LITR]
G_LITR = nx.Graph()
G_LITR.add_edges_from(e_LITR)

GR = [G_APMS,G_LIT,G_Y2H,G_LITR]

Nombre = {G_APMS: "G_APMS",
          G_LITR: "G_LITR",
          G_LIT: "G_LIT",
          G_Y2H: "G_Y2H"}

DICT = {}


for gr in GR:
    #Eigenvector------------------------------------------------------
    G = gr.copy()
    cen = nx.eigenvector_centrality(G)
    cen_sort = sorted(cen.items(),key=lambda x:x[1],reverse = True)
    # (X,Y) = lo que se va a plotear
    X = []
    Y = []
    N = len(G.nodes())
    T = len(max(nx.connected_component_subgraphs(G), key=len))
    for i in range(1,N):
        x = i/N
        G.remove_node(cen_sort[i][0])
        gig = max(nx.connected_component_subgraphs(G), key=len)
        y = len(gig)/T
        X.append(x)
        Y.append(y)
    X_eigenvector = X
    Y_eigenvector = Y
    #Betweenness------------------------------------------------------
    G = gr.copy()
    cen = nx.betweenness_centrality(G)
    cen_sort = sorted(cen.items(),key=lambda x:x[1],reverse = True)
    # (X,Y) = lo que se va a plotear
    X = []
    Y = []
    N = len(G.nodes())
    T = len(max(nx.connected_component_subgraphs(G), key=len))
    for i in range(1,N):
        x = i/N
        G.remove_node(cen_sort[i][0])
        gig = max(nx.connected_component_subgraphs(G), key=len)
        y = len(gig)/T
        X.append(x)
        Y.append(y)
    X_betweenness = X
    Y_betweenness = Y
    # Random------------------------------------------------------
    G = gr.copy()
    nodes = set(G.nodes())
    N = len(nodes)
    X = [0]*(N-1)
    Y = [0]*(N-1)
    M = 5
    t0 = time.time()
    T = len(max(nx.connected_component_subgraphs(G), key=len))
    for j in range(M):
        G = gr.copy()
        nod = copy.deepcopy(nodes)
        for i in range(1,N):
            n = random.choice(list(nod))
            nod.remove(n)
            x = i/N
            G.remove_node(n)
            gig = max(nx.connected_component_subgraphs(G), key=len)
            y = len(gig)/T
            X[i-1] += x/M
            Y[i-1] += y/M

    t = time.time()
    print(t-t0)
    X_random = X
    Y_random = Y
    #Degree------------------------------------------------------
    G = gr.copy()
    nk = list(G.degree())
    nk.sort(key = lambda x: x[1], reverse = True)
    t0 = time.time()
    T = len(max(nx.connected_component_subgraphs(gr), key=len))
    X,Y = [],[]
    N = len(G)
    for i in range(1,N):
        x = i/N
        G.remove_node(nk[i][0])
        g = max(nx.connected_component_subgraphs(G), key=len)
        y = len(g)/T
        X.append(x)
        Y.append(y)
    t = time.time()
    print(t-t0)
    X_degree = X
    Y_degree = Y
    #DICT------------------------------------------------------
    DICT[gr] = {}
    DICT[gr]["X_degree"] = X_degree
    DICT[gr]["X_random"] = X_random
    DICT[gr]["X_betweenness"] = X_betweenness
    DICT[gr]["X_eigenvector"] = X_eigenvector
    DICT[gr]["Y_degree"] = Y_degree
    DICT[gr]["Y_random"] = Y_random
    DICT[gr]["Y_betweenness"] = Y_betweenness
    DICT[gr]["Y_eigenvector"] = Y_eigenvector

# DICT es un diccionario cuyas keys son los grafos
# y los values son diccionarios {"x": x}, en donde x
# son las listas que se usan para plotear la figura 3
for gr in DICT.keys():
	for s in DICT[gr].keys():
		