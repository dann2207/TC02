import random
import copy
import time
import math
import networkx as nx
import matplotlib.pyplot as plt

G_APMS = nx.read_edgelist('yeast_AP-MS.txt')
G = G_APMS.copy()

nk = list(G.degree())
nk.sort(key = lambda x: x[1], reverse = True)

t0 = time.time()
T = len(max(nx.connected_component_subgraphs(G_APMS), key=len))
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