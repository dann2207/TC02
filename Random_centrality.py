# Random
import random
import copy
import time
import networkx as nx
import matplotlib.pyplot as plt

G_APMS = nx.read_edgelist('yeast_AP-MS.txt')
G = G_APMS.copy()

nodes = set(G.nodes())
N = len(nodes)

X = [0]*(N-1)
Y = [0]*(N-1)


M = 1
t0 = time.time()
for j in range(M):
    G = G_APMS.copy()
    nod = copy.deepcopy(nodes)
    for i in range(1,N):
        n = random.choice(list(nod))
        nod.remove(n)
        x = i/N
        G.remove_node(n)
        gig = max(nx.connected_component_subgraphs(G), key=len)
        y = len(gig)/len(G_APMS)
        X[i-1] += x/M
        Y[i-1] += y/M

t = time.time()
print(t-t0)

plt.plot(X,Y)
plt.show()