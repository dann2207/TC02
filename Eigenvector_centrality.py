# Eigenvector:
import networkx as nx
import matplotlib.pyplot as plt
G_APMS = nx.read_edgelist('yeast_AP-MS.txt')
G = G_APMS.copy()
cen = nx.eigenvector_centrality(G)
# cen[x] = centralidad del nodo x (un numerito)
# Ahora, vamos a ordenar los nodos de mayor a menor centralidad
# Nos armamos, a partir de cen, una lista de tuplas
# cen_sort=[(n0,c0),(n1,c1),(n2,c2),(n3,c3),(n4,c4),...]
# con c0 > c1 > c2 > c3 > c4 > ...
cen_sort = sorted(cen.items(),key=lambda x:x[1],reverse = True)

# (X,Y) = lo que se va a plotear
X = []
Y = []
N = len(G.nodes())
for i in range(1,N):
    x = i/N
    G.remove_node(cen_sort[i][0])
    gig = max(nx.connected_component_subgraphs(G), key=len)
    y = len(gig)/len(G_APMS)
    X.append(x)
    Y.append(y)

X_eigenvector = X
Y_eigenvector = Y

plt.plot(X,Y)
plt.show()
