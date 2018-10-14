# Inciso (C)
# Tabla 3 Zotenko
import networkx as nx
import matplotlib.pyplot as plt
import random

G_APMS = nx.read_edgelist('yeast_AP-MS.txt')
G_LIT  = nx.read_edgelist('yeast_LIT.txt')
G_Y2H  = nx.read_edgelist('yeast_Y2H.txt')
#G_LITR = nx.read_edgelist('yeast_LIT_Reguly.txt')

N_APMS = G_APMS.nodes()
N_LIT  = G_LIT.nodes()
N_Y2H  = G_Y2H.nodes()
#N_LITR = G_LITR.nodes()

#--------------------------------------------------------
# Ahora veamos si algun par de redes comparte proteinas
# Tomamos, por ejemplo, las redes AP-MS y Y2H

#N_set_APMS = {n for n in N_APMS}
#N_set_Y2H = {n for n in N_Y2H}

#N_set_APMS.intersection(N_set_Y2H)

# Con esto vimos que hay muchas proteinas en comun
#--------------------------------------------------------

# CAMBIAR EL CODIGO: ES IMPORTANTE QUE LOS NODOS "AL AZAR" QUE SACAMOS
# SEAN NO-ESENCIALES. EN EL PAPER DE ZOTENKO SE SACARON TODOS LOS NODOS
# ESENCIALES. LO QUE ESTAMOS HACIENDO ACA ES NO SACAR TODOS, SINO ALGUNOS,
# Y NOS QUEDARIA UNA TABLA EN DONDE EN VEZ DE NUMEROS TENDRIAMOS GRAFICOS.
# ¿De Donde?
# HAY UNA COMPLICACION CON CALCULAR EL DESVIO ESTANDAR DE ESTOS GRAFICOS,
# YA QUE NO HAY UNA UNICA FORMA DE ELEGIR ALGUNOS NODOS ESENCIALES AL AZAR.

# EN ESTE CODIGO VAMOS A SACAR TODOS LOS NODOS ESENCIALES

from funciones import ldata
PaHe = ldata("Essential_ORFs_paperHe.txt")

# La lista de los 1160 nodos esenciales es la siguiente
ess_list = [x[1] for x in PaHe if len(x)>1]
# Por ejemplo, un elemento de esta lista es:
# PaHe[40][1] = 'YBR049C'

# Vamos a hacer las cosas en el siguiente orden:
# 1) AP-MS
# 2) LIT
# 3) Y2H
# 4) LIT_Reguly


# ------------ Red APMS --------------------------------------------

# Nodos esenciales de esta red:
ess_APMS_set = set(ess_list).intersection(set(G_APMS.nodes()))
# Sacamos los nodos esenciales:
for n in ess_APMS_set:
	G_APMS.remove_node(n)

# Ahora, nos fijamos en la componente gigante superviviente
giant = max(nx.connected_component_subgraphs(G_APMS), key=len)

# Finalmente, calculamos la proporcion
# (Seria la columna "Essential" de la tabla 3)
r_APMS = len(giant)/len(G_APMS)

# Ahora, hay que hacer la otra columna...
# Recordemos que ess_APMS_set es el conjunto de
# nodos esenciales de la red. A partir de este
# conjunto, calculamos la distribución de grados:
from funciones import kdist
kdist_APMS = kdist(ess_APMS_set,G_APMS)

# Ahora, pasamos a lo mas dificil de todo este codigo,
# que es el criterio

# En kdist_APMS hay tuplas (k,nk)
# De ahi sacamos los grados k de los nodos esenciales
# Asi como tenemos la red, tenemos que buscar, para
# cada k, otros nk nodos con grado k (ninguno de esos
# otros puede ser un nodo esencial)
# Habra ciertos k para los que esto se pueda hacer,
# pero puede pasar que para algun k no haya tantos
# nodos no esenciales con ese grado. Si eso pasa tenemos
# que buscar nodos con grados cercanos, y sacar esos nodos.
# Pero


















# Shortest-path:
SP = {}
for x in G_APMS.nodes():
	G = G_APMS.remove_node(x)
	for a in G.nodes():
		for b in G.nodes():
			shortpaths = nx.all_shortest_paths(G,a,b)
			for path in shortpaths:
				if x in path:
					c += 1
		G.remove_node(a)
	SP[x] = c
# Hasta aca, D[x] = cantidad de shortest-paths
# que pasan por el nodo x

SP_sort = sorted(SP.items(),key=lambda x:x[1])


# Eigenvector:
import scipy.sparse.linalg
w,v = scipy.sparse.linalg.eigs()


# ------------ Red Y2H -----------------
# Nodos esenciales de esta red:
ess_Y2H = set(ess_list).intersection(set(G_Y2H.nodes()))
for n in ess_Y2H:
	G_Y2H.remove_node(n)

Gc = max(nx.connected_component_subgraphs(G_Y2H), key=len)
r_Y2H = len(Gc)/len(G_Y2H)


# ------------ Red LIT -----------------
# Nodos esenciales de esta red:
ess_LIT = set(ess_list).intersection(set(G_LIT.nodes()))
for n in ess_LIT:
	G_LIT.remove_node(n)

Gc = max(nx.connected_component_subgraphs(G_LIT), key=len)
r_LIT = len(Gc)/len(G_LIT)

