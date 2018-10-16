import networkx as nx
import numpy as np
import matplotlib.pylab as plt


def ldata(archive):
        f = open(archive)
        data = []
        for line in f:
            col = line.split("\t")
            col = [x.strip() for x in col]
            data.append(col)
        return data 
    
# Armo las tablas en el programa
#
    
LIT=ldata("yeast_LIT.txt")
APMS=ldata("yeast_AP-MS.txt")
Y2H=ldata("yeast_Y2H.txt")
LITR=ldata("yeast_LIT_Reguly.txt")
PaHe=ldata("Essential_ORFs_paperHe.txt")
del(PaHe[:2])
del(PaHe[1156:])

e_LITR=[[x[0],x[1]]  for x in LITR]
##Armo el grafo, igual que antes, si es red LITR hay que usar  G.add_edges_from(e_LITR)
G=nx.Graph()
#G=nx.read_edgelist('yeast_Y2H.txt')
G.add_edges_from(e_LITR)


#Armo la lista de nodos esenciales
nodos_esenciales=[]

for i in range(len(PaHe)):
    nodos_esenciales.append(PaHe[i][1])
nodos_esenciales1=[]
for i in nodos_esenciales:
    if i not in nodos_esenciales1:
        nodos_esenciales1.append(i)
        
##Asigno la esencialidad o no a cada nodo, es más fácil asignar primero a todos como no esenciales y luego los esenciales
esenciales_LIT=set(nodos_esenciales1).intersection(set(G.nodes()))

for nodo in G.nodes():
    G.node[nodo]["Tipo"]="No-Esencial"

for nodo_esencial in esenciales_LIT:
    G.node[nodo_esencial]["Tipo"]="Esencial"

#armo mi lista de pares (nodo, grado)
lista_nodos_grados=list(G.degree())

#Armo una lista sólo con los grados y luego aplicando set() me quedo sólo con los que no se repiten
lista_grados=[]
for nodo_grado in lista_nodos_grados:
    lista_grados.append(nodo_grado[1])
lista_grados=set(lista_grados)


##Asigno a cada nodo su grado
for nodo_grado in lista_nodos_grados:
    G.node[nodo_grado[0]]["Grado"]=nodo_grado[1]
    
##Busco pares de nodos no adyacentes entre si que compartan 3 o más vecinos (para la lista y2h tomaré si comparten 1 o más)
pares_nodos_no_ady_entre_si=[]

for i,nodo in enumerate(G.nodes()):
    for j,nodo1 in enumerate(G.nodes()):
        if ((j>i) and (nodo1 not in G.neighbors(nodo)) and (len(set(G.neighbors(nodo)).intersection(set(G.neighbors(nodo1))))>=3)):
            pares_nodos_no_ady_entre_si.append([nodo, nodo1])
pares_nodos_no_ady_entre_si

#Cuento cuántos de esos pares son del mismo tipo, sean esencial o no esencial
pares_no_ady_esenciales=[]
pares_no_ady_no_esenciales=[]

for par in pares_nodos_no_ady_entre_si:
    if G.node[par[0]]["Tipo"]=="Esencial" and G.node[par[1]]["Tipo"]=="Esencial":
        pares_no_ady_esenciales.append(par)
        
for par in pares_nodos_no_ady_entre_si:
    if G.node[par[0]]["Tipo"]=="No-Esencial" and G.node[par[1]]["Tipo"]=="No-Esencial":
        pares_no_ady_no_esenciales.append(par)
        
##Muestro los pares bsucados
print("Número total de pares", len(pares_nodos_no_ady_entre_si))
print(len(pares_no_ady_esenciales))
print(len(pares_no_ady_no_esenciales))
print("número total de pares del mismo tipo", len(pares_no_ady_esenciales)+len(pares_no_ady_no_esenciales))

##Armo un nuevo grafo sin la característica de esencial o no en cada nodo, de nuevo, si la red es LITR lo armo a partir de la lista e_LIT
G1=nx.Graph()
#G1=nx.read_edgelist('yeast_Y2H.txt')
G1.add_edges_from(e_LITR)


#Datos obtenidos del ajuste del punto tanterior para valores de alfa y beta que usaré para las probabilidades de cada nodo de ser esencial según su grado

#i=0 - RED LIT; i=1 RED APMS, i=2 RED Y2H, i=3 RED LITR
alfa_redes=(0.07881, 0.06168, 0.01465, 0.04423)
beta_redes=(0.25495, 0.16844, 0.18216, 0.07375)

#Armo los arrays de probabilidades de ser esenciales o no, uso arrays para operar más fácil
proba_nodo_esencial_k1=[1-((1-alfa_redes[3])*(1-beta_redes[3])**k) for k in lista_grados]
proba_nodo_esencial_k1=np.array(proba_nodo_esencial_k1)
proba_nodo_no_esencial_k1=1-proba_nodo_esencial_k1




#Me creo copias de las listas usadas antes para el nuevo grafo G1
pares_nodos_no_ady_entre_si1=pares_nodos_no_ady_entre_si
lista_grados1=list(lista_grados)
proba_esencial_k1=[]
proba_no_esencial_k1=[]

#Asigno grado a cada nodo
for nodo_grado in lista_nodos_grados:
    G1.node[nodo_grado[0]]["Grado"]=nodo_grado[1]

#Armo listas de probabilidad según grado (probabilidad de ser esencial y de no serlo)
for i in range(len(lista_grados1)):
    proba_esencial_k1.append([lista_grados1[i], proba_nodo_esencial_k1[i]])

for i in range(len(lista_grados1)):
    proba_no_esencial_k1.append([lista_grados1[i], proba_nodo_no_esencial_k1[i]])
    
##Le asigno a cada nodo, según su grado, su probabilidad de ser esencial o no
for nodo in G1.nodes():
    for k_proba_e in proba_esencial_k1:
        if G1.node[nodo]["Grado"]==k_proba_e[0]:
            G1.node[nodo]["P_E"]=k_proba_e[1]

for nodo in G1.nodes():
    for k_proba_no_e in proba_no_esencial_k1:
        if G1.node[nodo]["Grado"]==k_proba_no_e[0]:
            G1.node[nodo]["P_NO_E"]=k_proba_no_e[1]


##Planteo la suma sobre todos los pares de nodos no adyacentes que compartan un número dado de vecinos en común y sean del mismo tipo
##usando las probabilidades antes asignadas
cant_nodos_ajuste=0

for par in pares_nodos_no_ady_entre_si1:
    cant_nodos_ajuste=cant_nodos_ajuste+(G1.node[par[0]]["P_E"]*G1.node[par[1]]["P_E"])+(G1.node[par[0]]["P_NO_E"]*G1.node[par[1]]["P_NO_E"])
print("Número de pares del mismo tipo según el ajuste lineal", cant_nodos_ajuste)

