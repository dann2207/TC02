# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 15:47:23 2018

@author: Fabio
"""

import networkx as nx
import math
import matplotlib.pylab as plt
import numpy as np

def ldata(archive):
        f=open(archive)
        data=[]
        for line in f:
            line=line.strip()
            col=line.split("\t") # Cambiamos las separaciones para que sean con tabs
            data.append(col)
        return data

    
# Armo las tablas en el programa
# Al ldata le cambiamos para que las separaciones sean por tabs y no por espacios
    
LIT=ldata("yeast_LIT.txt")
APMS=ldata("yeast_AP-MS.txt")
Y2H=ldata("yeast_Y2H.txt")
LITR=ldata("yeast_LIT_Reguly.txt")
PaHe=ldata("Essential_ORFs_paperHe.txt")


# Con esto armo una lista de enlaces para la lista del LIT Reguly
# Porque el archivo LIT Reguly pasado a lista tenía más de dos columnas 
# y no se lo podía pasar al comando G.add_edges


e_LITR=[[x[0],x[1]]  for x in LITR]
del(e_LITR[0])


# Con el siguiente comando armo la lista de proteínas esenciales

del(PaHe[len(PaHe)-1])
del(PaHe[len(PaHe)-1])
del(PaHe[len(PaHe)-1])
del(PaHe[len(PaHe)-1])

Ess=[x[1] for x in PaHe]
del(Ess[0])
del(Ess[0])

# Armé una lista que es los datasets (DS) y otra que es los nombres 
# de los datasets (DSl), donde la l es por labels

DS=[LIT,Y2H,APMS,e_LITR]
DSl=["LIT","Y2H","APMS","e_LITR"]
#
#
#for i in range(len(DS)):
#    
#    # Armo el grafo para cada dataset
#
#    G=nx.Graph()
#    G.add_edges_from(DS[i])
#    
#    
#    # Datos para la tabla 1
#    
#    # Numero de nodos
#    N=len(G.nodes())
#    print("El numero de nodos para la red {} es:{}".format(DSl[i],N))
#    
#    # Numero de enlaces
#    E=len(G.edges())
#    print("El numero de enlaces para la red {} es:{}".format(DSl[i],E))
#    
#    # Grado medio
#    Grad= [x[1]  for x in G.degree()]
#    Gm=np.mean(Grad)
#    print("El valor del grado medio para la red {} es:{}".format(DSl[i],Gm))
#    
#    # Grado de Clustering medio
#    Cl=nx.average_clustering(G)
#    print("El coeficiente de Clustering promedio para la red {} es:{}".format(DSl[i],Cl))
#
#
## Tabla 2 Zotenko
#
#l=0
#
## Los siguientes tres for y el if son para hacer todas las combinatorias necesarias
## para armar la tabla 1 de Zotenko, donde cada dataset se compara con los otros tres
#
#for i in range(len(DS)):
#    for j in range(len(DS)):
#        if i!=j:
#            G=nx.Graph()
#            G.add_edges_from(DS[i])
#
#            L=nx.Graph()
#            L.add_edges_from(DS[j])
#            
#            # Armo las listas temporales y las curo
#            
#            temp1=[x for x in G.edges()]
#       
#            temp2=[x for x in L.edges()]         
#            l=l+1
#            
#            # Calculo del solapamiento
#            fc=0
#            
#            # el elemento tenlace1 lo uso para asegurarme que si un enlace está
#            # escrito al revés en una lista respecto de la otra, lo cuente igual
#            for enlace1 in temp1:
#                tenlace1=[enlace1[1],enlace1[0]]
#                for enlace2 in temp2:
#                    if enlace1==enlace2:
#                        fc=fc+(1/len(G.edges()))
#                    if tenlace1==enlace2:
#                        fc=fc+(1/len(G.edges()))
#                        
#                if enlace1==temp1[len(temp1)-1]:
#                    print("Cálculo número:{}".format(l))
#                    print("La fracción de enlaces comunes entre {} y {} es: {}".format(DSl[i],DSl[j],fc))
#                    

# Figura 1.a

G=nx.Graph()
G.add_edges_from(DS[0])

# Armo la lista de nodos con el grado asociado, la lsita NG (Nodos y Grados)

NG=[[x[0],x[1]] for x in G.degree()]

# Armo el criterio para que al ordenar, ordeno la lista NG segun grados

def tomarsegundo (elem):
    return elem[1]

NG.sort(key=tomarsegundo)


fr=np.arange(0.05,1.05,0.05) # Porcentajes de Cutoff del Hub


#for i in fr:
gradmin=math.floor((NG[len(NG)-1][1])*0.1)
HUBS=NG
print(gradmin)
for grado in HUBS:
    if grado[1] <= gradmin:
        print(grado)
        del(HUBS[0])
print(gradmin)
print(HUBS)