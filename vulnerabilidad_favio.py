# Inciso (C)
# tabla 3 y figura 3 (Zotenko)
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np
import math


def ldata(archive):
        f = open(archive)
        data = []
        for line in f:
            col = line.split("\t")
            col = [x.strip() for x in col]
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

Ness=[x.split("-") for x in Ess]
Ness=[x[0] for x in Ness]

# Armé una lista que es los datasets (DS) y otra que es los nombres 
# de los datasets (DSl), donde la l es por labels

DS=[LIT,Y2H,APMS,e_LITR]
DSl=["LIT","Y2H","APMS","LITR"]


G_APMS = nx.Graph()
G_APMS.add_edges_from(DS[2])
G_LIT  = nx.Graph()
G_LIT.add_edges_from(DS[0])
G_Y2H  = nx.Graph()
G_Y2H.add_edges_from(DS[1])
G_LITR = nx.Graph()
G_LITR.add_edges_from(DS[3])

Gl=["LIT","Y2H","APMS","LITR"]
Grafos=[G_LIT,G_Y2H,G_APMS,G_LITR]
#SAS=[]
#Sis=[]

mar=["r","b","k","g"]


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

PaHe = ldata("Essential_ORFs_paperHe.txt")

# La lista de los 1160 nodos esenciales es la siguiente
ess_list = [x[1] for x in PaHe if len(x)>1]
# Por ejemplo, un elemento de esta lista es:
# PaHe[40][1] = 'YBR049C'

## Vamos a hacer las cosas en el siguiente orden:
## 1) AP-MS
## 2) LIT
## 3) Y2H
## 4) LIT_Reguly
#
#
## ------------ Red APMS --------------------------------------------
#
#for i in range(len(Grafos)):
#    # Nodos esenciales de esta red:
#    G = Grafos[i].copy()
#    ess_set = set(ess_list).intersection(set(G.nodes()))
#    # Sacamos los nodos esenciales:
#    for n in ess_set:
#        G.remove_node(n)
#
#    # Ahora, nos fijamos en la componente gigante superviviente
#    giant = max(nx.connected_component_subgraphs(G), key=len)
#    
#    # Finalmente, calculamos la proporcion
#    # (Seria la columna "Essential" de la tabla 3)
#    r_es = len(giant)/len(G)
#    Sis.append(r_es)
#    
#    
#    # Volvemos a definir G y aca no paso nada
#    G = Grafos[i].copy()
#
#    # Ahora, hay que hacer la otra columna...
#    # Recordemos que ess_APMS_set es el conjunto de
#    # nodos esenciales de la red. A partir de este
#    # conjunto, calculamos la distribución de grados:
#    from funciones import kdist
#    KDI = kdist(ess_set,G)
#    kd = dict(KDI)
#
#    # Ahora, pasamos a lo mas dificil de todo este codigo,
#    # que es el criterio
#    
#    # En kdist_APMS hay tuplas (k,nk)
#    # De ahi sacamos los grados k de los nodos esenciales
#    # Asi como tenemos la red, tenemos que buscar, para
#    # cada k, otros nk nodos con grado k (ninguno de esos
#    # otros puede ser un nodo esencial)
#    # Habra ciertos k para los que esto se pueda hacer,
#    # pero puede pasar que para algun k no haya tantos
#    # nodos no esenciales con ese grado. Si eso pasa tenemos
#    # que buscar nodos con grados cercanos, y sacar esos nodos.
#    # Pero, al ir avanzando hacia otros grados cercanos, puede
#    # que nos topemos con alguno de los k
#
#    # Para empezar, sea K el conjunto de grados de 
#    # los nodos esenciales
#    K = {x[0] for x in KDI}
#
#    # Sea K_buenos los k tales que haya suficientes
#    # (nk o mas) nodos no esenciales con grado k,
#    # y K_malos los k restantes
#    from funciones import esbueno
#    K_buenos = {k for k in K if esbueno(k,ess_set,G) == True}
#    K_malos = K.difference(K_buenos)
#
#    # Con los K_buenos no hay problema
#    # Pero si partimos de un k_malo, tenemos que ir hacia k's
#    # cercanos, acumulando sus nodos, hasta tener nk_malo (o mas)
#    # nodos. Pero al ir hacia otros k's, nos podemos topar
#    # con un k_bueno o un k_malo
#    # Tenemos que decidir que hacer si se da alguno de esos
#    # dos casos (con criterio...)
#    
#    # Nuestro k_malo es k1. Si nos topamos con otro k_malo,
#    # digamos k2, ahora vamos a tener que buscar nk2 nodos mas,
#    # o sea, ibamos buscando nk1 nodos, no los encontramos,
#    # y ahora que nos topamos con k2 tenemos que encontrar
#    # nk1 + nk2 nodos en total
#    # Si, antes de llegar a los nk1 + nk2 nodos, nos volvemos
#    # a topar con un k_malo = k3, vamos a tener que seguir
#    # avanzando, hasta llegar a tener nk1 + nk2 + nk3 nodos
#    # (como minimo). Y asi siguiendo...
#    
#    # Lo anterior en algun momento se termina (creo)
#    
#    # Que hacemos si, antes de juntar los nodos, llegamos
#    # a un k_bueno?
#    # Propongo que si eso pasa, saquemos los nodos que faltan
#    # de los esenciales (con grados (K_malos) k1, k2, etc ),
#    # "al azar". Voy a pensar en como hacer esto solo si en
#    # el codigo vemos que sale mal y hay que hacerlo si o si
#
#    # Armamos un diccionario
#    # en las que las keys son los grados k
#    # y los values son los k's mas proximos por debajo
#    from funciones import ktup2
#    k_k = ktup2(K)
#
#    # F[k] = lista de nodos NO esenciales con grado k
#    from funciones import nod_k
#    K_todos = {G.degree(n) for n in G.nodes()}
#    nodk = nod_k(G)
#    F = {k: list(set(nodk[k]).difference(ess_set)) for k in K_todos}
#
#
#    H = {}
#    
#    # Completemos a H con los k_buenos
#    for k in K_buenos:
#        H[k] = F[k]
#    
#    
#    for k in K_malos:
#        nk = kd[k]
#        N = [] # lista de nodos acumulada
#        for c in range(len(F[k])):
#            N.append(F[k][c])
#        j = 1
#        while len(N)<nk and k-j>0:
#            if k-j in K_buenos and len(H[k-j])>kd[k-j]:
#                N.append(H[k-j][0])
#                del(H[k-j][0])
##                a=nk-len(N)
##                print("Faltan {} elementos".format(a))
##                print("Nodos NE de grado {}".format(k))
##                print("Estoy tomando nodos de grado {}".format(k-j))
##                print("Puedo tomar {} nodos mas".format(len(H[k-j])-kd[k-j]))
#            else:
#                j+=1
#        #¿Que hacemos si len(N) no es mayor que nk?
#        # Pidámole prestado nodos a las listas del dict anterior
#        # "Recemos" para que el while haya terminado
#        # con len(N)>nk
#        H[k] = N
#    # Las keys de H son los k_malos
#    # los values son las listas (si todo salio bien, ya que nuestro 
#    # "rezo" pudo no haber sido escuchado) de nodos no esenciales
#    # con "grado parecido", de las que vamos a sacar
#    # nk nodos al azar varias veces (o todas las veces posible)
#    # De ahi va a salir un valor medio y una dispersion en
#    # ese ratio r
#
#    l=False
#    
#    for k in H.keys():
#        if len(H[k])<kd[k]:
#            print("Para el grado {} te faltan nodos".format(k))
#            l=True
#        if l==True:
#            print("Hay grados que no tienen los nodos suficientes")
#
#
#    # Me parece que seria demasiado pesado si optamos por
#    # sacar de todas las formas posibles
#    # Saquemos unos 200 (ponele)
#    
#    
#    # 1 ensayo consiste en sacar nk nodos para los 
#    # k_esenciales, y medir r
#    # Queremos hacer M ensayos
#    # y guardar los r que van resutando en una lista R
#    # Despues, sobre esa lista R vamos a poder sacar 
#    # el valor medio y la varianza (por fin)
#    M = 500
#    R = []
#    # Me falta revisar mejor lo siguiente
#    for u in range(M):
#        for k in K:
#            nk = kd[k]
#            rs = random.sample(H[k],nk)
#            for n in rs:
#                G.remove_node(n)
#        gig = max(nx.connected_component_subgraphs(G), key=len)
#        r = len(gig)/len(Grafos[i])
#        R.append(r)
#        G = Grafos[i].copy()
#        if len(R)==500:
#            fin=np.mean(R)
#            fan=np.std(R)
#            A=[fin,fan]
#            SAS.append(A)


# Ahora, nos fijamos en la componente gigante superviviente
#giant = max(nx.connected_component_subgraphs(G), key=len)

## Centralidad Betwenness
## Creo que funca de una
#p=0
#BEW=nx.betweenness_centrality(Grafos[p])
#
#CBW=[[x[0],x[1]] for x in BEW.items()]
#CBW.sort(key=lambda x:x[1],reverse=True)
#
#veces=math.floor(len(Grafos[p].nodes())*0.3)
#TO=len(max(nx.connected_component_subgraphs(Grafos[p]), key=len))
#
#Y=[]
#X=[]
#
#for v in range(veces):
#    Grafos[p].remove_node(CBW[v][0])
#    giant = max(nx.connected_component_subgraphs(Grafos[p]), key=len)
#    x=1-len(Grafos[p].nodes())/len(CBW)
#    y=len(giant)/TO
#    Y.append(y)
#    X.append(x)
#
#plt.plot(X,Y,"{}".format(mar[p]),label="{}".format(DSl[p]))
#
#plt.grid(True)
#plt.rcParams["figure.figsize"] = [10,5]
#plt.xlabel("Fraccion de nodos removidos")
#plt.ylabel("Fracción del tamaño de la componente gigante")
#plt.title("Impacto de la remoción de nodos en función de la centralidad")
#plt.legend()
#plt.show()

G_APMS = nx.Graph()
G_APMS.add_edges_from(DS[2])
G_LIT  = nx.Graph()
G_LIT.add_edges_from(DS[0])
G_Y2H  = nx.Graph()
G_Y2H.add_edges_from(DS[1])
G_LITR = nx.Graph()
G_LITR.add_edges_from(DS[3])

#Centralidad por grado
p=0

NG=[[x[0],x[1]] for x in Grafos[p].degree()]

def tomarsegundo (elem):
        return elem[1]

NG.sort(key=tomarsegundo,reverse=True)


veces=math.floor(len(Grafos[p].nodes())*0.3)
TO=len(max(nx.connected_component_subgraphs(Grafos[p]), key=len))

Y=[]
X=[]

for v in range(veces):
    Grafos[p].remove_node(NG[v][0])
    giant = max(nx.connected_component_subgraphs(Grafos[p]), key=len)
    x=1-len(Grafos[p].nodes())/len(NG)
    y=len(giant)/TO
    Y.append(y)
    X.append(x)

plt.plot(X,Y,"{}".format(mar[p]),label="{}".format(DSl[p]))

plt.grid(True)
plt.rcParams["figure.figsize"] = [10,5]
plt.xlabel("Fraccion de nodos removidos")
plt.ylabel("Fracción del tamaño de la componente gigante")
plt.title("Impacto de la remoción de nodos en función de la centralidad")
plt.legend()
plt.show()












#for i in range(len(Grafos)):
#    # Nodos esenciales de esta red:
#    G = Grafos[i].copy()
#    ess_set = set(ess_list).intersection(set(G.nodes()))


## Eigenvector:
#import scipy.sparse.linalg
#w,v = scipy.sparse.linalg.eigs()
#
#
## ------------ Red Y2H -----------------
## Nodos esenciales de esta red:
#ess_Y2H = set(ess_list).intersection(set(G_Y2H.nodes()))
#for n in ess_Y2H:
#    G_Y2H.remove_node(n)
#
#Gc = max(nx.connected_component_subgraphs(G_Y2H), key=len)
#r_Y2H = len(Gc)/len(G_Y2H)
#
#
## ------------ Red LIT -----------------
## Nodos esenciales de esta red:
#ess_LIT = set(ess_list).intersection(set(G_LIT.nodes()))
#for n in ess_LIT:
#    G_LIT.remove_node(n)
#
#Gc = max(nx.connected_component_subgraphs(G_LIT), key=len)
#r_LIT = len(Gc)/len(G_LIT)
#