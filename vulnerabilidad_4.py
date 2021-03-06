# Inciso (C)
# tabla 3 y figura 3 (Zotenko)
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np

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

# Armé una lista que es los datasets (DS) y otra que es los nombres 
# de los datasets (DSl), donde la l es por labels

DS=[LIT,Y2H,APMS,e_LITR]
DSl=["LIT","Y2H","APMS","LITR"]


G_APMS=nx.Graph()
G_APMS.add_edges_from(DS[2])
G_LIT =nx.Graph()
G_LIT.add_edges_from(DS[0])
G_Y2H =nx.Graph()
G_Y2H.add_edges_from(DS[1])
G_LITR =nx.Graph()
G_e_LITR.add_edges_from(DS[3])


Grafos=dict{}
Grafos["G_LIT"]=G_LIT
Grafos["G_APMS"]=G_APMS
Grafos["G_Y2H"]=G_Y2H
Grafos["G_LITR"]=G_e_LITR


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
G = G_APMS.copy()
ess_APMS_set = set(ess_list).intersection(set(G.nodes()))
# Sacamos los nodos esenciales:
for n in ess_APMS_set:
    G.remove_node(n)

# Ahora, nos fijamos en la componente gigante superviviente
giant = max(nx.connected_component_subgraphs(G), key=len)

# Finalmente, calculamos la proporcion
# (Seria la columna "Essential" de la tabla 3)
r_APMS = len(giant)/len(G_APMS)

# Volvemos a definir G y aca no paso nada
G = G_APMS.copy()

# Ahora, hay que hacer la otra columna...
# Recordemos que ess_APMS_set es el conjunto de
# nodos esenciales de la red. A partir de este
# conjunto, calculamos la distribución de grados:
from funciones import kdist
kdist_APMS = kdist(ess_APMS_set,G_APMS)
kd = dict(kdist_APMS)

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
# Pero, al ir avanzando hacia otros grados cercanos, puede
# que nos topemos con alguno de los k

# Para empezar, sea K el conjunto de grados de 
# los nodos esenciales
K = {x[0] for x in kdist_APMS}

# Sea K_buenos los k tales que haya suficientes
# (nk o mas) nodos no esenciales con grado k,
# y K_malos los k restantes
from funciones import esbueno
K_buenos = {k for k in K if esbueno(k,ess_APMS_set,G) == True}
K_malos = K.difference(K_buenos)

# Con los K_buenos no hay problema
# Pero si partimos de un k_malo, tenemos que ir hacia k's
# cercanos, acumulando sus nodos, hasta tener nk_malo (o mas)
# nodos. Pero al ir hacia otros k's, nos podemos topar
# con un k_bueno o un k_malo
# Tenemos que decidir que hacer si se da alguno de esos
# dos casos (con criterio...)

# Nuestro k_malo es k1. Si nos topamos con otro k_malo,
# digamos k2, ahora vamos a tener que buscar nk2 nodos mas,
# o sea, ibamos buscando nk1 nodos, no los encontramos,
# y ahora que nos topamos con k2 tenemos que encontrar
# nk1 + nk2 nodos en total
# Si, antes de llegar a los nk1 + nk2 nodos, nos volvemos
# a topar con un k_malo = k3, vamos a tener que seguir
# avanzando, hasta llegar a tener nk1 + nk2 + nk3 nodos
# (como minimo). Y asi siguiendo...

# Lo anterior en algun momento se termina (creo)

# Que hacemos si, antes de juntar los nodos, llegamos
# a un k_bueno?
# Propongo que si eso pasa, saquemos los nodos que faltan
# de los esenciales (con grados (K_malos) k1, k2, etc ),
# "al azar". Voy a pensar en como hacer esto solo si en
# el codigo vemos que sale mal y hay que hacerlo si o si

# Armamos un diccionario
# en las que las keys son los grados k
# y los values son los k's mas proximos por debajo
from funciones import ktup2
k_k = ktup2(K)

# F[k] = lista de nodos NO esenciales con grado k
from funciones import nod_k
K_todos = {G.degree(n) for n in G.nodes()}
nodk = nod_k(G)
F = {k: list(set(nodk[k]).difference(ess_APMS_set)) for k in K_todos}

k_max = max(G.degree(n) for n in G.nodes())
for k in range(k_max):
    if k not in F.keys():
        F[k]=[]

H = {}

# Completemos a H con los k_buenos
for k in K_buenos:
    H[k] = F[k]


for k in K_malos:
    nk = kd[k]
    N = [] # lista de nodos acumulada
    for c in range(len(F[k])):
        N.append(F[k][c])
    j = 1
    while len(N)<nk and k-j>0:
        if k-j in K_buenos and len(H[k-j])>kd[k-j]:
            N.append(H[k-j][0])
            del(H[k-j][0])
            a=nk-len(N)
            #print("Faltan {} elementos".format(a))
            #print("Nodos NE de grado {}".format(k))
            #print("Estoy tomando nodos de grado {}".format(k-j))
            #print("Puedo tomar {} nodos mas".format(len(H[k-j])-kd[k-j]))
        else:
            j+=1
    #¿Que hacemos si len(N) no es mayor que nk?
    # Pidámole prestado nodos a las listas del dict anterior
    # "Recemos" para que el while haya terminado
    # con len(N)>nk
    H[k] = N
# Las keys de H son los k_malos
# los values son las listas (si todo salio bien, ya que nuestro 
# "rezo" pudo no haber sido escuchado) de nodos no esenciales
# con "grado parecido", de las que vamos a sacar
# nk nodos al azar varias veces (o todas las veces posible)
# De ahi va a salir un valor medio y una dispersion en
# ese ratio r



for k in H.keys():
    if len(H[k])<kd[k]:
        print("Para el grado {} te faltan nodos".format(k))


# Me parece que seria demasiado pesado si optamos por
# sacar de todas las formas posibles
# Saquemos unos 200 (ponele)


# 1 ensayo consiste en sacar nk nodos para los 
# k_esenciales, y medir r
# Queremos hacer M ensayos
# y guardar los r que van resutando en una lista R
# Despues, sobre esa lista R vamos a poder sacar 
# el valor medio y la varianza (por fin)
M = 20
R = []
# Me falta revisar mejor lo siguiente
for i in range(M):
    l = []
    for k in K:
        nk = kd[k]
        rs = random.sample(H[k],nk)
        for x in rs:
            l.append(x)
    for n in l:
        G.remove_node(n)
    gig = max(nx.connected_component_subgraphs(G), key=len)
    r = len(gig)/len(G_APMS)
    R.append(r)
    G = G_APMS.copy()

# SIIII salio bien!







# Shortest-path:
SP = {}
for x in G_APMS.nodes():
    G = G_APMS.remove_node(x)
    for a in G.nodes():
        for b in G.nodes():
            shortpaths = nx.all_shortest_paths(G_APMS,a,b)
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
