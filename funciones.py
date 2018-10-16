
# Funcion que lee un archivo de datos y devuelve una
# matriz que contiene estos datos
def ldata(archive):
        f = open(archive)
        data = []
        for line in f:
            col = line.split("\t")
            col = [x.strip() for x in col]
            data.append(col)
        return data 



# La siguiente funcion calcula la distribucion de grados 
# agarrando una lista de nodos y un grafo,
# y devuelve una lista de tuplas de tipo int
# en las que el primer elemento de cada tupla es el grado,
# y el segundo elemento es la cantidad de nodos de la lista con ese grado
# Hace falta especificar G ya que hay proteinas
# repetidas en varias redes distintas
def kdist(node_list,G):
    D = {}
    for n in node_list:
        k = G.degree(n)
        if k not in D.keys():
            D[k] = 1
        else:
            D[k] += 1
    tuple_list = sorted(D.items(),key=lambda x:x[0])
    return tuple_list



# ktup[i] = (k_i, k mas proximo (puede ser una lista de 2 elem))
# K puede ser lista o set
def ktup(K):
    K = list(K)
    tup = []
    for k in K:
        l = []
        d = min([np.abs(k-j) for j in K if j!=k])
        for j in K:
            if np.abs(k-j) == d:
                l.append(j)
        tup.append((k,l))
    return tup

def ktup2(k):
    k = sorted(k)
    k.insert(0,0) # mete un cero a la izq
    l = []
    for i in range(len(k)-1):
        l.append(k[i])
    return dict(zip(k[1:],l))

# nod_k es un dicionario 
# "k: lista de nodos con grado k"
# Por alguna razon no funciono lo de abajo
#def nod_k(G):
    #D = {}
    #for n in G.nodes():
        #if G.degree(n) not in D.keys():
            #k = G.degree(n)
            #D[k] = [n]
        #else:
            #D[k].append(n)
    #return D

def nod_k(G):
    D = {}
    K = {G.degree(n) for n in G.nodes()}
    for k in K:
        D[k] = [x for x in G.nodes() if k==G.degree(x)]
    return D

#
def esbueno(k,ess_list,G):
    a = len([x for x in ess_list if k==G.degree(x)])
    # a = cantidad de nodos esenciales con grado k
    b = len(nod_k(G)[k])
    # b = cantidad de nodos con grado k
    if 2*a <= b:
        return True
    else:
        return False

