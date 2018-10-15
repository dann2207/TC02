
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
# y el segundo elemento es la cantidad de nodos con ese grado
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
