import networkx as nx
import numpy as np
import matplotlib.pylab as plt
from scipy import optimize

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
DS=["yeast_LIT.txt", "yeast_AP-MS.txt", "yeast_Y2H.txt", "yeast_LIT_Reguly.txt"]

PaHe=ldata("Essential_ORFs_paperHe.txt")
##Borro lo que sobra del archivo
del(PaHe[:2])
del(PaHe[1156:])

##Armo la lista de LITR que tiene más columnas que las demás, para

e_LITR=[[x[0],x[1]]  for x in LITR]

##Armo una lista de los nombres de nodos esenciales
nodos_esenciales=[]

for i in range(len(PaHe)):
    nodos_esenciales.append(PaHe[i][1])
nodos_esenciales1=[]
for i in nodos_esenciales:
    if i not in nodos_esenciales1:
        nodos_esenciales1.append(i)
        

##Para generar el grafo de la red LITR uso G.add_edges_from(e_LITR), para los demás el comando de abajo

G=nx.Graph()
#G=nx.read_edgelist(DS[2]) ##Hay que comentar alguna de estas dos líneas para elegir el grafo a armar
G.add_edges_from(e_LITR)

print(G.number_of_selfloops())



    
    
##Asigno la esencialidad o no a cada nodo, es más fácil asignar primero a todos como no esenciales y luego los esenciales
esenciales_=set(nodos_esenciales1).intersection(set(G.nodes()))

for nodo in G.nodes():
    G.node[nodo]["Tipo"]="No-Esencial"

for nodo_esencial in esenciales_:
    G.node[nodo_esencial]["Tipo"]="Esencial"
    
        
#armo mi lista de pares (nodo, grado)
lista_nodos_grados=list(G.degree())

#Armo una lista sólo con los grados y luego aplicando set() me quedo sólo con los que no se repiten
for nodo_grado in lista_nodos_grados:
    G.node[nodo_grado[0]]["Grado"]=nodo_grado[1]
    
lista_grados=[]
for nodo_grado in lista_nodos_grados:
    lista_grados.append(nodo_grado[1])
lista_grados=set(lista_grados)

## Para cada valor de grado K me dijo cuantos enlaces esenciales hay
cuantos=[] ##Lista que tendrá la cantidad de nodos de grado k que son esenciales
for k in lista_grados:
    for nodo in G.nodes():
        if G.node[nodo]["Grado"]==k and G.node[nodo]["Tipo"]=="Esencial":
            cuantos.append(k)
            
lista_grados=list(lista_grados)#la vuelvo a transformar en lista



##Armo una lista de cuántos nodos esenciales hay dado su grado k, luego lo transformo en array para operar más fácil
cuantos_nodos_esenciales_k=[]
for k in lista_grados:
    cuantos_nodos_esenciales_k.append(cuantos.count(k))

cuantos_nodos_esenciales_k=np.array(cuantos_nodos_esenciales_k)


#Armo una lista de todos los nodos con grados en mi lista de grados
nodos_grado_k=[]
cuantos_nodos_grado_k=[]
for k in lista_grados:
    for nodo in G.nodes():
        if G.node[nodo]["Grado"]==k:
            nodos_grado_k.append(k)
            
##Cuento cuántos nodos en total hay para cada grado, luego lo transformo en array
for l in lista_grados:
    cuantos_nodos_grado_k.append(nodos_grado_k.count(l))

cuantos_nodos_grado_k=np.array(cuantos_nodos_grado_k)

##me saco de encima los valores de k>9 porque no son representativos estadisticamente, en la red Y2H me quedo con nodos hasta grado 10
cuantos_nodos_esenciales_k9=[]
cuantos_nodos_grado_k9=[]
lista_grados9=[]

for l in range(8):##Si es la red Y2H usar rango (10)
    cuantos_nodos_esenciales_k9.append(cuantos_nodos_esenciales_k[l])
cuantos_nodos_esenciales_k9=np.array(cuantos_nodos_esenciales_k9)
    
for l in range(8):
    cuantos_nodos_grado_k9.append(cuantos_nodos_grado_k[l])
cuantos_nodos_grado_k=np.array(cuantos_nodos_grado_k)

##Calculo la probabilidad de un nodo de ser esencial dado su grado k
proba_nodo_esencial_k9=(cuantos_nodos_esenciales_k9/cuantos_nodos_grado_k9)


for l in range(8):
    lista_grados9.append(lista_grados[l])
lista_grados9=np.array(lista_grados9)


##Calculo la probabilidad de un nodo de NO SER esencial dado su grado k
uno_menos_pe = 1-proba_nodo_esencial_k9



#Lo llevo a la forma que quiero representar
ln_uno_menos_pe= np.log(uno_menos_pe)

##Tengo que plotear proba_nodo_esencial_k vs lista de grados y ajustar por una lineal


x=lista_grados9
y=ln_uno_menos_pe

fitfunc = lambda p, x: p[0]*x+p[1]
p0 = [1, 1] 

errfunc = lambda p, x, y: fitfunc(p, x) - y 
out= optimize.leastsq(errfunc, p0[:], args=(x, y), full_output=1)
plt.xlabel("Grado (k)")
plt.ylabel("Ln(1-PE)")
plt.plot(x, y, "*", x, fitfunc(out[0], x), "k-", linewidth=1)
plt.plot(x, fitfunc(out[0], x), "k-", label = "Ln(1-PE)=-0.01476k -0.2010 ")
plt.legend()

#Obtengo los valores de pendiente y ordenada al origen del ajuste
a=out[0][0]
b=out[0][1]


err_a=np.sqrt(out[1][0][0])*out[0][0]
err_b=np.sqrt(out[1][1][1])
print("a=",out[0])
print("b=",out[1])

err_alfa=np.exp(a)*err_a
err_beta=np.exp(b)*err_b

##DA errores muy grandes en b
#Con los valores del a y b calculo los correspondientes alfa y beta del modelo de He
alfa=1-np.exp(a)
beta=1-np.exp(b)
print("Alfa=", alfa)
print("Beta=", beta)

#Expreso los valores porcentuales
alfa_porciento=alfa*100
beta_porciento=beta*100
print("Alfa%=", alfa_porciento)
print("Beta%=", beta_porciento)
