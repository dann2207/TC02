# Inciso (C)
import networkx as nx
import matplotlib.pyplot as plt
import random

G_APMS = nx.read_edgelist('yeast_AP-MS.txt')
G_LIT  = nx.read_edgelist('yeast_LIT.txt')
G_LITR = nx.read_edgelist('yeast_LIT_Reguly.txt')
G_Y2H  = nx.read_edgelist('yeast_Y2H.txt')

N_APMS = G_APMS.nodes()
#N_LIT  = G_LIT.nodes()
#N_LITR = G_LITR.nodes()
N_Y2H  = G_Y2H.nodes()

# Ahora veremos si algun par de redes comparte proteinas

N_set_APMS = {n for n in N_APMS}
N_set_Y2H = {n for n in N_Y2H}

N_set_APMS.intersection(N_set_Y2H)

# Con esto vimos que hay muchas proteinas en comun

# CAMBIAR EL CODIGO: ES IMPORTANTE QUE LOS NODOS "AL AZAR" QUE SACAMOS
# SEAN NO-ESENCIALES. EN EL PAPER DE ZOTENKO SE SACARON TODOS LOS NODOS
# ESENCIALES. LO QUE ESTAMOS HACIENDO ACA ES NO SACAR TODOS, SINO ALGUNOS,
# Y NOS QUEDARIA UNA TABLA EN DONDE EN VEZ DE NUMEROS TENDRIAMOS GRAFICOS.
# HAY UNA COMPLICACION CON CALCULAR EL DESVIO ESTANDAR DE ESTOS GRAFICOS,
# YA QUE NO HAY UNA UNICA FORMA DE ELEGIR ALGUNOS NODOS ESENCIALES AL AZAR.


# Ahora, pasamos a eliminar nodos esenciales de las redes
# Elegimos un numero M de nodos a sacar
# Sacamos, al azar, M nodos esenciales
# Luego, nos fijamos la distribucion de grados
# i.e. que grados aparecen entre los nodos que sacamos
# y cuantos nodos tienen cada uno de esos grados
# A continuacion, usando esto ultimo, sacams nodos al azar
# no necesariamente esenciales, que tengan la misma distribucion de grados
# Finalmente calculamos nuestra medida de impacto en ambos,
# que en este caso consiste en la fraccion de nodos en
# la componente mas grande

from funciones import ldata
PaHe = ldata("Essential_ORFs_paperHe.txt")

# La lista de los 1160 nodos esenciales es la siguiente
ess_list = [x[1] for x in PaHe]
# Por ejemplo, un elemento de esta lista es:
# PaHe[40][1] = 'YBR049C'

# Vamos a sacar 200 nodos
M = 200

# Tomamos una muestra al azar de M nodos esenciales
# y la guardamos en la siguiente lista
ess_list_sample = random.sample(ess_list,M)

# Ahora importamos una funcion que agarre la ultima
# lista y devuelva los grados con su numero de ocupacion
# i.e. dos listas de igual longitud (o una lista de tuplas)
from funciones import kdist
