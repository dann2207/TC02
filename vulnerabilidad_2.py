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

# 