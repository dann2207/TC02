# Si aun no estan los csv, hay que correr el codigo
# "Centralidad.py" antes de correr este
import pandas as pd
import matplotlib.pyplot as plt

df_APMS = pd.read_csv('G_APMS.csv', header=0, index_col=0)
df_LIT = pd.read_csv('G_LIT.csv', header=0, index_col=0)
df_Y2H = pd.read_csv('G_Y2H.csv', header=0, index_col=0)
df_LITR = pd.read_csv('G_LITR.csv', header=0, index_col=0)

DF = [df_APMS, df_LITR, df_LIT, df_Y2H]
Titulos = ["Red AP-MS","Red LIT","Red Y2H","Red LIT Reguly"]

ylabel = "Fraccion del tamaño de la \n componente gigante"
xlabel = "Fraccion de nodos removidos"
title = "Impacto de la remocion de nodos en funcion de la centralidad"
i = 0
for df in DF:
    i+=1
    plt.subplot(2,2,i)
    plt.plot(df.loc['X_random'],df.loc['Y_random'],label='Random')
    plt.plot(df.loc['X_degree'],df.loc['Y_degree'],label='Degree')
    plt.plot(df.loc['X_betweenness'],df.loc['Y_betweenness'],label='Betweenness')
    plt.plot(df.loc['X_eigenvector'],df.loc['Y_eigenvector'],label='Eigenvector')
    plt.title(Titulos[i-1])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid(True)
plt.suptitle(title,size=16)
plt.show()


