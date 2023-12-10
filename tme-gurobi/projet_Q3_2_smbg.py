#!/usr/bin/python

# Copyright 2013, Gurobi Optimization, Inc.


from gurobipy import *
import numpy as np
import matplotlib.pyplot as plt
import time

temps = []

np.random.seed(91)
for n in [2,5,10]:
    for p in [5,10,15,20]:
        
        u = np.random.randint(1,21,size=n*p).reshape((n,p))
        c = np.random.randint(1,11,size=p)*10
        w = sorted(np.random.randint(1,n,n),reverse=True)
        budget = c.sum()/2
        range_n,range_p = range(n),range(p)
        
        start = time.time()

        m = Model("mogplex")

        x = np.array(
            [m.addVar(vtype=GRB.BINARY, name="x%d" % (j+1))
            for j in range_p]
        )

        z = np.array(
            [m.addVar(vtype=GRB.CONTINUOUS, lb=0, name="z%d" % (i+1))
            for i in range_n]
        )

        r = np.array(
            [m.addVar(vtype=GRB.CONTINUOUS, lb=- GRB.INFINITY,
                    name="r%d" % (k+1)) for k in range_n]
        )

        b = np.array(
            [[m.addVar(vtype=GRB.CONTINUOUS, lb=0, name="b%d_%d" % (i+1, k+1))
            for k in range_n] for i in range_n]
        )

        for i in range_n:
            # ajout de la contrainte  zi=  somme de uij * xj sur j
            m.addConstr(z[i] - quicksum(u[i][j] * x[j]
                        for j in range_p) == 0, "Contrainte z %d" % (i+1))

        
        print("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        print(n,p)
        print(len(c))
        print(len(w))
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
        
        for k in range_n:
            for i in range_n:
                # ajout de la contrainte rk- bik <= zi
                m.addConstr(r[k] - b[i][k] - z[i] <= 0,
                            "Contrainte de b%d%d" % (i+1, k+1))

        # ajout de contrainte de budget (somme de xj*cj <= budget)
        m.addConstr(quicksum(c[j]*x[j] for j in range_p) <=
                    budget, "Contrainte affecation d'objet xj")

        # maj du modele pour integrer les nouvelles variables
        m.update()

        # definition de l'objectif
        obj = LinExpr()
        obj = 0
        for k in range_n:
            obj += w[k]*((k+1)*r[k] - quicksum(b[i][k] for i in range_n))

        m.setObjective(obj, GRB.MAXIMIZE)

        # Resolution
        m.optimize()
        
        end = time.time()
        
        temps.append((n,p,end-start))

temps = np.array(temps)
plt.figure()
plt.plot(temps[:,0],temps[:,2])
plt.title("Temps en fonction de n")
plt.show()
tp = sorted(temps,key=lambda npt:npt[1])
tp = np.array(tp)
plt.plot(tp[:,1],tp[:,2])
plt.title("Temps en fonction de p")
plt.show()
tnp = sorted(temps,key=lambda npt:npt[0]+npt[1])
tnp = np.array(tnp)
plt.plot(tnp[:,0]+tnp[:,1],tnp[:,2])
plt.title("Temps en fonction de n+p")
plt.show()
tnp = sorted(temps,key=lambda npt:npt[0]*npt[1])
tnp = np.array(tnp)
plt.plot(tnp[:,0]*tnp[:,1],tnp[:,2])
plt.title("Temps en fonction de n*p")
plt.show()
