#!/usr/bin/python

# Copyright 2013, Gurobi Optimization, Inc.


from gurobipy import *
import numpy as np
import matplotlib.pyplot as plt
import time


x2 = []
z2 = []
temps = []

for n in range(5,20,5):
    
    # u = [325, 225, 210, 115, 75, 50]
    p = 5*n
    w_dep = sorted(np.random.randint(1,400,size=n))
    w =     w = [ w_dep[i]-w_dep[i+1] for i in range(len(w_dep)-1)] + [w_dep[-1]] # les w primes
    u = np.random.randint(1,400,size=n*p).reshape((n,p))
    range_n = range(n)
    range_p = range(p)
    
    print("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
    print("n:",n,"p:",p,"w:",w)
    # print(u)
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
    
    start = time.time()
    m = Model("mogplex")

    x = np.array(
        [[m.addVar(vtype=GRB.BINARY, name="x%d_%d" % (i+1, j+1))
          for j in range_p] for i in range_n]
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
        # ajout de la contrainte  zi=  somme de uij * xij sur j
        m.addConstr(z[i] - quicksum(u[i,j] * x[i][j]
                    for j in range_p) == 0, "Contrainte z %d" % (i+1))

    for k in range_n:
        for i in range_n:
            # ajout de la contrainte rk- bik <= zi
            m.addConstr(r[k] - b[i][k] - z[i] <= 0,
                        "Contrainte de b%d%d" % (i+1, k+1))

    for j in range_p:
        # ajout de contrainte de (somme de xij <= 1)
        m.addConstr(quicksum(x[i][j] for i in range_n) <=
                    1, "Contrainte affecation d'objet %d" % (j+1))

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
    x2.append([[x[i][j].x for j in range_p] for i in range_n])
    z2.append([z[i].x for i in range_n])
    
    end = time.time()
    temps.append(end-start)

# x2 = np.array(x2)
# z2 = np.array(z2)

plt.figure()
plt.plot([i for i in range(5,20,5)],temps)
plt.title("Evolution du temps en fonction de n et p")
plt.show()
