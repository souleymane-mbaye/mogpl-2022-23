#!/usr/bin/python

# Copyright 2013, Gurobi Optimization, Inc.


from gurobipy import *
import numpy as np

u = [[19, 6, 17, 2], [2, 11, 4, 18]]
c = [40, 50, 60, 50]
budget = 100
w_dep = [[2, 1], [10, 1], [1/2, 1/2]]
w = [[1, 1], [9, 1], [0, 1/2]]
n, p = 2, 4
range_n,range_p = range(n),range(p)

x2,z2 = [],[]

for indice_w in range(len(w)):
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
        obj += w[indice_w][k]*((k+1)*r[k] - quicksum(b[i][k] for i in range_n))

    m.setObjective(obj, GRB.MAXIMIZE)

    # Resolution
    m.optimize()
    x2.append([x[j].x for j in range_p])
    z2.append([z[i].x for i in range_n])

x2 = np.array(x2)
z2 = np.array(z2)

for i in range(len(w)):
    print('\n------------------------------')
    print('  w = ', w_dep[i])
    print('  z = ', z2[i])
    print('  x = ', x2[i])
