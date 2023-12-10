#!/usr/bin/python

# Copyright 2013, Gurobi Optimization, Inc.


from gurobipy import *
import numpy as np

u = [325, 225, 210, 115, 75, 50]
w_dep = [[3, 2, 1], [10, 3, 1], [1/3, 1/3, 1/3]]
w = [[1, 1, 1], [7, 2, 1], [0, 0, 1/3]]
n, p = 3, 6
range_n = range(n)
range_p = range(p)

x2 = []
z2 = []

for indice_w in range(len(w)):
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
        m.addConstr(z[i] - quicksum(u[j] * x[i][j]
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
        obj += w[indice_w][k]*((k+1)*r[k] - quicksum(b[i][k] for i in range_n))

    m.setObjective(obj, GRB.MAXIMIZE)

    # Resolution
    m.optimize()
    x2.append([[x[i][j].x for j in range_p] for i in range_n])
    z2.append([z[i].x for i in range_n])

x2 = np.array(x2)
z2 = np.array(z2)

for i in range(len(w)):
    print('\n------------------------------')
    print('  w = ', w_dep[i])
    print('  z = ', z2[i])
    print('  x = \n', x2[i])
