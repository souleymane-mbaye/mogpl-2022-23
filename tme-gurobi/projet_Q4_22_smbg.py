#!/usr/bin/python

# Copyright 2013, Gurobi Optimization, Inc.


from gurobipy import *
import numpy as np

t = [ #   a         b       c       d          e        f       g
      [( 0, 0), ( 5, 3), (10, 4), ( 2, 6), ( 0, 0), ( 0, 0), (0, 0)],   # a
      [( 0, 0), ( 0, 0), ( 4, 2), ( 1, 3), ( 4, 6), ( 0, 0), (0, 0)],   # b
      [( 0, 0), ( 0, 0), ( 0, 0), ( 0, 0), ( 3, 1), ( 1, 2), (0, 0)],   # c
      [( 0, 0), ( 0, 0), ( 1, 4), ( 0, 0), ( 0, 0), ( 3, 5), (0, 0)],   # d
      [( 0, 0), ( 0, 0), ( 0, 0), ( 0, 0), ( 0, 0), ( 0, 0), (1, 1)],   # e
      [( 0, 0), ( 0, 0), ( 0, 0), ( 0, 0), ( 0, 0), ( 0, 0), (1, 1)],   # f
]
senario = 2 # le nombre de senarios
n = 7   # le nombre de sommets
range_n = range(n)
range_s = range(senario)
lw = [[1,1],[0,1/2]]

x2 = []
vf = []

for w in lw:    # un seul w
    
    m = Model("mogplex")

    dx = [[j for j in range_n if t[i][j][0]!=0] for i in range_n[:-1]]  # tableau des successeurs
    x = [[m.addVar(vtype=GRB.BINARY, name="x%d_%d" % (i+1, j+1))
          for j in dx[i]] for i in range_n[:-1]]   # arcs de poids non nul
    

    z = [m.addVar(vtype=GRB.CONTINUOUS, lb=- GRB.INFINITY, ub=0, name="z%d" % (i+1))
         for i in range_s]
    
    r = [m.addVar(vtype=GRB.CONTINUOUS, lb=- GRB.INFINITY,
                  name="r%d" % (k+1)) for k in range_s]
    
    b = [[m.addVar(vtype=GRB.CONTINUOUS, lb=0, name="b%d_%d" % (i+1, k+1))
          for k in range_s] for i in range_s]

    arcs = []
    for i in range_n[:-1]:
        for j in range(len(dx[i])):
            arcs.append((i,j))
    for s in range_s:
        # ajout de la contrainte  zi=  -somme de xij * tsij
        m.addConstr(z[s] + quicksum(x[i][j]*t[i][dx[i][j]][s] for (i,j) in arcs) 
                    == 0, "Contrainte z %d" % (s+1))

    for k in range_s:
        for i in range_s:
            # ajout de la contrainte rk-bik <= zi
            m.addConstr(r[k] - b[i][k] - z[i] <= 0,
                        "Contrainte de b%d%d" % (i+1, k+1))
    
    print("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
    print(w)
    print(arcs)
    print(dx)
    for i in range_n:
        t1 = [(j,dx[j].index(i)) for j in range_n[:-1] if i in dx[j]]
        print(t1)
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")

    # ajout de la contrainte  somme de x0j = 1 sur j, un seul arc sortant de la somme
    m.addConstr(quicksum(x[0][j]
                for j in range(len(x[0]))) == 1, "Contrainte un seul arc sortant de la source") # pas d'arc vers a

    for i in range_n[1:-1]: # que sur les sommets internes
        # ajout de contrainte de conservation (somme(j) de xji = somme(j) xij )
        m.addConstr(quicksum(x[i][j] for j in range(len(x[i]))) - quicksum(x[j][dx[j].index(i)] for j in range_n[:-1] if i in dx[j]) == 0
                    , "Contrainte de conservation sommet %d" % (i+1))

    for i in range_n[1:-1]:
        # ajout de la contrainte  de cyclage, un seul arc entrant au plus, somme(i) xij <= 1
        m.addConstr( quicksum(x[j][dx[j].index(i)] for j in range_n[:-1] if i in dx[j]) <= 1
                    , "Contrainte cyclage sommet %d" % (i+1))
    
    # maj du modele pour integrer les nouvelles variables
    m.update()

    # definition de l'objectif
    obj = LinExpr()
    obj = 0
    for k in range_s:
        obj += w[k]*((k+1)*r[k] - quicksum(b[i][k] for i in range_s))

    m.setObjective(obj, GRB.MAXIMIZE)

    # Resolution
    m.optimize()
    x1 = [[(i,dx[i][j]) for j in range(len(x[i])) if x[i][j].x!=0] for i in range_n[:-1]]
    x2.append(x1)
    v1,v2 = 0,0
    for tij in x1:
        if len(tij) > 0:
            i,j = tij[0][0],tij[0][1]
            v1 += t[i][j][0]
            v2 += t[i][j][1]
    vf.append((v1,v2))
# x2 = np.array(x2)

for i in range(len(lw)):
    print('\n------------------------------')
    print('  w = ', lw[i])
    print('  tij = ', vf[i])
    print('  x = \n', x2[i])
