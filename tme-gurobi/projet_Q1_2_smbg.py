#!/usr/bin/python

# Copyright 2013, Gurobi Optimization, Inc.


from gurobipy import *

liste_xet_zet=[]
for k in range(1,7) :
    nbcont=6 
    nbvar=7

    # Range of plants and warehouses
    lignes = range(nbcont)
    colonnes = range(nbvar)

    # Matrice des contraintes
    a = [
        [1,-1,0,0,0,0,0],
        [1,0,-1,0,0,0,0],
        [1,0,0,-1,0,0,0],
        [1,0,0,0,-1,0,0],
        [1,0,0,0,0,-1,0],
        [1,0,0,0,0,0,-1]
        ]

    # Second membre
    b = [4,7,1,3,9,2]

    # Coefficients de la fonction objectif
    c = [k, -1, -1, -1, -1, -1, -1]

    m = Model("mogplex")     
                
    # declaration variables de decision
    x = []
    # r va de -inf à +inf
    x.append(m.addVar(vtype=GRB.CONTINUOUS, lb=-float("inf"), ub=float("inf"), name="r"))

    for i in range(1, nbvar):
        #les bik vont de -inf à 0
        x.append(m.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=float("inf"), name="b%d" % (i+1)))

    # maj du modele pour integrer les nouvelles variables
    m.update()

    #obj = LinExpr()
    obj =0
    for j in colonnes:
        obj += c[j] * x[j]
                
    # definition de l'objectif
    m.setObjective(obj,GRB.MAXIMIZE)

    # Definition des contraintes
    for i in lignes:
        m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) <= b[i], "Contrainte%d" % i)

    # Resolution
    m.optimize()

    
    # print("")                
    # print('Solution optimale:')
    xet = []
    for j in colonnes:
        # print('x%d'%(j+1), '=', x[j].x)
        xet.append(x[j].x)
    # print("")
    # print('Valeur de la fonction objectif :', m.objVal)
    
    liste_xet_zet.append((xet,m.objVal))

print(liste_xet_zet)
