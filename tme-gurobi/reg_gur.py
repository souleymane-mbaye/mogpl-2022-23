#!/usr/bin/python

# Copyright 2013, Gurobi Optimization, Inc.


from gurobipy import *



nbcont=6 
nbvar=14

# Range of plants and warehouses
lignes = range(nbcont)
colonnes = range(nbvar)

# Matrice des contraintes
a = [[ 1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   4,  1],
     [ 0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  17,  1],
     [ 0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  37,  1],
     [ 0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  55,  1],
     [ 0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  88,  1],
     [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  -1, 96,  1]]

# Second membre
b = [11, 25, 46, 48, 65, 71]

# Coefficients de la fonction objectif
c = [1,1, 1,1, 1,1 ,1,1 ,1,1 ,1,1 ,0,0]

m = Model("mogplex")     
        
# declaration variables de decision
x = []
for i in colonnes:
    x.append(m.addVar(vtype=GRB.CONTINUOUS, lb=0, name="x%d" % (i+1)))

# maj du modele pour integrer les nouvelles variables
m.update()

obj = LinExpr();
obj =0
for j in colonnes:
    obj += c[j] * x[j]
        
# definition de l'objectif
m.setObjective(obj,GRB.MINIMIZE)

# Definition des contraintes
for i in lignes:
    m.addConstr(quicksum(a[i][j]*x[j] for j in colonnes) == b[i], "Contrainte%d" % i)

# Resolution
m.optimize()


print("")                
print('Solution optimale:')
for j in colonnes:
    print('x%d'%(j+1), '=', x[j].x)
print("")
print('Valeur de la fonction objectif :', m.objVal)

   




from numpy import arange, array, ones, linalg
from pylab import plot, show

xi = array([4,17,37,55,88,14]) # changing the point 96,71 to 14,97
A  = array([xi,ones(6)])

# linearly generated sequence
# y = [11,25,46,48,65,71]
y = [11,25,46,48,65,97] # changing the point 96,71 to 14,97

# obtaining the parameters
w = x[12].x, x[13].x

# plotting the line
line = w[0]*xi + w[1] # regression line
plot(xi, line, 'r-', xi, y, 'o')
show()