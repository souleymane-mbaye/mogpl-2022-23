from numpy import arange, array, ones, linalg
from pylab import plot, show

# xi = array([4,17,37,55,88,96])
xi = array([4,17,37,55,88,14]) # changing the point 96,71 to 14,97
A  = array([xi,ones(6)])

# linearly generated sequence
# y = [11,25,46,48,65,71]
y = [11,25,46,48,65,97] # changing the point 96,71 to 14,97

# obtaining the parameters
w = linalg.lstsq(A.T,y)[0]

# plotting the line
line = w[0]*xi + w[1] # regression line
plot(xi, line, 'r-', xi, y, 'o')
show()