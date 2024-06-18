parameters = {} # This is a dictionary stating parameter names and values you would like to consider to find the vectors that determine the kernel of the Jacobian matrices of the system

var = [] # List of variables passed as strings. Only provide the little letters. The code will add the variables in the bulk later

diffmatrix = [] # List of lists defining the full diffusion matrix of the system considering bulk and surface

kinetics = [] # List containing the kinetics of the full system. This includes the linear part in the bulk

lvalues = [] # List of numbers representing the values of l to be considered

tol = 1e-7 # Tolerance to check bifurcations and to check whether it is suitable to define a vector in terms of a particular submatrix of the non-invertible Jacobian matrix of the system

find_eq = 'n' # Boolean variable 'y/n' to specify whether you would like Python to obtain the base state for the parameter values provided in the variable 'parameters'

phiunit = 'n' # Boolean variable 'y/n' to specify whether you would like the script to normalize phi, the vector that spans the kernel of the non-invertible Jacobian matrix'

thirdcoef = 'y' # Boolean variable 'y/n' to specify whether you would like the script to compute second and third order coefficients

crosscoef = 'y' # Boolean variable 'y/n' to specify whether you would like the script to compute the coefficient associated with the unfolding C11

crosspar = 'c' # Parameter of the system with respect to which you would like to get the coefficient associated with the unfolding

plot2d = 'y' # Boolean variable 'y/n' to speficy whether you would like the script to save the variables to plot the bifurcation curves in Mathematica

parameters_on_axes = ['', ''] # List with the names of the parameters you would like to put in your bifurcation diagram. The first (respectively, second) one will be on the x-axis (respectively, y-axis)

names_of_parameters = [] # List with the names you would like to appear in the bifurcation diagram. If these coincide with the names given to the parameters, then this list can be empty

intervalx = [] # List with two entries representeing the minimum and maximum of the interval to consider in the x-axis in the bifurcation diagram

intervaly = [] # List with two entries representeing the minimum and maximum of the interval to consider in the y-axis in the bifurcation diagram

lines_to_search = {} # Dictionary to specify the lines you would like Mathematica to search bifurcation curves. This variable tells Mathematica to fix these parameter values separately and look for solutions
                     # to the bifurcation condition with respect to the other parameter. 