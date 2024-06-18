from IPython import get_ipython
get_ipython().run_line_magic('reset', '-sf')

import os
import shutil
from sys import exit
from sympy import *
import numpy as np
init_printing()
from sympy.solvers import solve
import math
import json
from mpmath import findroot
import copy

# The script runs the file 'name_of_model.py'.

modelname = input('Enter the name of the system you would like to analyze: ')

if not os.path.isdir(modelname):
    print('The directory of that system was not found. Create it first and place ' + 
          'the data file inside named in the same way.')
    exit()
else:
    os.chdir(os.getcwd() + '\\' + modelname)

try:
    exec(open(modelname + '.py').read())
except:
    print('The file ' + modelname + '.py could not be run')
    exit()
    
try:
    var
except:
    print('Variables were not provided')
    exit()
    
nvar = len(var)

# The script runs the file 'functions.py'

try:
    exec(open(os.path.dirname(os.path.realpath(__file__)) + '\\functions.py').read())
except:
    print('File functions.py is not in the same folder as the script you are running')
    exit()

file = open('List of Variables.txt', 'w')

# The code defines variables in the bulk by capital letter and variables on the surface with lower case

newvar = [None]*2*nvar
for varnum in range(nvar):
    try:
        exec(var[varnum].upper() + ' = symbols(var[varnum].upper(), real = True)')
        exec(var[varnum].lower() + ' = symbols(var[varnum].lower(), real = True)')
        newvar[varnum] = eval(var[varnum].upper())
        newvar[nvar+varnum] = eval(var[varnum].lower())
    except:
        print('The script could not define ' + var[varnum] + ' as a variable')
        exit()
var = Matrix(newvar)

for varnum in range(len(var)):
    file.write(latex(var[varnum]) + '\n')
file.close()
    
# The script defines the parameters as symbols in sympy.

try:
    parameters
except:
    print('Parameters were not provided. The script will assume that there are no parameters.')
    parameters=[]

npar = len(parameters)

if npar>0:
    newparameters=dict()
    for key in parameters.keys():
        try:
            exec(key + ' = symbols(key, real = True)')
            newparameters[eval(key)] = parameters[key]
        except:
            print('The script could not define your variable ' + key + ' as a variable')
            exit()
    parameters = newparameters
    
# It defines the diffusion matrix as a symbolic matrix in terms of the parameters of the system.

try:
    diffmatrix
except:
    print('Diffusion matrix was not provided.')
    exit()
    
diffmatrix_for_mathematica = copy.deepcopy(diffmatrix)
try:
    for row in range(2*nvar):
        for col in range(2*nvar):
            diffmatrix[row][col] = eval(str(diffmatrix[row][col]).replace('^','**'))
            diffmatrix_for_mathematica[row][col] = latex(diffmatrix[row][col])
except:
    print('The diffusion matrix is not a function of the parameters of the system')
    exit()

diffmatrix = Matrix(diffmatrix)

file = open('Diffusion matrix for Mathematica.json','w')
json.dump(diffmatrix_for_mathematica,file)
file.close()

# Kinetics

# The script defines the list of kinetics as a symbolic matrix in terms of the parameters of the system.

try:
    kinetics
except:
    print('Kinetics were not provided.')
    exit()

for functionnumber in range(2*nvar):
    try:
        kinetics[functionnumber] = eval(str(kinetics[functionnumber]).replace('^', '**'))
    except:
        print('The expression ' + kinetics[functionnumber] + 'is not a function of the parameters of your system')

kinetics = Matrix(kinetics)

# The script writes the list of kinetics functions in LaTeX in a text file.

file=open('Kinetics.txt','w')
for functionnumber in range(2*nvar):
    file.write(latex(kinetics[functionnumber]) + '\n')
file.close()

jacobianmat = kinetics.jacobian(var)

# The script checks whether B, K and D_U are diagonal

if np.count_nonzero(jacobianmat[0:nvar, 0:nvar] - np.diag(np.diagonal(jacobianmat[0:nvar, 0:nvar])))>0:
    print('B is not diagonal')
    exit()
if np.count_nonzero(jacobianmat[nvar:2*nvar, 0:nvar] - np.diag(np.diagonal(jacobianmat[nvar:2*nvar, 0:nvar])))>0:
    print('K is not diagonal')
    exit()
if np.count_nonzero(diffmatrix[0:nvar, 0:nvar] - np.diag(np.diagonal(diffmatrix[0:nvar, 0:nvar])))>0:
    print('The diffusion matrix in the bulk is not diagonal')
    exit()

# The script defines symbols and variables that will be used to find the coefficients of the normal form

rNF = symbols('rNF')

BNF = - jacobianmat[0:nvar, 0:nvar]

KNF = jacobianmat[nvar:2*nvar, 0:nvar]

omegaNF = matrix('omegaNF', nvar, eye(nvar))

for varnum in range(nvar):
    omegaNF.actualcoord[varnum, varnum] = sqrt(Mul(BNF[varnum, varnum], Pow(diffmatrix[varnum, varnum], - 1)))
    
i0omegar = zeros(nvar)
i0primeomegaR = zeros(nvar)

for varnum in range(nvar):
    i0omegar[varnum, varnum] = i_n(0, Mul(omegaNF.dummy[varnum, varnum], rNF))
    i0primeomegaR[varnum, varnum] = i_n(0, Mul(omegaNF.dummy[varnum, varnum], R), derivative = True)
    
i0omegaR = i0omegar.subs(rNF, R)

diagonal0 = matrix('diag', nvar, Add(Mul(Matrix(diffmatrix[0:nvar, 0:nvar]), i0primeomegaR, omegaNF.dummy),
                                     Mul(KNF, i0omegaR)))

S0 = matrix('S0', nvar, Mul(diagonal0.dummy.inv(), KNF))

Ustar = Mul(i0omegar, S0.dummy, Matrix(var[nvar:2*nvar]))

# The script finds a base state after evaluating the kinetics if requested.

if find_eq=='y':
    kineticsevaluated = Matrix(kinetics[nvar:2*nvar])
    for varnum in range(nvar):
        kineticsevaluated = kineticsevaluated.subs(var[varnum], Ustar[varnum])
    for varnum in range(nvar):
        kineticsevaluated = kineticsevaluated.subs(S0.dummy[varnum, varnum], S0.actualcoord[varnum, varnum])
    for varnum in range(nvar):
        kineticsevaluated = kineticsevaluated.subs(diagonal0.dummy[varnum, varnum], diagonal0.actualcoord[varnum, varnum])
    for varnum in range(nvar):
        kineticsevaluated = kineticsevaluated.subs(omega[varnum,varnum], omegaeval[varnum, varnum])
    kineticsevaluated = kineticsevaluated.subs(r, R)
    for key in parameters:
        kineticsevaluated = kineticsevaluated.subs(key, parameters[key])
    eq = solve(kineticsevaluated, var)

    while eq==[]:
        print('This script could not find an equilibrium point for the given parameter values shown below: ')
        print(parameters)
        parameterchange = input('Would you like to change a parameter value? [y/n] ')
        while parameterchange!='y' and parameterchange!='n':
            parameterchange = input('You must enter your answer in the format [y/n]. Would you like to change a parameter value? ')
        if parameterchange=='y':
            whichparam = input('Which parameter would you like to change? ')
            while True:
                try:
                    whichparam = eval(whichparam)
                    while whichparam not in parameters.keys():
                        whichparam = input('That is not a parameter of the system. Which parameter would you like to change? ')
                    break
                except:
                    whichparam = input('That is not a parameter of the system. Which parameter would you like to change? ')
            if whichparam in parameters.keys():
                parameters[whichparam] = input('Enter a value of ' + whichparam + ': ')
                while not isfloat(parameters[whichparam]):
                    parameters[whichparam] = input('What you entered before is not a number. Enter a value of ' + whichparam + ': ')
                parameters[whichparam] = eval(parameters[whichparam])
            kineticsevaluated=kinetics
            for key in parameters:
                kineticsevaluated = kineticsevaluated.subs(key, parameters[key])
            eq = solve(kineticsevaluated,var)
        else:
            break
    if isinstance(eq,list) and len(eq)>1:
        print('Your system has ' + str(len(eq)) + ' equilibria for the given parameter values given by:')
        for eqnum in range(len(eq)):
            print(str(eqnum+1) + '.- ' + str(eq[eqnum][nvar:2*nvar]) + '\n')
        eqnumber = input('Enter the number of the equilibrium point you want to consider: ')
        while not eqnumber.isnumeric() or int(eqnumber)==0:
            eqnumber = input('What you entered before is not a positive integer. Enter the number of the equilibrium point' +
                             ' you want to consider: ')
        eqnumber = int(eqnumber)
        eq = Matrix(list(eq[eqnumber-1]))        
    elif isinstance(eq,list) and len(eq)==1:
        print('Your system has only one equilibrium point given by:')
        print(eq[0])
        eq = Matrix(eq[0])
    
else:
    eq = []
                    
# The script obtains derivatives of the nonlinear vector field up to order three

negativeRHS = Vector('negativeRHS')

firstorderderivatives = list()
secondorderderivatives = list()
thirdorderderivatives = list()
for counter1 in range(nvar, 2*nvar):
    firstorderderivatives.append(diff(Add(Matrix(kinetics[nvar:2*nvar]),
                                          Mul(KNF, Matrix(var[nvar:2*nvar]))), var[counter1]))
    secondorderderivatives.append(list())
    thirdorderderivatives.append(list())
    for counter2 in range(nvar, 2*nvar):
        secondorderderivatives[counter1 - nvar].append(diff(firstorderderivatives[counter1 - nvar],
                                                            var[counter2]))
        thirdorderderivatives[counter1 - nvar].append(list())
        for counter3 in range(nvar, 2*nvar):
            thirdorderderivatives[counter1 - nvar][counter2 - nvar].append(diff(
                secondorderderivatives[counter1 - nvar][counter2 - nvar], var[counter3]))

# firstorderderivatives does not contain the elements with K

x = symbols('x')

# The script saves the values of l that are considered in the study

file = open('Values of l.txt', 'w')
for counter in range(len(lvalues)):
    file.write(str(lvalues[counter]) + '\n')
file.close()

sNF = symbols('sNF')

# The calculation below is done for each value of l provided

for lNF in lvalues:
    phiNF = Vector('phiNF')
    psiNF = Vector('psiNF')
    
    if not os.path.isdir('l=' + str(lNF)):
        os.mkdir('l=' + str(lNF))
    os.chdir(os.getcwd() + '\\l=' + str(lNF))
    
    # The script creates matrices of coefficients that will be used to obtain third-order coefficients
    
    ilomegar = matrix('ilomegar', nvar)
    ipomegaR = []
    ipomegaR_eval = []
    ipprimeomegaR = []
    ipprimeomegaR_eval = []
    
    ipuneval = []
    ipprimeuneval = []
    
    diagonalp = []
    diagonalp_eval = []
    Sp = []
    
    coefmatp = []
    
    for pNF in range(2*lNF + 1):
        ipomegaR.append(matrix('iomegaR' + str(pNF), nvar))
        ipprimeomegaR.append(matrix('iprimeomegaR' + str(pNF), nvar))
        
        ipuneval.append(i_n(pNF, rNF))
        ipprimeuneval.append(i_n(pNF, rNF, derivative = True))
        
        for varnum in range(nvar):
            ipomegaR[pNF].actualcoord[varnum, varnum] = ipuneval[pNF].subs(rNF,
                                                                           Mul(omegaNF.dummy[varnum, varnum], R))
            ipprimeomegaR[pNF].actualcoord[varnum, varnum] = ipprimeuneval[pNF].subs(rNF,
                                                                                     Mul(omegaNF.dummy[varnum, varnum],
                                                                                         R))
            
        diagonalp.append(matrix('diag' + str(pNF), nvar, Add(Mul(Matrix(diffmatrix[0:nvar, 0:nvar]),
                                                                 ipprimeomegaR[pNF].dummy, omegaNF.dummy),
                                                             Mul(KNF, ipomegaR[pNF].dummy))))
        
        Sp.append(Mul(diagonalp[pNF].dummy.inv(), KNF))
        
        coefmatp.append(matrix('coefmatp', nvar, Add(Mul(KNF, ipomegaR[pNF].dummy, Sp[pNF]),
                                                     Matrix(jacobianmat[nvar:2*nvar, nvar:2*nvar]),
                                                     Mul(- 1, pNF, pNF + 1, Pow(R, - 2),
                                                         Matrix(diffmatrix[nvar:2*nvar,
                                                                           nvar:2*nvar])))))
        
        diagonalp_eval.append(evaluation_dict(diagonalp[pNF], 'matrix'))
        ipomegaR_eval.append(evaluation_dict(ipomegaR[pNF], 'matrix'))
        ipprimeomegaR_eval.append(evaluation_dict(ipprimeomegaR[pNF], 'matrix'))

    for varnum in range(nvar):
        ilomegar.actualcoord[varnum, varnum] = ipuneval[lNF].subs(rNF, Mul(omegaNF.dummy[varnum, varnum], rNF))        
        
    omegaNF_eval = evaluation_dict(omegaNF, 'matrix')
        
    # The term - K is hidden inside jacobianmat
    
    # The script obtains and saves the equation for a Turing bifurcation
    # The Jacobian matrix can be complicated algebraically so the script computes the determinant of a
    # dummy matrix to simplify the process.
    
    criticalmatdet = coefmatp[lNF].dummy.det()
    
    for row in range(nvar):
        for col in range(nvar):
            criticalmatdet = criticalmatdet.subs(coefmatp[lNF].dummy[row, col],
                                                 coefmatp[lNF].actualcoord[row, col])
    for varnum in range(nvar):
        criticalmatdet = criticalmatdet.subs(diagonalp_eval[lNF]).subs(ipomegaR_eval[lNF])\
            .subs(ipprimeomegaR_eval[lNF]).subs(omegaNF_eval)
            
    # If the script finds at least one steady state, it evaluates the determinant and the Jacobian
    # matrix at the parameters provided to see whether there is a bifurcation for the value of l the code
    # is dealing with.
    
    determinanteval = criticalmatdet
    if eq!=[]:
        for varnum in range(nvar):
            determinanteval = determinanteval.subs(var[nvar + varnum], eq[nvar + varnum])
        for key in parameters:
            determinanteval = determinanteval.subs(key, parameters[key])
        tol=1e-6
        if abs(N(determinanteval))<tol:
            print('There is a Turing bifurcation for l = ' + str(lNF) + '.')
        
    file = open('Determinant l=' + str(lNF) + '.txt', 'w')
    file.write(latex(criticalmatdet))
    file.close()
    
    if thirdcoef=='y' or crosscoef=='y':
        getout = 0
        
        # The script looks for an invertible (n - 1) x (n - 1) submatrix of the Jacobian matrix of the system
        # to obtain a vector that spans its kernel
        
        if eq!=[]:
            for row in range(nvar):
                for col in range(nvar):
                    submatrixrows = list(range(nvar))
                    submatrixcols = list(range(nvar))
                    submatrixrows.remove(row)
                    submatrixcols.remove(col)
                    invertiblesubmatrix = criticalmat.extract(submatrixrows, submatrixcols)
                    submatrixeval = invertiblesubmatrix
                    for varnum in range(nvar):
                        submatrixeval = submatrixeval.subs(diagonalp_eval[lNF]).subs(ipomegaR_eval[lNF])\
                            .subs(ipprimeomegaR_eval[lNF]).subs(omegaNF_eval)
                    submatrixeval = submatrixeval.subs(parameters)
                    if eq!=[]:
                        for varnum in range(nvar):
                            submatrixeval = submatrixeval.subs(var[nvar + varnum], eq[nvar + varnum])
                    if abs(N(submatrixeval.det()))>tol:
                        phiNF.actualcoord[col] = 1
                        criticalrow = row
                        criticalcol = col
                        getout = 1
                        break
                if getout==1:
                    break
        else:
            submatrixrows = list(range(1, nvar))
            submatrixcols = list(range(1, nvar))
            invertiblesubmatrix = coefmatp[lNF].actualcoord.extract(submatrixrows, submatrixcols)
            phiNF.actualcoord[0] = 1
            criticalrow = 0
            criticalcol = 0
            
        coefsubmatrix = matrix('dummysubmatrix', nvar - 1, invertiblesubmatrix)
        
        # The script finds the vector \phi_1^1 that spans the kernel of the Jacobian matrix using dummy matrices.
        
        auxiliaryterm, = linsolve(Add(Mul(coefsubmatrix.dummy, Matrix(phiNF.actualcoord).extract(submatrixcols, [0])), \
                                      coefmatp[lNF].dummy.extract(submatrixrows, [criticalcol])),
                                  list(Matrix(phiNF.actualcoord).extract(submatrixcols, [0])))
            
        phiNF.actualcoord[0:criticalcol] = auxiliaryterm[0:criticalcol]
        phiNF.actualcoord[criticalcol + 1:nvar] = auxiliaryterm[criticalcol:nvar - 1]
        
        phiNF.actualcoord = Matrix(phiNF.actualcoord)
        
        for row in range(nvar):
            for col in range(nvar):
                phiNF.actualcoord = phiNF.actualcoord.subs(coefmatp[lNF].dummy[row, col],
                                                           coefmatp[lNF].actualcoord[row, col])
                if row<nvar-1 and col<nvar-1:
                    phiNF.actualcoord = phiNF.actualcoord.subs(coefsubmatrix.dummy[row, col],
                                                               coefsubmatrix.actualcoord[row, col])
                    
        if phiunit=='y':
            phiNF.actualcoord = Mul(Pow(phiNF.actualcoord.dot(phiNF.actualcoord), - 1), phiNF.actualcoord)
        
        print('phiNF ready')
        
        # The script finds the vector \psi that spans the Kernel of the adjoint of the Jacobian matrix
        
        Tdummycoefsubmatrix = transpose(coefsubmatrix.dummy)
        Tdummycoefmatl = transpose(coefmatp[lNF].dummy)
        
        psiNF.actualcoord[criticalrow] = 1
        
        auxiliaryterm, = linsolve(Add(Mul(Tdummycoefsubmatrix, Matrix(psiNF.actualcoord).extract(submatrixrows,
                                                                                                 [0])),
                                      Tdummycoefmatl.extract(submatrixcols, [criticalrow])),
                                  list(Matrix(psiNF.actualcoord).extract(submatrixrows, [0])))
            
        psiNF.actualcoord[0:criticalrow] = auxiliaryterm[0:criticalrow]
        psiNF.actualcoord[criticalrow+1:nvar] = auxiliaryterm[criticalrow:nvar-1]
        
        psiNF.actualcoord = Matrix(psiNF.actualcoord)
        
        for row in range(nvar):
            for col in range(nvar):
                psiNF.actualcoord = psiNF.actualcoord.subs(coefmatp[lNF].dummy[row, col],
                                                           coefmatp[lNF].actualcoord[row, col])
                if row<nvar - 1 and col<nvar - 1:
                    psiNF.actualcoord = psiNF.actualcoord.subs(coefsubmatrix.dummy[row, col],
                                                               coefsubmatrix.actualcoord[row, col])
        
        print('psiNF ready')
        
        phiNF_eval = evaluation_dict(phiNF)
        psiNF_eval = evaluation_dict(psiNF)
        
        integral = I_l()
        
        print('Hellish integral done')
    
    if thirdcoef=='y':
        DS_phiphi = secondorderapplied(phiNF, phiNF)
            
        for mNF in range(0, lNF + 1):
            
            # If l is even, the script finds second-order coefficients
            
            if lNF%2==0:
                
                # It obtains the set of two grouped numbers that add up to m
                
                list_of_q = add_2_up_to_m(lNF, mNF)
                
                for part_list in list_of_q:
                    part_list.sort()
                    
                    q_1 = part_list[0]
                    q_2 = part_list[1]
                    
                    dl = d_p(lNF, mNF, lNF, q_1, q_2)
                    
                    C2 = Mul(Pow(R, 2), Pow(integral, - 1), dl, psiNF.dummy.dot(DS_phiphi))
                    
                    if q_1!=q_2:
                        C2 = Mul(2, C2)
                    
                    C2 = C2.subs(phiNF_eval).subs(psiNF_eval).subs(diagonalp_eval[lNF]).subs(ipomegaR_eval[lNF])\
                        .subs(ipprimeomegaR_eval[lNF]).subs(omegaNF_eval)
                
                    file = open('C2 q1=' + str(q_1) + ', q2=' + str(q_2) + ', m=' + str(mNF) + '.txt','w')
                    file.write(latex(C2))
                    file.close()
                
                    print('Second order coefficient ready')
                    
            # The script continues to get third-order coefficients
                
            TS_phiphiphi = thirdorderapplied(phiNF, phiNF, phiNF)
                    
            # It obtains the set of three grouped numbers that add up to m
            
            list_of_q = add_3_up_to_m(lNF, mNF)
            
            for part_list in list_of_q:
                
                bigsum = 0
                
                for combnum in range(len(set(part_list))):
                    if len(set(part_list))==3:
                        q_1 = part_list[(0 + combnum)%3]
                        q_2 = part_list[(1 + combnum)%3]
                    
                    elif len(set(part_list))==2:
                        index = non_repeated_element(part_list)
                        
                        if combnum==0:
                            q_1 = part_list[(index + 0)%3]
                            q_2 = part_list[(index + 1)%3]
                        else:
                            q_1 = part_list[(index + 1)%3]
                            q_2 = part_list[(index + 2)%3]
                            
                    else:
                        q_1 = part_list[0]
                        q_2 = part_list[1]
                    
                    up2 = []
                    
                    DS_phiup2 = []
                
                    for pNF in range(abs(q_1 + q_2), 2*lNF + 1):
                        
                        # For each value of p, the code obtains the solutions to the second order equations
                        
                        dp = d_p(lNF, mNF, pNF, q_1, q_2)
                        
                        up2.append(Vector('up2^' + str(pNF)))
                        
                        if dp==0:
                            up2[- 1].dummy = Matrix([0]*nvar)
                            up2[- 1].actualcoord = Matrix([0]*nvar)
                        else:
                            negativeRHS.actualcoord = Mul(dp, DS_phiphi)
                            if lNF%2==0 and pNF==lNF:
                                up2[- 1].actualcoord = critical_linearsolver(up2[- 1], negativeRHS,
                                                                             criticalcol, coefsubmatrix,
                                                                             submatrixrows, submatrixcols)
                            else:
                                up2[- 1].actualcoord = linearsolver(up2[- 1], negativeRHS, coefmatp[pNF])
                                
                            up2[- 1].actualcoord = up2[- 1].actualcoord.subs(diagonalp_eval[pNF])\
                                .subs(ipomegaR_eval[pNF]).subs(ipprimeomegaR_eval[pNF])
                    
                        DS_phiup2.append(secondorderapplied(phiNF, up2[pNF - abs(q_1 + q_2)]))
                    
                    print('Second order ready')
                        
                    if len(set(part_list))==3:
                        q_3 = part_list[(2 + combnum)%3]
                        
                    elif len(set(part_list))==2:
                        if combnum==0:
                            q_3 = part_list[(index + 2)%3]
                        else:
                            q_3 = part_list[(index + 0)%3]
                            
                    else:
                        q_3 = part_list[2]
                        
                    summation = first_term_C3(lNF, mNF, q_1, q_2, q_3)
                    
                    for pNF in range(abs(q_1 + q_2), 2*lNF + 1):
                        for varnum in range(nvar):
                            summation = summation.subs(up2[pNF - abs(q_1 + q_2)].dummy[varnum],
                                                       up2[pNF - abs(q_1 + q_2)].actualcoord[varnum])
                    
                    print('Hard part ready')
                    
                    if q_1!=q_2:
                        summation = Mul(2, summation)
                    
                    bigsum = Add(bigsum, summation)
                    
                # The script computes the third order coefficient for the value of q_1, q_2, q_3 it is dealing
                # with
                    
                legendreintegral = fourth_integral(lNF, mNF, q_1, q_2, q_3)
                
                second_summand = Mul(TS_phiphiphi, legendreintegral)
                
                if len({q_1, q_2, q_3})==3:
                    second_summand = Mul(6, second_summand)
                elif len({q_1, q_2, q_3})==2:
                    second_summand = Mul(3, second_summand)
                    
                # The script computes the third order coefficient and saves it into a .txt file
                
                C3 = Mul(Pow(R, 2), Pow(integral, - 1), psiNF.dummy.dot(Add(Mul(2, bigsum), second_summand)))
                
                C3 = C3.subs(phiNF_eval).subs(psiNF_eval).subs(diagonalp_eval[lNF]).subs(ipomegaR_eval[lNF])\
                    .subs(ipprimeomegaR_eval[lNF]).subs(omegaNF_eval)
                    
                sorting = [q_1, q_2, q_3]
                sorting.sort()
                    
                file = open('C3 q1=' + str(sorting[0]) + ', q2=' + str(sorting[1]) + ', q3='\
                            + str(sorting[2]) + ', m=' +
                            str(mNF) + '.txt', 'w')
                file.write(latex(C3))
                file.close()
                        
                print('Third order coefficient ready')
    
    # If requested, the code computes the terms necessary for Mathematica to get the cross-order coefficient
    
    if crosscoef=='y':
        try:
            crosspar = eval(crosspar)
        except:
            if lNF==lvalues[0]:
                print('The Cross parameter was not provided in the right way')
        
        firsthellishterm = 0
        
        if diff(BNF, crosspar)!=zeros(nvar):
            for varnum in range(nvar):
                firsthellishterm = Add(firsthellishterm,
                                       Mul(- diff(B[varnum, varnum], crosspar), Pow(Sp[lNF][varnum, varnum], 2),
                                           phiNF.dummy[varnum], psiNF.dummy[varnum],
                                           Pow(omegaNF.dummy[varnum, varnum], - 3),
                                           integrate(Mul(Pow(ipomegar[lNF].actualcoord[varnum, varnum].
                                                             subs(r, Mul(s, Pow(omegaNF.dummy[varnum, varnum], - 1))),
                                                             2), Pow(s, 2)), (s, 0, Mul(omegaNF.dummy[varnum, varnum], R)))))
        
        if diff(KNF, crosspar)!=zeros(nvar):
            firsthellishterm = Add(firsthellishterm, Mul(Pow(R, 2),
                                                         psiNF.dummy.dot(Mul(- diff(KNF, crosspar),
                                                                             Add(phiNF.dummy,
                                                                                 Mul(- 1, ipomegaR[lNF].dummy,
                                                                                     Sp[lNF], phiNF.dummy))))))
        
        secondhellishterm = 0
        
        if diff(Matrix(diffmatrix[0:nvar, 0:nvar]), crosspar)!=zeros(nvar):
            for varnum in range(nvar):
                secondhellishterm = Add(secondhellishterm,
                                        Mul(diff(diffmatrix[varnum, varnum], crosspar), Pow(Sp[lNF][varnum, varnum], 2),
                                            phiNF.dummy[varnum], psiNF.dummy[varnum], Pow(omegaNF.dummy[varnum, varnum], - 3),
                                            integrate(Mul(Pow(ipomegar[lNF].actualcoord[varnum, varnum].
                                                              subs(r, Mul(s, Pow(omegaNF.dummy[varnum, varnum], -1))), 2),
                                                          Pow(s, 2)), (s, 0, Mul(omegaNF.dummy[varnum, varnum], R)))))
                
        if diff(Matrix(diffmatrix[nvar:2*nvar, nvar:2*nvar]), crosspar)!=zeros(nvar):
            secondhellishterm = Add(secondhellishterm,
                                    Mul(Pow(R, 2),
                                        psiNF.dummy.dot(Mul(diff(diffmatrix[nvar:2*nvar, nvar:2*nvar], crosspar),
                                                            Matrix(phiNF.dummy)))))
        
        extraterm = Mul(Pow(R, 2), Pow(integral, - 1))
        
        termtodiff = crossorderapplied(phiNF)
        termtodiff = psiNF.dummy.dot(termtodiff)
        
        # The script computes and saves the variables needed to compute C11 if requested
            
        C11 = Mul(Pow(integral, - 1), Add(firsthellishterm, Mul(- lNF, lNF + 1, Pow(R, - 2), secondhellishterm)))
        
        C11 = C11.subs(phiNF_eval).subs(psiNF_eval).subs(diagonalp_eval[lNF]).subs(ipomegaR_eval[lNF])\
            .subs(ipprimeomegaR_eval[lNF]).subs(omegaNF_eval)
            
        extraterm = extraterm.subs(phiNF_eval).subs(psiNF_eval).subs(diagonalp_eval[lNF])\
            .subs(ipomegaR_eval[lNF]).subs(ipprimeomegaR_eval[lNF]).subs(omegaNF_eval)
        
        phiNF.actualcoord = phiNF.actualcoord.subs(diagonalp_eval[lNF]).subs(ipomegaR_eval[lNF])\
            .subs(ipprimeomegaR_eval[lNF]).subs(omegaNF_eval)
            
        psiNF.actualcoord = psiNF.actualcoord.subs(diagonalp_eval[lNF]).subs(ipomegaR_eval[lNF])\
            .subs(ipprimeomegaR_eval[lNF]).subs(omegaNF_eval)
            
        
        file = open('Simple part of cross-order coefficient l=' + str(lNF) + '.txt', 'w')
        file.write(latex(C11))
        file.close()
        
        file = open('Kernels.txt', 'w')
        for varnum in range(nvar):
            if varnum<nvar - 1:
                file.write(latex(phiNF.actualcoord[varnum]) + ',')
            else:
                file.write(latex(phiNF.actualcoord[varnum]) + '\n')
        for varnum in range(nvar):
            if varnum<nvar - 1:
                file.write(latex(psiNF.actualcoord[varnum]) + ',')
            else:
                file.write(latex(psiNF.actualcoord[varnum]))
        file.close()
        
        file = open('Factor to get cross-order coefficient l=' + str(lNF) + '.txt', 'w')
        file.write(latex(extraterm))
        file.close()
        
        file = open('Term to differentiate to get cross-order coefficient l=' + str(lNF) + '.txt', 'w')
        file.write(latex(crosspar) + '\n')
        file.write(latex(termtodiff))
        file.close()
    
    os.chdir(os.path.dirname(os.getcwd()))
    
# The sript saves text files with all the information required by Mathematica to find and plot the
# bifurcation diagram
    
if plot2d=='y':
    try:
        for parnum in range(2):
            parameters_on_axes[parnum] = eval(parameters_on_axes[parnum])
    except:
        print('The variable párameters_on_axes is not well defined.')
        exit()

    try:
        file = open('Initial conditions for Turing bifurcation curves.txt','w')
        
        for key in lines_to_search.keys():
            if isinstance(lines_to_search[key],list):
                for initialsolnum in range(len(lines_to_search[key])):
                    if isfloat(str(lines_to_search[key][initialsolnum])):
                        file.write(latex(eval(key)) + ',' + 
                                   latex(eval(str(lines_to_search[key][initialsolnum]))) + '\n')
            else:
                if isfloat(str(lines_to_search[key])):
                    file.write(latex(eval(key)) + ',' + latex(eval(str(lines_to_search[key]))) + '\n')
                
        file.close()    
    except:
        print('The variable lines_to_search is not well defined.')
        exit()
    
    try:
        auxpar = dict()
        
        for key in parameter_functions.keys():
            auxpar[eval(key)] = eval(parameter_functions[key])
        parameter_functions = auxpar
    except:
        parameter_functions = dict()
    
    file = open('Fixed parameter values.txt','w')
    for key in parameters.keys():
        if key not in parameters_on_axes:
            if key not in parameter_functions.keys():
                file.write(latex(key) + ',' + latex(parameters[key]) + '\n')
            else:
                file.write(latex(key) + ',' + latex(parameter_functions[key]) + '\n')
            
    file.close()
        
    file=open('Parameters on axes.txt','w')
    file.write(latex(parameters_on_axes[0]) + ',' + latex(intervalx[0]) + ',' + latex(intervalx[1]) + '\n')
    file.write(latex(parameters_on_axes[1]) + ',' + latex(intervaly[0]) + ',' + latex(intervaly[1]) + '\n')
    file.close()

    try:
        if len(names_of_parameters)==0:
            for parnum in range(2):
                names_of_parameters[parnum] = latex(parameters_on_axes[parnum])
    except:
        names_of_parameters = parameters_on_axes
        for parnum in range(2):
            names_of_parameters[parnum] = latex(names_of_parameters[parnum])
        
    file = open('Actual names of parameters.txt','w')
    file.write(names_of_parameters[0] + ',' + names_of_parameters[1])
    file.close()
    
    if crosscoef=='y':
        file = open('Crosspar.txt','w')
        file.write(latex(crosspar))
        file.close()
    
    if not os.path.isfile('Plotter.nb'):
        shutil.copyfile(os.path.dirname(os.path.realpath(__file__)) + '\\Plotter.nb', 'Plotter.nb')
    
    print('The variables to plot the bifurcation diagram in Mathematica were correctly saved')