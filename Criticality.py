from IPython import get_ipython
get_ipython().magic('reset -sf')

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

try:
    exec(open(os.path.dirname(os.path.realpath(__file__)) + '\\functions.py').read())
except:
    print('File functions.py is not in the same folder as the script you are running')
    exit()

file = open('List of Variables.txt','w')

newvar = [None]*2*nvar
for varnum in range(nvar):
    try:
        exec(var[varnum].upper() + ' = symbols(var[varnum].upper(), real=True)')
        exec(var[varnum].lower() + ' = symbols(var[varnum].lower(), real=True)')
        newvar[varnum]=eval(var[varnum].upper())
        newvar[nvar+varnum]=eval(var[varnum].lower())
    except:
        print('The script could not define ' + var[varnum] + ' as a variable')
        exit()
var = Matrix(newvar)

for varnum in range(len(var)):
    file.write(latex(var[varnum]) + '\n')
file.close()
    
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
            exec(key + ' = symbols(key, real=True)')
            newparameters[eval(key)]=parameters[key]
        except:
            print('The script could not define your variable ' + key + ' as a variable')
            exit()
    parameters = newparameters
    
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

try:
    kinetics
except:
    print('Kinetics were not provided.')
    exit()

for functionnumber in range(2*nvar):
    try:
        kinetics[functionnumber] = eval(str(kinetics[functionnumber]).replace('^','**'))
    except:
        print('The expression ' + kinetics[functionnumber] + 'is not a function of the parameters of your system')

kinetics = Matrix(kinetics)

file=open('Kinetics.txt','w')
for functionnumber in range(2*nvar):
    file.write(latex(kinetics[functionnumber]) + '\n')
file.close()

jacobianmat = kinetics.jacobian(var)

if np.count_nonzero(jacobianmat[0:nvar,0:nvar] - np.diag(np.diagonal(jacobianmat[0:nvar, 0:nvar])))>0:
    print('A is not diagonal')
    exit()
if np.count_nonzero(jacobianmat[nvar:2*nvar,0:nvar] - np.diag(np.diagonal(jacobianmat[nvar:2*nvar, 0:nvar])))>0:
    print('K is not diagonal')
    exit()
if np.count_nonzero(diffmatrix[0:nvar,0:nvar] - np.diag(np.diagonal(diffmatrix[0:nvar, 0:nvar])))>0:
    print('The diffusion matrix in the bulk is not diagonal')
    exit()

# Equilibria

r = symbols('r')

B = - jacobianmat[0:nvar, 0:nvar]

K = jacobianmat[nvar:2*nvar, 0:nvar]

omega = zeros(nvar)
omegaeval = zeros(nvar)
for varnum in range(nvar):
    omega[varnum, varnum] = symbols('omega_' + str(varnum+1))
    omegaeval[varnum, varnum] = sqrt(Mul(B[varnum, varnum], Pow(diffmatrix[varnum, varnum], -1)))
    
i0omegar = zeros(nvar)
i0primeomegaR = zeros(nvar)

for varnum in range(nvar):
    i0omegar[varnum, varnum] = i_n(0,Mul(omega[varnum, varnum], r))
    i0primeomegaR[varnum, varnum] = i_n(0,Mul(omega[varnum, varnum], R), derivative=True)
    
i0omegaR = i0omegar.subs(r, R)

diagonal0 = matrix('diag', nvar, Add(Mul(Matrix(diffmatrix[0:nvar, 0:nvar]), i0primeomegaR, omega), Mul(K, i0omegaR)))

S0 = matrix('S0', nvar, Mul(diagonal0.dummy.inv(), K))

Ustar = Mul(i0omegar, S0.dummy, Matrix(var[nvar:2*nvar]))

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

# Turing conditions and Normal form
                    
negativeRHS = Vector('negativeRHS')

firstorderderivatives = list()
secondorderderivatives = list()
thirdorderderivatives = list()
for counter1 in range(nvar, 2*nvar):
    firstorderderivatives.append(diff(Add(Matrix(kinetics[nvar:2*nvar]), Mul(K, Matrix(var[nvar:2*nvar]))), var[counter1]))
    secondorderderivatives.append(list())
    thirdorderderivatives.append(list())
    for counter2 in range(nvar, 2*nvar):
        secondorderderivatives[counter1-nvar].append(diff(firstorderderivatives[counter1-nvar], var[counter2]))
        thirdorderderivatives[counter1-nvar].append(list())
        for counter3 in range(nvar, 2*nvar):
            thirdorderderivatives[counter1-nvar][counter2-nvar].append(diff(
                secondorderderivatives[counter1-nvar][counter2-nvar], var[counter3]))

# firstorderderivatives does not contain the elements with K

x = symbols('x')

file = open('Values of l.txt','w')
for counter in range(len(lvalues)):
    file.write(str(lvalues[counter]) + '\n')
file.close()

s = symbols('s')

for lNF in lvalues:
    phiNF = Vector('phi^NF')
    psiNF = Vector('psi^NF')
    
    if not os.path.isdir('l=' + str(lNF)):
        os.mkdir('l=' + str(lNF))
    os.chdir(os.getcwd() + '\\l=' + str(lNF))
    
    ilomegar = zeros(nvar)
    ilprimeomegaR = zeros(nvar)
    
    i_nuneval = i_n(lNF, r)
    i_nprimeuneval = i_n(lNF, r, derivative = True)

    for varnum in range(nvar):
        ilomegar[varnum, varnum] = i_nuneval.subs(r, Mul(omega[varnum, varnum], r))
        ilprimeomegaR[varnum, varnum] = i_nprimeuneval.subs(r, Mul(omega[varnum, varnum], R))
        
    ilomegar = matrix('ilomegar', nvar, ilomegar)
    ilprimeomegaR = matrix('ilprimeomegaR', nvar, ilprimeomegaR)
    
    ilomegaR = matrix('ilomegaR', nvar, ilomegar.actualcoord.subs(r, R))
    
    diagonall = matrix('diagl', nvar, Add(Mul(Matrix(diffmatrix[0:nvar, 0:nvar]), ilprimeomegaR.dummy, omega),
                                          Mul(K, ilomegaR.dummy)))
    
    Sl = Mul(diagonall.dummy.inv(), K)
    
    criticalmat = Add(Mul(K, ilomegaR.dummy, Sl), jacobianmat[nvar:2*nvar, nvar:2*nvar], \
                      Mul(-1, lNF, lNF + 1, Pow(R, -2), diffmatrix[nvar:2*nvar, nvar:2*nvar]))
        
    # The term -K is hiddein inside jacobianmat
    
    coefmatl = matrix('coefmatl', nvar, criticalmat)
    
    criticalmatdet = coefmatl.dummy.det()
    
    for row in range(nvar):
        for col in range(nvar):
            criticalmatdet = criticalmatdet.subs(coefmatl.dummy[row,col], coefmatl.actualcoord[row,col])
    for varnum in range(nvar):
        criticalmatdet = criticalmatdet.subs(diagonall.dummy[varnum, varnum], diagonall.actualcoord[varnum, varnum])
    for varnum in range(nvar):
        criticalmatdet = criticalmatdet.subs(ilomegaR.dummy[varnum, varnum], ilomegaR.actualcoord[varnum, varnum])
        criticalmatdet = criticalmatdet.subs(ilprimeomegaR.dummy[varnum, varnum],
                                             ilprimeomegaR.actualcoord[varnum, varnum])
    for varnum in range(nvar):
        criticalmatdet = criticalmatdet.subs(omega[varnum, varnum], omegaeval[varnum, varnum])
    
    determinanteval = criticalmatdet
    if eq!=[]:
        for varnum in range(nvar):
            determinanteval = determinanteval.subs(var[nvar + varnum], eq[nvar + varnum])
        for key in parameters:
            determinanteval = determinanteval.subs(key, parameters[key])
        tol=1e-6
        if abs(N(determinanteval))<tol:
            print('There is a Turing bifurcation when l=' + str(lNF) + '.')
        
    file = open('Determinant l=' + str(lNF) + '.txt','w')
    file.write(latex(criticalmatdet))
    file.close()
    
    if thirdcoef=='y' or crosscoef=='y':
        getout=0
        
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
                        submatrixeval = submatrixeval.subs(diagonall.dummy[varnum, varnum],
                                                           diagonall.actualcoord[varnum, varnum])
                    for varnum in range(nvar):
                        submatrixeval = submatrixeval.subs(ilomegaR.dummy[varnum, varnum],
                                                           ilomegaR.actualcoord[varnum, varnum])
                        submatrixeval = submatrixeval.subs(ilprimeomegaR.dummy[varnum, varnum],
                                                           ilprimeomegaR.actualcoord[varnum, varnum])
                    for varnum in range(nvar):
                        submatrixeval = submatrixeval.subs(omega[varnum,varnum], omegaeval[varnum,varnum])
                    for key in parameters:
                        submatrixeval = submatrixeval.subs(key, parameters[key])
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
            submatrixrows = list(range(1,nvar))
            submatrixcols = list(range(1,nvar))
            invertiblesubmatrix = criticalmat.extract(submatrixrows, submatrixcols)
            phiNF.actualcoord[0] = 1
            criticalrow = 0
            criticalcol = 0
            
        coefsubmatrix = matrix('dummysubmatrix', nvar-1, invertiblesubmatrix)
        
        auxiliaryterm, = linsolve(Add(Mul(coefsubmatrix.dummy, Matrix(phiNF.actualcoord).extract(submatrixcols, [0])), \
                                      coefmatl.dummy.extract(submatrixrows, [criticalcol])),
                                  list(Matrix(phiNF.actualcoord).extract(submatrixcols, [0])))
            
        phiNF.actualcoord[0:criticalcol] = auxiliaryterm[0:criticalcol]
        phiNF.actualcoord[criticalcol+1:nvar] = auxiliaryterm[criticalcol:nvar-1]
        
        phiNF.actualcoord = Matrix(phiNF.actualcoord)
        
        for row in range(nvar):
            for col in range(nvar):
                phiNF.actualcoord = phiNF.actualcoord.subs(coefmatl.dummy[row, col], coefmatl.actualcoord[row, col])
                if row<nvar-1 and col<nvar-1:
                    phiNF.actualcoord = phiNF.actualcoord.subs(coefsubmatrix.dummy[row, col],
                                                               coefsubmatrix.actualcoord[row, col])
                    
        if phiunit=='y':
            phiNF.actualcoord = Mul(Pow(phiNF.actualcoord.dot(phiNF.actualcoord), -1), phiNF.actualcoord)
        
        print('phiNF ready')
        
        Tdummycoefsubmatrix = transpose(coefsubmatrix.dummy)
        Tdummycoefmatl = transpose(coefmatl.dummy)
        
        psiNF.actualcoord[criticalrow] = 1
        
        auxiliaryterm, = linsolve(Add(Mul(Tdummycoefsubmatrix, Matrix(psiNF.actualcoord).extract(submatrixrows, [0])), \
                                      Tdummycoefmatl.extract(submatrixcols, [criticalrow])),
                                  list(Matrix(psiNF.actualcoord).extract(submatrixrows, [0])))
            
        psiNF.actualcoord[0:criticalrow] = auxiliaryterm[0:criticalrow]
        psiNF.actualcoord[criticalrow+1:nvar] = auxiliaryterm[criticalrow:nvar-1]
        
        psiNF.actualcoord = Matrix(psiNF.actualcoord)
        
        for row in range(nvar):
            for col in range(nvar):
                psiNF.actualcoord = psiNF.actualcoord.subs(coefmatl.dummy[row, col], coefmatl.actualcoord[row, col])
                if row<nvar-1 and col<nvar-1:
                    psiNF.actualcoord = psiNF.actualcoord.subs(coefsubmatrix.dummy[row, col],
                                                               coefsubmatrix.actualcoord[row, col])
        
        print('psiNF ready')
        
        integral = I_l()
        
        print('Hellish integral done')
    
    if thirdcoef=='y':
        DS_phiphi = secondorderapplied(phiNF.dummy, phiNF.dummy)
        
        if lNF%2==1:
            TS_phiphiphi = thirdorderapplied(phiNF.dummy, phiNF.dummy, phiNF.dummy)
            
        for mNF in range(0, lNF + 1):
            if lNF%2==0:
                list_of_q = add_2_up_to_m(lNF, mNF)
                
                for part_list in list_of_q:
                    q_1 = part_list[0]
                    q_2 = part_list[1]
                    
                    dl = d_p(lNF, mNF, lNF, q_1, q_2)
                    
                    C2 = Mul(Pow(R, 2), Pow(integral, -1), dl, psiNF.dummy.dot(DS_phiphi))
                    
                    if q_1!=q_2:
                        C2 = Mul(2, C2)
                    
                    for varnum in range(nvar):
                        C2 = C2.subs(phiNF.dummy[varnum], phiNF.actualcoord[varnum])
                        C2 = C2.subs(psiNF.dummy[varnum], psiNF.actualcoord[varnum])
                    for varnum in range(nvar):
                        C2 = C2.subs(diagonall.dummy[varnum, varnum], diagonall.actualcoord[varnum, varnum])
                    for varnum in range(nvar):
                        C2 = C2.subs(ilomegar.dummy[varnum, varnum], ilomegar.actualcoord[varnum, varnum])
                        C2 = C2.subs(ilomegaR.dummy[varnum, varnum], ilomegaR.actualcoord[varnum, varnum])
                        C2 = C2.subs(ilprimeomegaR.dummy[varnum, varnum], ilprimeomegaR.actualcoord[varnum, varnum])
                    for varnum in range(nvar):
                        C2 = C2.subs(omega[varnum, varnum], omegaeval[varnum, varnum])
                
                    file = open('C2 q1=' + str(q_1) + ', q2=' + str(q_2) + ', m=' + str(mNF) + '.txt','w')
                    file.write(latex(C2))
                    file.close()
                
                    print('Second order coefficient ready')
                    
            else:
                Nlm = Mul(2, math.factorial(lNF + abs(mNF)), Pow(2*lNF + 1, -1),
                          Pow(math.factorial(lNF - abs(mNF)), -1))
                
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
                    
                        for p in range(abs(q_1 + q_2), 2*lNF + 1):
                            ipomegaR = zeros(nvar)
                            ipprimeomegaR = zeros(nvar)
                
                            ipuneval = i_n(p, r)
                            ipprimeuneval = i_n(p, r, derivative = True)
                
                            for varnum in range(nvar):
                                ipomegaR[varnum, varnum] = ipuneval.subs(r, Mul(omega[varnum, varnum], R))
                                ipprimeomegaR[varnum, varnum] = ipprimeuneval.subs(r, Mul(omega[varnum, varnum], R))
                                
                            ipomegaR = matrix('ipomegaR', nvar, ipomegaR)
                            ipprimeomegaR = matrix('ipprimeomegaR', nvar, ipprimeomegaR)
                            
                            diagonalp = matrix('diagp', nvar, Add(Mul(Matrix(diffmatrix[0:nvar,0:nvar]),
                                                                      ipprimeomegaR.dummy, omega),
                                                                  Mul(K, ipomegaR.dummy)))
                            
                            Sp = Mul(diagonalp.dummy.inv(), K)
                            
                            criticalmatp = matrix('coefmatp', nvar, Add(Mul(K, ipomegaR.dummy, Sp),
                                                                        Matrix(jacobianmat[nvar:2*nvar, nvar:2*nvar]),
                                                                        Mul(-1, p, p + 1, Pow(R, -2)),
                                                                            Matrix(diffmatrix[nvar:2*nvar,
                                                                                              nvar:2*nvar]))))
                            
                            dp = d_p(lNF, mNF, p, q_1, q_2)
                            
                            up2.append(Vector('up2^' + str(p)))
                            
                            if dp==0:
                                up2[-1].dummy = Matrix([0]*nvar)
                                up2[-1].actualcoord = Matrix([0]*nvar)
                            else:
                                negativeRHS.actualcoord = Mul(dp, DS_phiphi)
                                up2[-1].actualcoord = linearsolver(up2[-1], negativeRHS, criticalmatp)
                                
                                for varnum in range(nvar):
                                    up2[-1].actualcoord = up2[-1].actualcoord.subs(diagonalp.dummy[varnum, varnum],
                                                                                   diagonalp.actualcoord[varnum, varnum])
                                for varnum in range(nvar):
                                    up2[-1].actualcoord = up2[-1].actualcoord.subs(ipomegaR.dummy[varnum, varnum],
                                                                                   ipomegaR.actualcoord[varnum, varnum])
                                    up2[-1].actualcoord = up2[-1].actualcoord.subs(ipprimeomegaR.dummy[varnum, varnum],
                                                                                   ipprimeomegaR.actualcoord[varnum, varnum])
                        
                            DS_phiup2.append(secondorderapplied(phiNF.dummy, up2[p - abs(q_1 + q_2)].dummy))
                        
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
                        
                        for p in range(abs(q_1 + q_2), 2*lNF + 1):
                            for varnum in range(nvar):
                                summation = summation.subs(up2[p - abs(q_1 + q_2)].dummy[varnum],
                                                           up2[p - abs(q_1 + q_2)].actualcoord[varnum])
                        
                        print('Hard part ready')
                        
                        if q_1!=q_2:
                            summation = Mul(2, summation)
                        
                        bigsum = Add(bigsum, summation)
                        
                    legendreintegral = fourth_integral(lNF, mNF, q_1, q_2, q_3)
                    
                    second_summand = Mul(TS_phiphiphi, legendreintegral)
                    
                    if len({q_1, q_2, q_3})==3:
                        second_summand = Mul(6, second_summand)
                    elif len({q_1, q_2, q_3})==2:
                        second_summand = Mul(3, second_summand)
                    
                    C3 = Mul(Pow(R, 2), Pow(integral, -1), Pow(Nlm, -1),
                             psiNF.dummy.dot(Add(Mul(2, bigsum), second_summand)))
                    
                    for varnum in range(nvar):
                        C3 = C3.subs(phiNF.dummy[varnum], phiNF.actualcoord[varnum])
                        C3 = C3.subs(psiNF.dummy[varnum], psiNF.actualcoord[varnum])
                    for varnum in range(nvar):
                        C3 = C3.subs(diagonall.dummy[varnum, varnum], diagonall.actualcoord[varnum, varnum])
                    for varnum in range(nvar):
                        C3 = C3.subs(ilomegar.dummy[varnum, varnum], ilomegar.actualcoord[varnum, varnum])
                        C3 = C3.subs(ilomegaR.dummy[varnum, varnum], ilomegaR.actualcoord[varnum, varnum])
                        C3 = C3.subs(ilprimeomegaR.dummy[varnum, varnum],
                                     ilprimeomegaR.actualcoord[varnum, varnum])
                    for varnum in range(nvar):
                        C3 = C3.subs(omega[varnum,varnum], omegaeval[varnum,varnum])
                        
                    file = open('C3 q1=' + str(q_1) + ', q2=' + str(q_2) + ', q3=' + str(q_3) + ', m=' +
                                str(mNF) + '.txt','w')
                    file.write(latex(C3))
                    file.close()
                            
                    print('Third order coefficient ready')
    
    if crosscoef=='y':
        try:
            crosspar = eval(crosspar)
            for varnum in range(len(equilibrium)):
                equilibrium[varnum] = eval(equilibrium[varnum])
        except:
            if lNF==lvalues[0]:
                print('The Cross parameter or equilibrium were not provided in the right way')
        
        # kineticsevaluated = Matrix(kinetics[nvar:2*nvar])
        # for varnum in range(nvar):
        #     kineticsevaluated = kineticsevaluated.subs(var[varnum],
        #                                                Mul(i0omegaR[varnum,varnum], S0.dummy[varnum,varnum],
        #                                                    var[nvar + varnum]))
        # for varnum in range(nvar):
        #     kineticsevaluated = kineticsevaluated.subs(omega[varnum, varnum], omegaeval[varnum, varnum])
        # for key in parameters:
        #     if key!=crosspar:
        #         kineticsevaluated = kineticsevaluated.subs(key, parameters[key])
        # equilibrium = solve(kineticsevaluated, var[nvar:2*nvar])[0]
                    
        SS_phi = crossorderapplied(phiNF.dummy, crosspar, equilibrium)
        
        firsthellishterm = 0
        
        if diff(B, crosspar)!=zeros(nvar):
            for varnum in range(nvar):
                firsthellishterm = Add(firsthellishterm,
                                        Mul(-diff(B[varnum, varnum], crosspar), Pow(Sl[varnum, varnum], 2),
                                            phiNF.dummy[varnum], psiNF.dummy[varnum], Pow(omega[varnum, varnum], -3),
                                            integrate(Mul(Pow(ilomegar.actualcoord[varnum, varnum].
                                                              subs(r, Mul(s, Pow(omega[varnum, varnum], -1))), 2),
                                                          Pow(s, 2)), (s, 0, omega[varnum, varnum] * R))))
        
        firsthellishterm = Add(firsthellishterm, Mul(Pow(R, 2), psiNF.dummy.dot(SS_phi)))
        
        if diff(K, crosspar)!=zeros(nvar):
            firsthellishterm = Add(firsthellishterm, Mul(Pow(R, 2),
                                                          psiNF.dummy.dot(Mul(-diff(K, crosspar),
                                                                              Add(phiNF.dummy,
                                                                                  Mul(-1, ilomegaR.dummy,
                                                                                      Sl, phiNF.dummy))))))
        
        secondhellishterm = 0
        
        if diff(Matrix(diffmatrix[0:nvar, 0:nvar]), crosspar)!=zeros(nvar):
            for varnum in range(nvar):
                secondhellishterm = Add(secondhellishterm,
                                        Mul(diff(diffmatrix[varnum, varnum], crosspar), Pow(Sl[varnum, varnum], 2),
                                            phiNF.dummy[varnum], psiNF.dummy[varnum], Pow(omega[varnum, varnum], -3),
                                            integrate(Mul(Pow(ilomegar.actualcoord[varnum, varnum].
                                                              subs(r, Mul(s, Pow(omega[varnum, varnum], -1))), 2),
                                                          Pow(s, 2)), (s, 0, omega[varnum, varnum] * R))))
                
        if diff(Matrix(diffmatrix[nvar:2*nvar, nvar:2*nvar]), crosspar)!=zeros(nvar):
            secondhellishterm = Add(secondhellishterm,
                                    Mul(Pow(R,2),
                                        psiNF.dummy.dot(Mul(diff(diffmatrix[nvar:2*nvar, nvar:2*nvar], crosspar),
                                                            Matrix(phiNF.dummy)))))
            
        C11 = Mul(Pow(integral, -1), Add(firsthellishterm, Mul(-lNF, lNF + 1, Pow(R, -2), secondhellishterm)))
        
        for varnum in range(nvar):
            C11 = C11.subs(phiNF.dummy[varnum], phiNF.actualcoord[varnum])
            C11 = C11.subs(psiNF.dummy[varnum], psiNF.actualcoord[varnum])
        for varnum in range(nvar):
            C11 = C11.subs(diagonall.dummy[varnum, varnum], diagonall.actualcoord[varnum, varnum])    
        for varnum in range(nvar):
            C11 = C11.subs(ilomegar.dummy[varnum, varnum], ilomegar.actualcoord[varnum, varnum])
            C11 = C11.subs(ilomegaR.dummy[varnum, varnum], ilomegaR.actualcoord[varnum, varnum])
            C11 = C11.subs(ilprimeomegaR.dummy[varnum, varnum],
                         ilprimeomegaR.actualcoord[varnum, varnum])
        for varnum in range(nvar):
            C11 = C11.subs(omega[varnum, varnum], omegaeval[varnum, varnum])
        
        file = open('Cross-order coefficient l=' + str(lNF) + '.txt','w')        
        file.write(latex(C11))
        file.close()
    
    os.chdir(os.path.dirname(os.getcwd()))
    
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