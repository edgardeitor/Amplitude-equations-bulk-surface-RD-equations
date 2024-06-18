class Vector:
    '''This class defines a vector with two components. A dummy version of it and
    its actual coordinates'''
    def __init__(self,name):
        self.dummy = []
        self.actualcoord = []
        for varnum in range(nvar):
            self.dummy.append(symbols(name + '_' + str(varnum)))
            self.actualcoord.append(symbols(name + '_' + str(varnum)))
        self.dummy = Matrix(self.dummy)
            
class matrix:
    '''This class defines a sympy Matrix with two components. A dummy version of it and
    its actual coordinates. Here you can provide the actual matrix to be saves into the
    latter one'''
    def __init__(self, name, ncol=nvar, refmatrix=zeros(nvar)):
        self.dummy = zeros(ncol)
        for row in range(ncol):
            for col in range(ncol):
                if refmatrix[row, col]!=0:
                    self.dummy[row, col] = symbols(name + '_' + str(row + 1) + '^' + str(col + 1))
        if refmatrix!=zeros(nvar):
            self.actualcoord = refmatrix
        else:
            self.actualcoord = self.dummy

def isfloat(value):
    '''This function tells you whether a string is a float number or not'''
    try:
        float(eval(value))
        return True
    except:
        return False

def i_n(ell, arg, derivative = False):
    '''This function returns a symbolic expression of the modified spherical Bessel function
    of the first kind i_l(x)'''
    z = symbols('z')
    gbackward = []
    auxiliarybackward = - Pow(z, -2)
    gbackward.append(Pow(z, -1))
    gforward = []
    gforward.append(Pow(z, -1))
    for counter in range(max(ell + 1, - ell)):
        if counter==0:
            gbackward.append(Add(auxiliarybackward, Mul(- 2*counter + 1, Pow(z, -1), gbackward[- 1])))
        else:
            gbackward.append(Add(gbackward[- 2],
                                 Mul(- 2*counter + 1, Pow(z, -1), gbackward[- 1])))
    if max(ell, - ell - 1)>0:
        gforward.append(- Pow(z, -2))
    if max(ell, - ell - 2)>0:
        for counter in range(1, max(ell, - ell - 1)):
            gforward.append(Add(gforward[- 2],
                                Mul(- 1, 2*counter + 1, Pow(z, - 1), gforward[- 1])))
    if ell>=0:
        modified_spherical_bessel = Add(Mul(gforward[-1], sinh(z)), Mul(gbackward[-1], cosh(z)))
    else:
        modified_spherical_bessel = Add(Mul(gbackward[-1], sinh(z)), Mul(gforward[-1], cosh(z)))
    if derivative==False:
        return simplify(modified_spherical_bessel.subs(z, arg))
    else:
        return simplify(diff(modified_spherical_bessel, z).subs(z, arg))

def legendre(ell, m, arg):
    '''This function returns a symbolic expression of the normalized associated Legendre polynomial
    of the first kind P_l^m(x)'''
    if abs(m)>ell:
        return 0
    x = symbols('x')
    leg = Pow(Add(Pow(x, 2), -1), ell)
    for dernum in range(ell + m):
        leg = diff(leg, x)
    leg = Mul(Pow(Add(1, Mul(-1, Pow(x, 2))), Mul(m, Pow(2, - 1))), leg)
    const = Pow(sqrt(integrate(Pow(leg, 2), (x, - 1, 1))), - 1)
    return (simplify(leg.subs(x, arg)), const)

def crossorderapplied(u):
    '''This function is a coded version of F_1'''
    SS = zeros(nvar, 1)
    for counter1 in range(nvar):
        SS = Add(SS, Mul(u.dummy[counter1], firstorderderivatives[counter1]))
    return SS
    
def secondorderapplied(u, v):
    '''This function is a coded version of F_2'''
    DS = zeros(nvar,1)
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            DS = Add(DS, Mul(u.dummy[counter1], v.dummy[counter2],
                             secondorderderivatives[counter1][counter2]))
    return Mul(Pow(2, -1), DS)
    
def thirdorderapplied(u, v, w):
    '''This function is a coded version of F_3'''
    TS = zeros(nvar, 1)
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            for counter3 in range(nvar):
                TS = Add(TS, Mul(u.dummy[counter1], v.dummy[counter2], w.dummy[counter3],
                                 thirdorderderivatives[counter1][counter2][counter3]))
    return Mul(Pow(6, -1), TS)

def dummyvareval(vector, negativeRHS, coefmat):
    '''This function evaluates all the dummy variables used to solve a linear
    system in a simple way'''
    for row in range(nvar):
        vector = vector.subs(negativeRHS.dummy[row], negativeRHS.actualcoord[row])
        for col in range(nvar):
            vector = vector.subs(coefmat.dummy[row, col], coefmat.actualcoord[row, col])
    return vector

def linearsolver(vector,negativeRHS,coefmat):
    '''This function solves a linear system with dummy variables and then uses
    the previous function to evaluate the dummy variables'''
    vector.actualcoord = linsolve(Add(Mul(coefmat.dummy, vector.dummy), negativeRHS.dummy), list(vector.dummy))
    vector.actualcoord = transpose(Matrix(list(vector.actualcoord)))
    vector.actualcoord = dummyvareval(vector.actualcoord, negativeRHS, coefmat)
    return vector.actualcoord

def critical_linearsolver(vector, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols):
    '''This function finds the solution to the equations that determine a solvability condition'''
    vector.dummy[criticalcol] = 0
    vector.actualcoord[criticalcol] = 0
    
    auxiliaryterm, = linsolve(Add(Mul(coefsubmatrix.dummy,Matrix(vector.actualcoord).extract(submatrixcols,[0])),
                                  negativeRHS.dummy.extract(submatrixrows,[0])),
                              list(Matrix(vector.actualcoord).extract(submatrixcols,[0])))
        
    vector.actualcoord[0:criticalcol] = auxiliaryterm[0:criticalcol]
    vector.actualcoord[criticalcol+1:nvar] = auxiliaryterm[criticalcol:nvar-1]
    
    vector.actualcoord = Matrix(vector.actualcoord)
    
    for row in range(nvar):
        vector.actualcoord = vector.actualcoord.subs(negativeRHS.dummy[row],negativeRHS.actualcoord[row])
        for col in range(nvar-1):
            if row<nvar-1:
                vector.actualcoord = vector.actualcoord.subs(coefsubmatrix.dummy[row,col],
                                                             coefsubmatrix.actualcoord[row,col])
                
    return vector.actualcoord

def add_2_up_to_m(ell, m):
    '''This function returns the set of pairs of numbers between -l and l that add up to m'''
    set_check = []
    non_repeated_list = []
    for q1 in range(-ell, ell + 1):
        q2 = m - q1
        if -ell <= q2 <= ell and {q1, q2} not in set_check:
            set_check.append({q1, q2})
            non_repeated_list.append([q1, q2])
    return non_repeated_list

def add_3_up_to_m(ell, m):
    '''This function returns the set of trios of numbers between -l and l that add up to m'''
    set_check = []
    non_repeated_list = []
    for q1 in range(-ell, ell + 1):
        for q2 in range(-ell, ell + 1):
            q3 = m - q1 - q2
            if -ell <= q3 <=ell and {q1, q2, q3} not in set_check:
                set_check.append({q1, q2, q3})
                non_repeated_list.append([q1, q2, q3])
    return non_repeated_list

def I_l():
    '''This function returns the value of the integral I_l in the denominator of the thir-order coefficients'''
    integral = 0
    for varnum in range(nvar):
        integral = Add(integral, Mul(Pow(Sp[lNF][varnum, varnum], 2), phiNF.dummy[varnum],
                                     psiNF.dummy[varnum], Pow(omegaNF.dummy[varnum, varnum], - 3),
                                     integrate(Mul(Pow(ilomegar.actualcoord[varnum, varnum].
                                                       subs(rNF, Mul(sNF, Pow(omegaNF.dummy[varnum, varnum], -1))), 2),
                                                   Pow(sNF, 2)), (sNF, 0, Mul(omegaNF.dummy[varnum, varnum], R)))))
        
    integral = Add(integral, Mul(Pow(R, 2), phiNF.dummy.dot(psiNF.dummy)))
    
    return integral

def d_p(lNF, mNF, p, q_1, q_2):
    '''This function computes d_{p, q_1, q_2}, which is the a coefficient in the expansion of the products
    of Legendre polynomials'''
    pol1, cons1 = legendre(lNF, q_1, x)
    pol2, cons2 = legendre(lNF, q_2, x)
    pol3, cons3 = legendre(p, q_1 + q_2, x)
    dp = Mul(cons1, cons2, cons3, integrate(Mul(pol1, pol2, pol3), (x, - 1, 1)))
    return dp

def first_term_C3(lNF, mNF, q_1, q_2, q_3):
    '''This function returns the first integral on the definition of the third-order coefficients'''
    summation = 0
    pol1, cons1 = legendre(lNF, q_3, x)
    pol2, cons2 = legendre(lNF, mNF, x)
    for p in range(abs(q_1 + q_2), 2*lNF + 1):
        pol3, cons3 = legendre(p, q_1 + q_2, x)
        summation = Add(summation, Mul(DS_phiup2[p - abs(q_1 + q_2)],
                                       cons1, cons2, cons3, integrate(Mul(pol1, pol2, pol3), (x, -1, 1))))
    
    return summation
        
def fourth_integral(lNF, mNF, q_1, q_2, q_3):
    '''This function returns the integral of the product of four associated Legendre polynomials of the
    first kind P_l^{q_1}(x), P_l^{q_2}(x), P_l^{q_3}(x) and P_l^{m}(x)'''
    pol1, cons1 = legendre(lNF, q_1, x)
    pol2, cons2 = legendre(lNF, q_2, x)
    pol3, cons3 = legendre(lNF, q_3, x)
    pol4, cons4 = legendre(lNF, mNF, x)
    integral = Mul(cons1, cons2, cons3, cons4, integrate(Mul(pol1, pol2, pol3, pol4), (x, -1, 1)))
    
    return integral

def non_repeated_element(lista):
    '''This function returns the non-repeated element in a list with three elements with only two different
    elements'''
    if len({lista[0], lista[1]})==1:
        return 2
    elif len({lista[0], lista[2]})==1:
        return 1
    else:
        return 0
    
def evaluation_dict(vector, type = ''):
    '''This function creates a dictionary used to evaluate dummy variables easily'''
    auxvec = Vector('dummy')
    if type=='matrix':
        auxvec.dummy = Matrix([vector.dummy[varnum, varnum] for varnum in range(nvar)])
        auxvec.actualcoord = Matrix([vector.actualcoord[varnum, varnum] for varnum in range(nvar)])
        actual_dict = dict(zip(auxvec.dummy, auxvec.actualcoord))
    else:
        actual_dict = dict(zip(vector.dummy, vector.actualcoord))
    return actual_dict